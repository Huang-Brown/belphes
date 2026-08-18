// build_templates.cpp -- step 2 of the belphes template pipeline.
/**
 *  Turns CMS NanoAOD into joint DeepJet templates for belphes.
 *
 *  For every (hadronFlavour, |eta| bin, pt bin) cell we store ONE 2-D
 *  histogram of (btagDeepFlavB, btagDeepFlavCvL).  Two axes are enough:
 *  DeepJet's grouped probabilities satisfy
 *
 *      B + C + L = 1,   CvB = C/(C+B),   CvL = C/(C+L)
 *
 *  so given (B, CvL) the rest follows exactly:
 *
 *      C   = CvL * (1 - B)
 *      CvB = C / (C + B)
 *
 *  Storing (B, CvL) rather than (B, CvB) is deliberate -- recovering CvB from
 *  (B, CvL) is numerically stable everywhere, whereas the reverse divides by
 *  (1 - CvB) and loses all precision for the high-B jets that matter most.
 *  Verified on open-data record 67727: max |residual| 4.3e-4 over 8924 jets,
 *  i.e. NanoAOD's float storage precision.
 *
 *  The pt/eta edges, flavour list and score resolution are written into the
 *  output file, so the sampling module reads the binning back from the
 *  templates instead of having it restated in a Delphes card.
 *
 *  Compile:
 *    g++ -O2 -Wall build_templates.cpp $(root-config --cflags --libs) -o build_templates
 *
 *  Run (see build_templates.sh for the usual invocation):
 *    ./build_templates -c binning.conf -o templates.root  INPUT.root [INPUT.root ...]
 */

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

struct Config
{
  std::vector<Double_t> ptBins;
  std::vector<Double_t> absEtaBins;
  std::vector<Int_t> flavors;
  Int_t scoreBinsB = 100;
  Int_t scoreBinsCvL = 100;
  // Optional explicit score-axis edges. When set, these override the uniform
  // bin counts above. CvB = C/(C+B) is nonlinear and most sensitive where B is
  // small, so concentrating bins near 0 buys accuracy for c and light jets
  // without paying for 100+ uniform bins everywhere.
  std::vector<Double_t> scoreEdgesB;
  std::vector<Double_t> scoreEdgesCvL;
  Double_t jetPtMin = 20.0;
  Double_t jetAbsEtaMax = 2.5;
  Long64_t minCellEntries = 500;
};

static bool ReadConfig(const std::string &path, Config &cfg)
{
  std::ifstream in(path);
  if(!in) { std::cerr << "ERROR: cannot open config " << path << "\n"; return false; }

  std::string line;
  while(std::getline(in, line))
  {
    if(auto h = line.find('#'); h != std::string::npos) line.erase(h);
    std::istringstream ls(line);
    std::string key;
    if(!(ls >> key)) continue;

    if(key == "PtBins")             { Double_t v; cfg.ptBins.clear();     while(ls >> v) cfg.ptBins.push_back(v); }
    else if(key == "AbsEtaBins")    { Double_t v; cfg.absEtaBins.clear(); while(ls >> v) cfg.absEtaBins.push_back(v); }
    else if(key == "Flavors")       { Int_t v;    cfg.flavors.clear();    while(ls >> v) cfg.flavors.push_back(v); }
    else if(key == "ScoreBinsB")    ls >> cfg.scoreBinsB;
    else if(key == "ScoreBinsCvL")  ls >> cfg.scoreBinsCvL;
    else if(key == "ScoreEdgesB")   { Double_t v; cfg.scoreEdgesB.clear();   while(ls >> v) cfg.scoreEdgesB.push_back(v); }
    else if(key == "ScoreEdgesCvL") { Double_t v; cfg.scoreEdgesCvL.clear(); while(ls >> v) cfg.scoreEdgesCvL.push_back(v); }
    else if(key == "JetPtMin")      ls >> cfg.jetPtMin;
    else if(key == "JetAbsEtaMax")  ls >> cfg.jetAbsEtaMax;
    else if(key == "MinCellEntries")ls >> cfg.minCellEntries;
    else std::cerr << "WARNING: unknown config key '" << key << "' ignored\n";
  }

  // A silently unsorted edge list would corrupt every lookup downstream.
  auto checkSorted = [](const std::vector<Double_t> &v, const char *name) {
    if(v.size() < 2)                                  { std::cerr << "ERROR: " << name << " needs >= 2 edges\n"; return false; }
    if(!std::is_sorted(v.begin(), v.end()))           { std::cerr << "ERROR: " << name << " is not ascending\n"; return false; }
    if(std::adjacent_find(v.begin(), v.end()) != v.end()) { std::cerr << "ERROR: " << name << " has duplicate edges\n"; return false; }
    return true;
  };
  if(!checkSorted(cfg.ptBins, "PtBins")) return false;
  if(!checkSorted(cfg.absEtaBins, "AbsEtaBins")) return false;
  if(cfg.flavors.empty())                     { std::cerr << "ERROR: Flavors is empty\n"; return false; }
  if(cfg.scoreBinsB < 1 || cfg.scoreBinsCvL < 1) { std::cerr << "ERROR: score bin counts must be >= 1\n"; return false; }

  // Score axes must span exactly [0,1]: the discriminants are defined there,
  // and anything outside would be silently dropped into over/underflow, which
  // TH2::GetRandom2 never samples.
  auto checkScore = [&](const std::vector<Double_t> &v, const char *name) {
    if(v.empty()) return true;
    if(!checkSorted(v, name)) return false;
    if(std::fabs(v.front()) > 1e-12 || std::fabs(v.back() - 1.0) > 1e-12)
    { std::cerr << "ERROR: " << name << " must run from 0 to 1\n"; return false; }
    return true;
  };
  if(!checkScore(cfg.scoreEdgesB, "ScoreEdgesB")) return false;
  if(!checkScore(cfg.scoreEdgesCvL, "ScoreEdgesCvL")) return false;
  return true;
}

//------------------------------------------------------------------------------

// Index of the cell containing x, or -1 if x is outside [edges.front(), edges.back()).
static Int_t FindCell(const std::vector<Double_t> &edges, Double_t x)
{
  if(x < edges.front() || x >= edges.back()) return -1;
  return Int_t(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
}

static TString HistName(Int_t flavor, Int_t iEta, Int_t iPt)
{
  return TString::Format("tmpl_f%d_eta%d_pt%d", flavor, iEta, iPt);
}

//------------------------------------------------------------------------------

int main(int argc, char **argv)
{
  std::string configPath = "binning.conf";
  std::string outPath = "templates.root";
  std::vector<std::string> inputs;

  for(int i = 1; i < argc; ++i)
  {
    std::string a = argv[i];
    if((a == "-c" || a == "--config") && i + 1 < argc)      configPath = argv[++i];
    else if((a == "-o" || a == "--output") && i + 1 < argc) outPath = argv[++i];
    else if(a == "-h" || a == "--help")
    {
      std::cout << "usage: " << argv[0] << " [-c binning.conf] [-o templates.root] <nanoaod.root> ...\n";
      return 0;
    }
    else inputs.push_back(a);
  }

  if(inputs.empty()) { std::cerr << "ERROR: no input files given\n"; return 1; }

  Config cfg;
  if(!ReadConfig(configPath, cfg)) return 1;

  const Int_t nEta = Int_t(cfg.absEtaBins.size()) - 1;
  const Int_t nPt = Int_t(cfg.ptBins.size()) - 1;
  const Int_t nFlav = Int_t(cfg.flavors.size());

  std::printf("==> config %s\n", configPath.c_str());
  std::printf("    grid          : %d flavour x %d |eta| x %d pt = %d cells\n",
              nFlav, nEta, nPt, nFlav * nEta * nPt);
  const Int_t nScoreB = cfg.scoreEdgesB.empty() ? cfg.scoreBinsB : Int_t(cfg.scoreEdgesB.size()) - 1;
  const Int_t nScoreCvL = cfg.scoreEdgesCvL.empty() ? cfg.scoreBinsCvL : Int_t(cfg.scoreEdgesCvL.size()) - 1;
  std::printf("    resolution    : %d (B%s) x %d (CvL%s) = %d bins/cell\n",
              nScoreB, cfg.scoreEdgesB.empty() ? "" : ", variable",
              nScoreCvL, cfg.scoreEdgesCvL.empty() ? "" : ", variable",
              nScoreB * nScoreCvL);
  std::printf("    jet selection : pt > %.1f, |eta| < %.2f\n", cfg.jetPtMin, cfg.jetAbsEtaMax);

  // ---- book one joint template per cell ------------------------------------
  // Indexed [flavour][eta][pt]; owned by the output file below.
  std::vector<std::vector<std::vector<TH2F *>>> hist(
    nFlav, std::vector<std::vector<TH2F *>>(nEta, std::vector<TH2F *>(nPt, nullptr)));

  TFile *fout = TFile::Open(outPath.c_str(), "RECREATE");
  if(!fout || fout->IsZombie()) { std::cerr << "ERROR: cannot create " << outPath << "\n"; return 1; }

  for(Int_t f = 0; f < nFlav; ++f)
    for(Int_t e = 0; e < nEta; ++e)
      for(Int_t p = 0; p < nPt; ++p)
      {
        TString name = HistName(cfg.flavors[f], e, p);
        TString title = TString::Format(
          "hadronFlavour %d, %.2f < |#eta| < %.2f, %.0f < p_{T} < %.0f GeV;btagDeepFlavB;btagDeepFlavCvL",
          cfg.flavors[f], cfg.absEtaBins[e], cfg.absEtaBins[e + 1], cfg.ptBins[p], cfg.ptBins[p + 1]);
        // Variable edges when given, uniform otherwise. TH2::GetRandom2
        // honours per-bin widths either way.
        if(!cfg.scoreEdgesB.empty() && !cfg.scoreEdgesCvL.empty())
          hist[f][e][p] = new TH2F(name, title,
                                   Int_t(cfg.scoreEdgesB.size()) - 1, cfg.scoreEdgesB.data(),
                                   Int_t(cfg.scoreEdgesCvL.size()) - 1, cfg.scoreEdgesCvL.data());
        else if(!cfg.scoreEdgesB.empty())
          hist[f][e][p] = new TH2F(name, title,
                                   Int_t(cfg.scoreEdgesB.size()) - 1, cfg.scoreEdgesB.data(),
                                   cfg.scoreBinsCvL, 0.0, 1.0);
        else if(!cfg.scoreEdgesCvL.empty())
          hist[f][e][p] = new TH2F(name, title,
                                   cfg.scoreBinsB, 0.0, 1.0,
                                   Int_t(cfg.scoreEdgesCvL.size()) - 1, cfg.scoreEdgesCvL.data());
        else
          hist[f][e][p] = new TH2F(name, title,
                                   cfg.scoreBinsB, 0.0, 1.0,
                                   cfg.scoreBinsCvL, 0.0, 1.0);
        hist[f][e][p]->SetDirectory(fout);
      }

  // Map hadronFlavour -> index, so the event loop stays a hash lookup.
  std::map<Int_t, Int_t> flavIndex;
  for(Int_t f = 0; f < nFlav; ++f) flavIndex[cfg.flavors[f]] = f;

  // ---- event loop ----------------------------------------------------------
  TChain chain("Events");
  for(const auto &in : inputs)
    if(chain.Add(in.c_str()) == 0) std::cerr << "WARNING: no tree added from " << in << "\n";

  std::printf("==> %d input file(s), %lld events\n", int(inputs.size()), chain.GetEntries());

  TTreeReader reader(&chain);
  TTreeReaderArray<Float_t> jetPt(reader, "Jet_pt");
  TTreeReaderArray<Float_t> jetEta(reader, "Jet_eta");
  TTreeReaderArray<Int_t> jetFlav(reader, "Jet_hadronFlavour");
  TTreeReaderArray<Float_t> jetB(reader, "Jet_btagDeepFlavB");
  TTreeReaderArray<Float_t> jetCvL(reader, "Jet_btagDeepFlavCvL");

  Long64_t nJets = 0, nKept = 0;
  Long64_t nOutOfGrid = 0, nBadFlavor = 0, nBadScore = 0;

  while(reader.Next())
  {
    for(size_t j = 0; j < jetPt.GetSize(); ++j)
    {
      ++nJets;

      const Double_t pt = jetPt[j];
      const Double_t absEta = std::fabs(jetEta[j]);
      if(pt < cfg.jetPtMin || absEta > cfg.jetAbsEtaMax) continue;

      auto itFlav = flavIndex.find(jetFlav[j]);
      if(itFlav == flavIndex.end()) { ++nBadFlavor; continue; }

      // NanoAOD writes -1 when a discriminator's denominator vanishes.
      const Double_t B = jetB[j], CvL = jetCvL[j];
      if(B < 0.0 || B > 1.0 || CvL < 0.0 || CvL > 1.0) { ++nBadScore; continue; }

      const Int_t iEta = FindCell(cfg.absEtaBins, absEta);
      const Int_t iPt = FindCell(cfg.ptBins, pt);
      if(iEta < 0 || iPt < 0) { ++nOutOfGrid; continue; }

      hist[itFlav->second][iEta][iPt]->Fill(B, CvL);
      ++nKept;
    }
  }

  std::printf("==> jets seen %lld, filled %lld\n", nJets, nKept);
  std::printf("    rejected: %lld outside grid, %lld unlisted flavour, %lld invalid score\n",
              nOutOfGrid, nBadFlavor, nBadScore);

  // ---- self-describing metadata -------------------------------------------
  // The module reads the binning from here; nothing restates it in a card.
  fout->cd();

  TVectorD vPt(nPt + 1), vEta(nEta + 1), vFlav(nFlav);
  for(Int_t i = 0; i <= nPt; ++i) vPt[i] = cfg.ptBins[i];
  for(Int_t i = 0; i <= nEta; ++i) vEta[i] = cfg.absEtaBins[i];
  for(Int_t i = 0; i < nFlav; ++i) vFlav[i] = cfg.flavors[i];
  vPt.Write("PtBins");
  vEta.Write("AbsEtaBins");
  vFlav.Write("Flavors");

  TNamed("Basis", "joint (btagDeepFlavB, btagDeepFlavCvL); CvB derived as "
                  "C=CvL*(1-B), CvB=C/(C+B)").Write();
  TNamed(TString("Provenance"), TString::Format(
    "built from %d NanoAOD file(s); first input: %s", int(inputs.size()), inputs[0].c_str())).Write();

  // ---- per-cell report -----------------------------------------------------
  Int_t nEmpty = 0, nThin = 0;
  Double_t minEntries = -1;
  TString minCell;
  std::vector<std::pair<Double_t, TString>> thin;
  for(Int_t f = 0; f < nFlav; ++f)
    for(Int_t e = 0; e < nEta; ++e)
      for(Int_t p = 0; p < nPt; ++p)
      {
        const Double_t n = hist[f][e][p]->GetEntries();
        if(n <= 0) ++nEmpty;
        else if(n < cfg.minCellEntries) { ++nThin; thin.emplace_back(n, hist[f][e][p]->GetName()); }
        if(minEntries < 0 || n < minEntries) { minEntries = n; minCell = hist[f][e][p]->GetName(); }
      }
  std::sort(thin.begin(), thin.end(),
            [](const std::pair<Double_t, TString> &a, const std::pair<Double_t, TString> &b)
            { return a.first < b.first; });

  fout->Write();
  fout->Close();
  delete fout;

  std::printf("==> wrote %s\n", outPath.c_str());
  std::printf("    cells empty            : %d\n", nEmpty);
  std::printf("    cells under %-10lld : %d\n", cfg.minCellEntries, nThin);
  std::printf("    sparsest cell          : %s (%.0f entries)\n", minCell.Data(), minEntries);

  const Int_t nScoreBins = nScoreB * nScoreCvL;

  if(nEmpty > 0)
    std::printf("\nWARNING: %d empty cell(s). TH2::GetRandom2 returns (0,0) for an empty\n"
                "         histogram, so PseudoDeepFlavScore refuses to load a file\n"
                "         containing any. Widen the bins or add input files.\n", nEmpty);

  if(nThin > 0)
  {
    std::printf("\nWARNING: %d cell(s) below MinCellEntries = %lld. These are not empty,\n"
                "         so they will load, but a cell with N entries spread over %d bins\n"
                "         samples as a comb of at most N spikes rather than a distribution.\n",
                nThin, cfg.minCellEntries, nScoreBins);
    std::printf("         Sparsest %d:\n", int(std::min<size_t>(thin.size(), 8)));
    for(size_t i = 0; i < thin.size() && i < 8; ++i)
      std::printf("           %-24s %6.0f entries (%.3f per bin)\n",
                  thin[i].second.Data(), thin[i].first, thin[i].first / nScoreBins);
    std::printf("         Fix by coarsening PtBins/AbsEtaBins in the forward, high-pt\n"
                "         corner, or by lowering the score resolution.\n");
  }

  return 0;
}
