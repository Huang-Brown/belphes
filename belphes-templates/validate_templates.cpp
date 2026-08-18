// validate_templates.cpp -- closure test for the belphes joint templates.
/**
 *  Answers three questions about a template file, per (flavour, |eta|, pt) cell:
 *
 *    1. Closure     -- does sampling the template reproduce the NanoAOD it was
 *                      built from (marginal means of B and CvL)?
 *    2. Correlation -- is rho(B, CvL) preserved?  This is what a joint template
 *                      buys over two independent 1-D templates, so the run also
 *                      samples the marginals independently for contrast.
 *    3. Derivation  -- does CvB, reconstructed as C = CvL*(1-B), CvB = C/(C+B),
 *                      match the CvB actually stored in NanoAOD?
 *
 *  Compile:
 *    g++ -O2 -Wall validate_templates.cpp $(root-config --cflags --libs) -o validate_templates
 *
 *  Run:
 *    ./validate_templates -t templates.root -f 5 INPUT.root [INPUT.root ...]
 */

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

// The one derivation the whole design rests on.  Kept here (and mirrored in the
// Delphes module) so both sides of the pipeline agree by construction.
static inline Double_t DeriveCvB(Double_t B, Double_t CvL)
{
  const Double_t C = CvL * (1.0 - B);
  const Double_t den = C + B;
  return den > 0.0 ? C / den : 0.0;
}

struct Stats
{
  Long64_t n = 0;
  Double_t sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;

  void Add(Double_t x, Double_t y)
  {
    ++n; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
  }
  Double_t MeanX() const { return n ? sx / n : 0; }
  Double_t MeanY() const { return n ? sy / n : 0; }
  Double_t Rho() const
  {
    if(n < 2) return 0;
    const Double_t vx = sxx / n - MeanX() * MeanX();
    const Double_t vy = syy / n - MeanY() * MeanY();
    if(vx <= 0 || vy <= 0) return 0;
    return (sxy / n - MeanX() * MeanY()) / std::sqrt(vx * vy);
  }
};

static Int_t FindCell(const std::vector<Double_t> &edges, Double_t x)
{
  if(x < edges.front() || x >= edges.back()) return -1;
  return Int_t(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
}

//------------------------------------------------------------------------------

int main(int argc, char **argv)
{
  std::string tmplPath = "templates.root";
  Int_t wantFlavor = 5;
  Long64_t nSample = 400000;
  std::vector<std::string> inputs;

  for(int i = 1; i < argc; ++i)
  {
    std::string a = argv[i];
    if((a == "-t" || a == "--templates") && i + 1 < argc) tmplPath = argv[++i];
    else if((a == "-f" || a == "--flavor") && i + 1 < argc) wantFlavor = atoi(argv[++i]);
    else if((a == "-n" || a == "--nsample") && i + 1 < argc) nSample = atoll(argv[++i]);
    else if(a == "-h" || a == "--help")
    {
      std::cout << "usage: " << argv[0]
                << " [-t templates.root] [-f FLAVOUR] [-n NSAMPLE] INPUT.root ...\n";
      return 0;
    }
    else inputs.push_back(a);
  }
  if(inputs.empty()) { std::cerr << "ERROR: no NanoAOD input given\n"; return 1; }

  TFile *ft = TFile::Open(tmplPath.c_str(), "READ");
  if(!ft || ft->IsZombie()) { std::cerr << "ERROR: cannot open " << tmplPath << "\n"; return 1; }

  // Binning comes back out of the template file -- nothing is restated here.
  TVectorD *vPt = dynamic_cast<TVectorD *>(ft->Get("PtBins"));
  TVectorD *vEta = dynamic_cast<TVectorD *>(ft->Get("AbsEtaBins"));
  if(!vPt || !vEta) { std::cerr << "ERROR: binning metadata missing from " << tmplPath << "\n"; return 1; }

  std::vector<Double_t> ptBins(vPt->GetNrows()), etaBins(vEta->GetNrows());
  for(Int_t i = 0; i < vPt->GetNrows(); ++i) ptBins[i] = (*vPt)[i];
  for(Int_t i = 0; i < vEta->GetNrows(); ++i) etaBins[i] = (*vEta)[i];
  const Int_t nPt = Int_t(ptBins.size()) - 1, nEta = Int_t(etaBins.size()) - 1;

  std::printf("==> templates %s : %d |eta| x %d pt, flavour %d\n",
              tmplPath.c_str(), nEta, nPt, wantFlavor);

  // ---- truth: accumulate per cell straight from NanoAOD ---------------------
  std::vector<Stats> truth(nEta * nPt);
  std::vector<Stats> truthCvB(nEta * nPt); // (derived, stored) pairs
  std::vector<Long64_t> truthN(nEta * nPt, 0);
  Double_t maxCvBResid = 0;

  TChain chain("Events");
  for(const auto &in : inputs) chain.Add(in.c_str());

  TTreeReader reader(&chain);
  TTreeReaderArray<Float_t> jetPt(reader, "Jet_pt");
  TTreeReaderArray<Float_t> jetEta(reader, "Jet_eta");
  TTreeReaderArray<Int_t> jetFlav(reader, "Jet_hadronFlavour");
  TTreeReaderArray<Float_t> jetB(reader, "Jet_btagDeepFlavB");
  TTreeReaderArray<Float_t> jetCvL(reader, "Jet_btagDeepFlavCvL");
  TTreeReaderArray<Float_t> jetCvB(reader, "Jet_btagDeepFlavCvB");

  while(reader.Next())
    for(size_t j = 0; j < jetPt.GetSize(); ++j)
    {
      if(jetFlav[j] != wantFlavor) continue;
      const Double_t B = jetB[j], CvL = jetCvL[j], CvB = jetCvB[j];
      if(B < 0 || B > 1 || CvL < 0 || CvL > 1 || CvB < 0 || CvB > 1) continue;

      const Int_t ie = FindCell(etaBins, std::fabs(jetEta[j]));
      const Int_t ip = FindCell(ptBins, jetPt[j]);
      if(ie < 0 || ip < 0) continue;

      const Int_t k = ie * nPt + ip;
      truth[k].Add(B, CvL);
      truthCvB[k].Add(DeriveCvB(B, CvL), CvB);
      ++truthN[k];
      maxCvBResid = std::max(maxCvBResid, std::fabs(DeriveCvB(B, CvL) - CvB));
    }

  std::printf("==> identity check on NanoAOD itself:\n");
  std::printf("    max |CvB_derived - CvB_stored| over all jets = %.2e", maxCvBResid);
  std::printf("   %s\n", maxCvBResid < 1e-2 ? "(storage precision -- identity holds)" : "(TOO LARGE)");

  // ---- pick the fattest cell so the comparison is statistically meaningful --
  Int_t best = -1;
  for(Int_t k = 0; k < nEta * nPt; ++k)
    if(best < 0 || truthN[k] > truthN[best]) best = k;
  if(best < 0 || truthN[best] == 0) { std::cerr << "ERROR: no jets of flavour " << wantFlavor << "\n"; return 1; }

  const Int_t ie = best / nPt, ip = best % nPt;
  TString hname = TString::Format("tmpl_f%d_eta%d_pt%d", wantFlavor, ie, ip);
  TH2F *h = dynamic_cast<TH2F *>(ft->Get(hname));
  if(!h) { std::cerr << "ERROR: " << hname << " not found\n"; return 1; }

  std::printf("==> densest cell %s : %.2f<|eta|<%.2f, %.0f<pt<%.0f, %lld jets\n",
              hname.Data(), etaBins[ie], etaBins[ie + 1], ptBins[ip], ptBins[ip + 1], truthN[best]);

  // ---- sample the joint, and the marginals independently for contrast ------
  TRandom3 rng(12345);
  TH1D *mB = h->ProjectionX("mB");
  TH1D *mCvL = h->ProjectionY("mCvL");

  Stats joint, indep;
  Long64_t badIndep = 0;
  Double_t x, y;
  for(Long64_t i = 0; i < nSample; ++i)
  {
    h->GetRandom2(x, y, &rng);
    joint.Add(x, y);

    // Two independent 1-D templates -- what running the module twice would give.
    indep.Add(mB->GetRandom(&rng), mCvL->GetRandom(&rng));
  }
  (void)badIndep;

  std::printf("\n%-22s %10s %10s %10s\n", "", "truth", "joint", "independent");
  std::printf("%-22s %10.4f %10.4f %10.4f\n", "mean B",   truth[best].MeanX(), joint.MeanX(), indep.MeanX());
  std::printf("%-22s %10.4f %10.4f %10.4f\n", "mean CvL", truth[best].MeanY(), joint.MeanY(), indep.MeanY());
  std::printf("%-22s %10.4f %10.4f %10.4f\n", "rho(B, CvL)", truth[best].Rho(), joint.Rho(), indep.Rho());

  const Double_t dRhoJoint = std::fabs(joint.Rho() - truth[best].Rho());
  const Double_t dRhoIndep = std::fabs(indep.Rho() - truth[best].Rho());
  std::printf("\n    correlation error: joint %.4f, independent %.4f\n", dRhoJoint, dRhoIndep);

  // ---- derived CvB from sampled jets vs truth CvB --------------------------
  Stats cvbCmp;
  for(Long64_t i = 0; i < nSample; ++i)
  {
    h->GetRandom2(x, y, &rng);
    cvbCmp.Add(DeriveCvB(x, y), 0.0);
  }
  std::printf("    mean CvB: truth %.4f, derived-from-sampled %.4f\n",
              truthCvB[best].MeanY(), cvbCmp.MeanX());

  ft->Close();
  delete ft;
  return 0;
}
