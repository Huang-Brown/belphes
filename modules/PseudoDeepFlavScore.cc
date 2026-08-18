/** \class PseudoDeepFlavScore
 *
 *  Assigns a jointly-sampled set of pseudo DeepJet discriminants
 *  (btagDeepFlavB, btagDeepFlavCvB, btagDeepFlavCvL) to each jet, drawn from
 *  templates built by belphes-templates/build_templates from CMS open data.
 *
 *  Why two axes and not three
 *  --------------------------
 *  DeepJet's grouped output probabilities satisfy B + C + L = 1, with
 *
 *      B   = probb + probbb + problepb
 *      CvB = C / (C + B)
 *      CvL = C / (C + L)
 *
 *  so the triplet has only two degrees of freedom.  Sampling (B, CvL) jointly
 *  and deriving
 *
 *      C = CvL * (1 - B),   CvB = C / (C + B)
 *
 *  reproduces all three exactly, with no possibility of an inconsistent
 *  combination.  Verified against open-data record 67727: the derived CvB
 *  agrees with the stored one to 4.3e-4, i.e. NanoAOD's float precision.
 *
 *  (B, CvL) is the right pair to store.  Going the other way -- sampling
 *  (B, CvB) and deriving CvL -- divides by (1 - CvB) and loses all precision
 *  for high-B jets, which is exactly the region b-tagging cares about.
 *
 *  Sampling B and CvL from separate 1-D templates would destroy their
 *  correlation, which is large: rho(B, CvL) runs 0.71 to 0.84 depending on
 *  flavour, and independent draws return rho ~ 0.
 *
 *  \author J. Huang - Brown U, Providence
 *
 */

#include "modules/PseudoDeepFlavScore.h"
ClassImp(PseudoDeepFlavScore)

#include "classes/DelphesClasses.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <TObjArray.h>
#include <TRandom3.h>
#include <TString.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

PseudoDeepFlavScore::PseudoDeepFlavScore() :
  fItJetInputArray(nullptr),
  fJetInputArray(nullptr),
  fRandom(nullptr),
  fNbinsPT(0),
  fNbinsAbsEta(0),
  fClampPt(kTRUE),
  fDefaultFlavor(0),
  fNScored(0),
  fNUnscored(0),
  fNNoFlavor(0)
{
}

//------------------------------------------------------------------------------

PseudoDeepFlavScore::~PseudoDeepFlavScore()
{
}

//------------------------------------------------------------------------------

void PseudoDeepFlavScore::Init()
{
  // Opening a TFile changes gDirectory; put it back before returning so that
  // any module initialised after this one still books into the output file.
  TDirectory *savedDir = gDirectory;

  const char *templateFile = GetString("TemplateFile", "");
  TFile *file = TFile::Open(templateFile, "READ");
  if(!file || file->IsZombie())
  {
    throw runtime_error(Form(
      "PseudoDeepFlavScore: cannot open template file '%s'", templateFile));
  }

  //----------*----------*----------

  // The binning lives in the template file, so the card cannot contradict it.
  TVectorD *vPt = dynamic_cast<TVectorD *>(file->Get("PtBins"));
  TVectorD *vEta = dynamic_cast<TVectorD *>(file->Get("AbsEtaBins"));
  TVectorD *vFlavor = dynamic_cast<TVectorD *>(file->Get("Flavors"));

  if(!vPt || !vEta || !vFlavor)
  {
    file->Close();
    delete file;
    if(savedDir) savedDir->cd();
    throw runtime_error(Form(
      "PseudoDeepFlavScore: '%s' has no binning metadata (PtBins / AbsEtaBins / "
      "Flavors). Rebuild it with belphes-templates/build_templates.", templateFile));
  }

  fPtBins.resize(vPt->GetNrows());
  for(Int_t i = 0; i < vPt->GetNrows(); ++i) fPtBins[i] = (*vPt)[i];
  fNbinsPT = Int_t(fPtBins.size()) - 1;

  fAbsEtaBins.resize(vEta->GetNrows());
  for(Int_t i = 0; i < vEta->GetNrows(); ++i) fAbsEtaBins[i] = (*vEta)[i];
  fNbinsAbsEta = Int_t(fAbsEtaBins.size()) - 1;

  if(fNbinsPT < 1 || fNbinsAbsEta < 1 ||
     !is_sorted(fPtBins.begin(), fPtBins.end()) ||
     !is_sorted(fAbsEtaBins.begin(), fAbsEtaBins.end()))
  {
    file->Close();
    delete file;
    if(savedDir) savedDir->cd();
    throw runtime_error("PseudoDeepFlavScore: template binning is empty or not ascending");
  }

  //----------*----------*----------

  // Detach each template from the file, then close it -- the standard Delphes
  // pattern (cf. TrackSmearing), so nothing depends on the file staying open.
  for(Int_t f = 0; f < vFlavor->GetNrows(); ++f)
  {
    const Int_t flavor = Int_t((*vFlavor)[f]);
    auto &grid = fTemplates[flavor];
    grid.resize(fNbinsAbsEta);

    for(Int_t i = 0; i < fNbinsAbsEta; ++i)
    {
      grid[i].resize(fNbinsPT, nullptr);
      for(Int_t j = 0; j < fNbinsPT; ++j)
      {
        TString name = TString::Format("tmpl_f%d_eta%d_pt%d", flavor, i, j);
        TH2 *hist = dynamic_cast<TH2 *>(file->Get(name));
        if(!hist)
        {
          file->Close();
          delete file;
          if(savedDir) savedDir->cd();
          throw runtime_error(Form(
            "PseudoDeepFlavScore: template %s missing from %s", name.Data(), templateFile));
        }

        // An empty template makes GetRandom2 hand back (0,0) for every jet in
        // the cell -- a legal-looking score that is silently wrong.  Refuse.
        if(hist->GetEntries() <= 0 || hist->Integral() <= 0.0)
        {
          file->Close();
          delete file;
          if(savedDir) savedDir->cd();
          throw runtime_error(Form(
            "PseudoDeepFlavScore: template %s is empty. Rebuild with coarser bins "
            "or more input files.", name.Data()));
        }

        hist->SetDirectory(nullptr);
        grid[i][j] = hist;
      }
    }
  }

  file->Close();
  delete file;

  //----------*----------*----------

  fClampPt = GetBool("ClampPt", true);
  fDefaultFlavor = GetInt("DefaultFlavor", 0);

  if(fTemplates.find(fDefaultFlavor) == fTemplates.end())
  {
    if(savedDir) savedDir->cd();
    throw runtime_error(Form(
      "PseudoDeepFlavScore: DefaultFlavor %d has no templates in %s",
      fDefaultFlavor, templateFile));
  }

  // A private stream: seeding this must not perturb every other module's
  // smearing, which is what touching gRandom would do.
  fRandom = new TRandom3(GetInt("RandomSeed", 0));

  fJetInputArray = ImportArray(GetString("JetInputArray", "JetEnergyScale/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  if(savedDir) savedDir->cd();
}

//------------------------------------------------------------------------------

void PseudoDeepFlavScore::Finish()
{
  if(fItJetInputArray)
  {
    delete fItJetInputArray;
    fItJetInputArray = nullptr;
  }

  for(auto &entry : fTemplates)
    for(auto &row : entry.second)
      for(TH2 *hist : row) delete hist;
  fTemplates.clear();

  if(fRandom)
  {
    delete fRandom;
    fRandom = nullptr;
  }

  // Quiet unless something is actually off, so a normal run stays clean.
  if(fNScored > 0 && fNNoFlavor == fNScored)
  {
    cout << "** WARNING: PseudoDeepFlavScore scored " << fNScored
         << " jets but none carried a b or c flavour. Is JetFlavorAssociation "
            "in the ExecutionPath before this module?" << endl;
  }
  if(fNUnscored > 0)
  {
    cout << "** INFO: PseudoDeepFlavScore left " << fNUnscored << " of "
         << (fNScored + fNUnscored) << " jets unscored (outside the template grid)." << endl;
  }
}

//------------------------------------------------------------------------------

void PseudoDeepFlavScore::Process()
{
  Candidate *jet;

  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    const Double_t pt = jetMomentum.Pt();
    const Double_t absEta = std::abs(jetMomentum.Eta());

    // |eta| outside the grid means outside tracker acceptance: genuinely
    // untaggable, so leave it flagged rather than extrapolating.
    Int_t etaBin = -1;
    if(absEta >= fAbsEtaBins.front() && absEta < fAbsEtaBins.back())
    {
      etaBin = Int_t(upper_bound(fAbsEtaBins.begin(), fAbsEtaBins.end(), absEta) - fAbsEtaBins.begin()) - 1;
    }

    // pt is different: a 2 TeV jet is real and taggable, the templates just
    // stop. Clamping to the edge bin beats declaring it untaggable.
    Int_t ptBin = -1;
    if(pt >= fPtBins.front() && pt < fPtBins.back())
    {
      ptBin = Int_t(upper_bound(fPtBins.begin(), fPtBins.end(), pt) - fPtBins.begin()) - 1;
    }
    else if(fClampPt)
    {
      ptBin = (pt < fPtBins.front()) ? 0 : fNbinsPT - 1;
    }

    if(etaBin < 0 || ptBin < 0)
    {
      jet->Jet_btagDeepFlavB = -1.0;
      jet->Jet_btagDeepFlavCvB = -1.0;
      jet->Jet_btagDeepFlavCvL = -1.0;
      ++fNUnscored;
      continue;
    }

    //----------*----------*----------

    // Delphes' Flavor is the highest parton PDG code in the cone; anything
    // that is not b or c shares the light template, as BTagging does.
    auto itFlavor = fTemplates.find(Int_t(jet->Flavor));
    if(itFlavor == fTemplates.end())
    {
      itFlavor = fTemplates.find(fDefaultFlavor);
      ++fNNoFlavor;
    }

    Double_t B = 0.0, CvL = 0.0;
    itFlavor->second[etaBin][ptBin]->GetRandom2(B, CvL, fRandom);

    // B + C + L = 1 makes the third discriminant exact, not another draw.
    const Double_t C = CvL * (1.0 - B);
    const Double_t den = C + B;
    const Double_t CvB = (den > 0.0) ? C / den : 0.0;

    jet->Jet_btagDeepFlavB = Float_t(B);
    jet->Jet_btagDeepFlavCvL = Float_t(CvL);
    jet->Jet_btagDeepFlavCvB = Float_t(CvB);
    ++fNScored;
  }
}

//------------------------------------------------------------------------------
