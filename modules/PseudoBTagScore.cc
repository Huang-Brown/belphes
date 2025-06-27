/** \class PseudoBTagScore
 *
 *  Give pseudo b-tag scores to each jet by sampling from a given distribution
 *  depending on the ground truth flavor, \eta, and PT of each jet.
 *
 *  \author J. Huang - Brown U, Providence
 *
 */

#include "modules/PseudoBTagScore.h"

#include "classes/DelphesFactory.h"
#include "classes/DelphesClasses.h"
#include <TFile.h>
#include <TH1.h>           // for TH1*
#include <TRandom.h>       // for gRandom
#include <stdexcept>       // for std::runtime_error
#include <algorithm>       // for std::upper_bound
#include <iterator>        // for std::distance
#include <vector>          // for std::vector
#include <cmath>           // std::abs

//------------------------------------------------------------------------------

// Constructor / Destructor
PseudoBTagScore::PseudoBTagScore() : 
  fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

PseudoBTagScore::~PseudoBTagScore() 
{
}

//------------------------------------------------------------------------------

void PseudoBTagScore::Init() 
{
  // load the files
  const char *file_b_dir    = GetString("Jet_btagDeepFlavB_file_b");
  const char *file_nonb_dir = GetString("Jet_btagDeepFlavB_file_nonb");

  fFile_b    = new TFile(file_b_dir,     "READ");
  fFile_nonb = new TFile(file_nonb_dir,  "READ");

  if(!fFile_b || fFile_b->IsZombie())
    throw std::runtime_error(Form(
      "PseudoBTagScore: cannot open the b file at %s", file_b_dir));

  if(!fFile_nonb || fFile_nonb->IsZombie())
    throw std::runtime_error(Form(
      "PseudoBTagScore: cannot open the non-b file at %s", file_nonb_dir));

  //----------*----------*----------

  // read bin edges - taken from ParticleDensity.cc
  ExRootConfParam paramPT = GetParam("PTBins");
  Int_t sizePT = paramPT.GetSize();
  fPtBins.reserve(sizePT);
  for(Int_t i = 0; i < sizePT; ++i) {
    fPtBins.push_back(paramPT[i].GetDouble());
  }
  fNbinsPT = sizePT - 1;

  ExRootConfParam paramAbsEta = GetParam("AbsEtaBins");
  Int_t sizeAbsEta = paramAbsEta.GetSize();
  fAbsEtaBins.reserve(sizeAbsEta);
  for(Int_t i = 0; i < sizeAbsEta; ++i) {
    fAbsEtaBins.push_back(paramAbsEta[i].GetDouble());
  }
  fNbinsAbsEta = sizeAbsEta - 1;

  //----------*----------*----------

  // load all the histograms
  fHists_b.resize(fNbinsAbsEta);
  fHists_nonb.resize(fNbinsAbsEta);

  for (Int_t i = 0; i < fNbinsAbsEta; ++i) {
    fHists_b[i].resize(fNbinsPT);
    fHists_nonb[i].resize(fNbinsPT);

    for (Int_t j = 0; j < fNbinsPT; ++j) {
      // The naming convention of each histogram in the .root file is  hist_eta[X]_pt[Y]
      TString nameB    = Form("hist_eta%d_pt%d", i, j);
      TString nameNonB = Form("hist_eta%d_pt%d", i, j);

      fHists_b[i][j]      = dynamic_cast<TH1*>(fFile_b->Get(nameB));
      fHists_nonb[i][j]   = dynamic_cast<TH1*>(fFile_nonb->Get(nameNonB));

      if (!fHists_b[i][j]) {
        throw std::runtime_error(Form(
          "PseudoBTagScore: histogram %s not found in b file", nameB.Data()));
      }
      if (!fHists_nonb[i][j]) {
        throw std::runtime_error(Form(
          "PseudoBTagScore: histogram %s not found in non-b file", nameNonB.Data()));
      }
    }
  }

  //----------*----------*----------

  // import input array 
  fJetInputArray   = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  //----------*----------*----------

  // OPTIONAL: seed the RNG
  UInt_t seed = GetInt("RandomSeed", 0u);
  if(seed != 0u) {
    gRandom->SetSeed(seed);
  }

}

//------------------------------------------------------------------------------

void PseudoBTagScore::Finish()
{
  // close off the input jet array
  if (fItJetInputArray) delete fItJetInputArray;

  // also need to close off all files
  if (fFile_b) 
  {
    fFile_b->Close();
    delete fFile_b;     fFile_b = nullptr;
  }
  if (fFile_nonb) 
  {
    fFile_nonb->Close();
    delete fFile_nonb;  fFile_nonb = nullptr;
  }
}

//------------------------------------------------------------------------------

void PseudoBTagScore::Process() 
{
  Candidate *jet; // Candidate is a Delphes class that can represent any object
  Double_t jet_pt,     jet_abseta;
  Int_t    jet_pt_bin, jet_abseta_bin;
  // loop over all input jets 
  fItJetInputArray->Reset();
  while ((jet = static_cast<Candidate *>(fItJetInputArray->Next()))) // while we pick up next jet from iterator
  {
    // obtain the pt and eta of the jet 
    const TLorentzVector &jetMomentum = jet->Momentum; // take 4-momentum of jet; TLorentzVector is outdated
    jet_pt    = jetMomentum.Pt();
    jet_abspt = std::abs(jetMomentum.Eta()); // we use abs eta because of axial symmetry in the detector

    // pre-set the bin
    jet_pt_bin     = -1;
    jet_abseta_bin = -1;

    //----------*----------*----------

    // find a suitable bin in pt and abseta for the jet 
    auto it_pt = std::upper_bound(fPtBins.begin(), fPtBins.end(), jet_pt);
    if (it_pt != fPtBins.begin() && it_pt != fPtBins.end()) {
      jet_pt_bin = std::distance(fPtBins.begin(), it_pt) - 1;
    }

    auto it_abseta = std::upper_bound(fAbsEtaBins.begin(), fAbsEtaBins.end(), jet_abseta);
    if (it_abseta != fAbsEtaBins.begin() && it_abseta != fAbsEtaBins.end()) {
      jet_abseta_bin = std::distance(fAbsEtaBins.begin(), it_abseta) - 1;
    }

    //----------*----------*----------

    // find the related histogram and sample from it
    if (jet_pt_bin>=0 && jet_abseta_bin>=0) {
      auto& histVec = (jet->Flavor==5 ? fHists_b : fHists_nonb);
      TH1*  hist    = histVec[jet_abseta_bin][jet_pt_bin];
      Float_t score = hist->GetRandom(); // Jet_btagDeepFlavB is set to be a Float_t
      jet->Jet_btagDeepFlavB = score;

    } else {
      jet->Jet_btagDeepFlavB = -1.0;

    }
  }
}


