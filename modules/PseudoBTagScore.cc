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
#include "classes/DelphesClasses.h"    // for Jet
#include <TFile.h>
#include <TH1D.h>

//------------------------------------------------------------------------------

// Constructor / Destructor
PseudoBTagScore::PseudoBTagScore() : 
  fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

PseudoBTagScore::~PseudoBTagScore() {
  if (fFile) {
    fFile->Close();
    delete fFile;
  }
}

//------------------------------------------------------------------------------

void PseudoBTagScore::Init() 
{
  // load the files
  const char* file_b_dir, file_nonb_dir;

  file_b_dir    = GetString("Jet_btagDeepFlavB_file_b");
  file_nonb_dir = GetString("Jet_btagDeepFlavB_file_nonb");

  fFile_b    = new TFile(file_b_dir,     "READ");
  fFile_nonb = new TFile(file_nonb_dir,  "READ");

  if(!fFile_b || fFile_b->IsZombie())
    throw std::runtime_error(Form("PseudoBTagScore: cannot open %s", file_b_dir));

  if(!fFile_nonb || fFile_nonb->IsZombie())
    throw std::runtime_error(Form("PseudoBTagScore: cannot open %s", file_nonb_dir));

  //----------------------------------------------------------------------------

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

  //----------------------------------------------------------------------------

  // import input array 
  fJetInputArray   = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

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
  Double_t pt, eta;

  // loop over all input jets 
  fItJetInputArray->Reset();
  while ((jet = static_cast<Candidate *>(fItJetInputArray->Next()))) // while we are able to pick up the next jet from the iterator
  {
    // obtain the pt and eta of the jet 
    const TLorentzVector &jetMomentum = jet->Momentum;
    pt = jetMomentum.Pt();
    eta = jetMomentum.Eta();

    // find a suitable bin for the jet 

    // apply the corresponding distribution 

    jet->Jet_btagDeepFlavB |= (gRandom->);
    

  }
}


