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
  // load the file 
  const char* Jet_btagDeepFlavB_file_b, Jet_btagDeepFlavB_file_nonb;
  Jet_btagDeepFlavB_file_b    = GetString("Jet_btagDeepFlavB_file_b");
  Jet_btagDeepFlavB_file_nonb = GetString("Jet_btagDeepFlavB_file_nonb");

  if(!Jet_btagDeepFlavB_file_b || Jet_btagDeepFlavB_file_b->IsZombie())
    throw std::runtime_error(Form("PseudoBTagScore: cannot open %s", Jet_btagDeepFlavB_file_b));

  if(!Jet_btagDeepFlavB_file_nonb || Jet_btagDeepFlavB_file_b->IsZombie())
    throw std::runtime_error(Form("PseudoBTagScore: cannot open %s", Jet_btagDeepFlavB_file_nonb));

  // read bin edges, cache TH1* pointers, seed RNG...

  //----------------------------------------------------------------------------

  // import input array 
  fJetInputArray = ImportArray(GetString("JetInputArray")); // need to decide default variable
  fItJetInputArray = fJetInputArray->MakeIterator();

}

//------------------------------------------------------------------------------

void PseudoBTagScore::Finish()
{
  // close off the input jet array
  if (fItJetInputArray) delete fItJetInputArray;

  // also need to close off all histograms
  if (Jet_btagDeepFlavB_file_b)    Jet_btagDeepFlavB_file_b->Close();
  if (Jet_btagDeepFlavB_file_nonb) Jet_btagDeepFlavB_file_nonb->Close();
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


