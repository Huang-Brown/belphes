#ifndef PseudoBTagScore_h
#define PseudoBTagScore_h

/** \class PseudoBTagScore
 *
 *  Give pseudo b-tag scores to each jet by sampling from a given distribution
 *  depending on the ground truth flavor, \eta, and PT of each jet.
 *
 *  \author J. Huang - Brown U, Providence
 *
 */

#include <vector>          // for std::vector
#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"   // for the Jet class

class TH1;                // forward-declare ROOT histogram
class TObjArray;

class PseudoBTagScore : public DelphesModule 
{
public:
  PseudoBTagScore();
  ~PseudoBTagScore();

  void Init();     ///< Load histogram from disk
  void Process();  ///< Sample one value per jet
  void Finish();   ///< Clean up

private:
  TIterator            *fItJetInputArray; //!
  const TObjArray      *fJetInputArray;   //!

  TFile                *fFile_b;          //!
  TFile                *fFile_nonb;       //!

  std::vector<Float_t>  fPtBins;          //!
  std::vector<Float_t>  fAbsEtaBins;      //!
  Int_t                 fNbinsPT;         //!
  Int_t                 fNbinsAbsEta;     //!

  std::vector<std::vector<TH1*>> fHists_b;    //!
  std::vector<std::vector<TH1*>> fHists_nonb; //!

  ClassDef(PseudoBTagScore, 1)
};

#endif

