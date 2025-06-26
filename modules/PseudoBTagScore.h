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

#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"   // for the Jet class

class TObjArray;
class DelphesFormula;
// class TH1D; // needed for distributions?

class PseudoBTagScore : public DelphesModule 
{
public:
  PseudoBTagScore();
  ~PseudoBTagScore();

  void Init();     ///< Load histogram from disk
  void Process();  ///< Sample one value per jet
  void Finish();   ///< Clean up

private:
  TIterator       *fItJetInputArray; //!
  const TObjArray *fJetInputArray; //!

  // TFile           *fFile;       ///< ROOT file holding histogram
  // TH1D            *fHistGauss;  ///< Histogram to sample N(0,1)

  ClassDef(PseudoBTagScore, 1)

  /* 
    For reference, BTagging.h included the following if statement for an 
    "fEfficiencyMap"
  */
  // #if !defined(__CINT__) && !defined(__CLING__)
  // std::map<std::vector<Float_t>, TH1F *> fEfficiencyMap; //!
  // #endif
};

#endif

