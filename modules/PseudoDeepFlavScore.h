#ifndef PseudoDeepFlavScore_h
#define PseudoDeepFlavScore_h

/** \class PseudoDeepFlavScore
 *
 *  Assigns a jointly-sampled set of pseudo DeepJet discriminants
 *  (btagDeepFlavB, btagDeepFlavCvB, btagDeepFlavCvL) to each jet.
 *
 *  Supersedes PseudoBTagScore, which sampled btagDeepFlavB alone from a
 *  b / non-b pair of 1-D histograms.  Two things change here: the truth
 *  flavour axis is split b / c / light, and B and CvL are drawn together
 *  from one 2-D template so their correlation survives.  CvB is then
 *  derived rather than sampled -- see the .cc for why that is exact.
 *
 *  The (pt, |eta|) binning and the flavour list are read from the template
 *  file itself, not from the Delphes card, so the two cannot disagree.
 *
 *  \author J. Huang - Brown U, Providence
 *
 */

#include <map>
#include <vector>

#include "classes/DelphesModule.h"

class TH2;
class TObjArray;
class TRandom3;

class PseudoDeepFlavScore: public DelphesModule
{
public:
  PseudoDeepFlavScore();
  ~PseudoDeepFlavScore();

  void Init();    ///< Load templates and their binning from disk
  void Process(); ///< Sample (B, CvL) per jet and derive CvB
  void Finish();  ///< Clean up

private:
  TIterator *fItJetInputArray; //!
  const TObjArray *fJetInputArray; //!

  TRandom3 *fRandom; //!< private stream; does not disturb gRandom

  std::vector<Double_t> fPtBins; //!
  std::vector<Double_t> fAbsEtaBins; //!
  Int_t fNbinsPT; //!
  Int_t fNbinsAbsEta; //!

  Bool_t fClampPt; //!< score jets beyond the pt grid using the edge bin
  Int_t fDefaultFlavor; //!< template flavour used for unlisted jet flavours

  /// fTemplates[flavour][absEtaBin][ptBin]; owned by this module.
  std::map<Int_t, std::vector<std::vector<TH2 *>>> fTemplates; //!

  Long64_t fNScored; //!
  Long64_t fNUnscored; //!
  Long64_t fNNoFlavor; //!

  ClassDef(PseudoDeepFlavScore, 1)
};

#endif
