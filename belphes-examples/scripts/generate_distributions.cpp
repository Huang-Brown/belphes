// File: generate_distributions.cpp
/*
    For our examples, make two .root files. File 1 "distrib_b.root" contains three 
    gaussians, and file 2 "distrib_nonb.root" contains three exponential distributions.

    Run by `root -l -q 'generate_distributions.cpp'`.
*/
#include <TRandom3.h>
#include <TH1D.h>
#include <TFile.h>
#include <TString.h>

void generate_distributions() {
  // Random engine
  TRandom3 rnd(0);

  // --- 1) b_distrib.root with three Gaussians σ = 1, 2, 3
  TFile *fb = TFile::Open("../example_distributions/distrib_b.root", "RECREATE");
  for (int i = 1; i <= 3; ++i) {
    // Create histogram named "hist_i"
    TH1D *h = new TH1D(
      TString::Format("hist_%d", i),
      TString::Format("Gaussian #sigma = %d; x; Entries", i),
      100,               // 100 bins
      -5.0 * i,          // x-min
      +5.0 * i           // x-max
    );
    // Fill with 100 000 entries from N(0,i)
    for (int n = 0; n < 100000; ++n)
      h->Fill(rnd.Gaus(0.0, i));
    h->Write();
  }
  fb->Close();

  // --- 2) nonb_distrib.root with three exponentials λ = 1, 2, 3
  // TRandom3::Exp(tau) generates PDF ∝ exp(–x/τ), so λ = 1/τ ⇒ τ = 1/λ
  TFile *fn = TFile::Open("../example_distributions/distrib_nonb.root", "RECREATE");
  for (int i = 1; i <= 3; ++i) {
    TH1D *h = new TH1D(
      TString::Format("hist_%d", i),
      TString::Format("Exponential #lambda = %d; x; Entries", i),
      100,   // bins
      0.0,   // x-min
      10.0   // x-max
    );
    double tau = 1.0 / double(i);
    for (int n = 0; n < 100000; ++n)
      h->Fill(rnd.Exp(tau));
    h->Write();
  }
  fn->Close();
}
