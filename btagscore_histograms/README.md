# B-tag Score Histograms

This folder has two histograms which contain the b-tag score distributions (in histograms)
of the `Jet_btagDeepFlavB` value from CMS open data: https://opendata.cern.ch/record/67727.
The two files are for the distribution of b-jets and non b-jets, respectively.

TODO: I should add a script later to convert .root files of events into histograms
suitable for sampling.

The naming convention of each histogram in the .root file is 
```
hist_eta[X]_pt[Y]
```
where X and Y indicate the n-th number of bins.