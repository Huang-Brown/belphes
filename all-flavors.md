# all-flavors: usage guide

Three tools, run in order. Step 1 pulls CMS Open Data, step 2 turns it into
templates, step 3 samples those templates inside Delphes.

```
fetch_opendata.sh  ->  build_templates  ->  PseudoDeepFlavScore
   NanoAOD .root        templates.root       Jet_btagDeepFlav{B,CvB,CvL}
```

Steps 1 and 2 live in `belphes-templates/`; step 3 is a Delphes module in
`modules/`. `belphes-templates/build_templates.sh` drives 1 and 2 together if
you don't want to run them separately.

## Environment

Everything needs ROOT on `PATH`. On this machine:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ~/CMSSW_15_0_0/src && cmsenv
cd /isilon/export/home/jhuan166/belphes
```

`~/CMSSW_15_0_0` (ROOT 6.32.11, gcc 12.3.1) is what this branch was built and
tested against. `~/CMSSW_14_0_10` and `~/Vcb/CMSSW_15_1_0_patch4` also exist
but are untried here.

Step 1 additionally needs `xrdcp` (XRootD client) and `python3`.

Do **not** run `./configure` before `make` — it regenerates the Makefile and
drops the hand-added `fmt` include/lib flags at the top.

---

## 1. `fetch_opendata.sh`

Queries the CERN Open Data API for a record, resolves its NanoAOD file list,
and downloads over XRootD into `data/`.

```bash
cd belphes-templates
./fetch_opendata.sh --list --limit 5 --smallest   # look before downloading
./fetch_opendata.sh --limit 4 --smallest          # ~23 MB, for testing
./fetch_opendata.sh --all                         # 68 files, 17.6 GB
```

| Option | Meaning |
|---|---|
| `--record N` | Open Data record id. Default `67727`. |
| `--limit N` | Download only the first N files of the selection. |
| `--smallest` | Order by file size ascending. Pair with `--limit` for tests. |
| `--all` | Take the whole record. |
| `--list` | Print the selection and exit; download nothing. |
| `--outdir DIR` | Where files land. Default `./data`. |

Notes:

- It **refuses to run without `--limit` or `--all`**, because the default
  record is 17.6 GB and an accidental full pull is expensive.
- Re-running **skips files already present at the published size**, so an
  interrupted fetch resumes rather than restarting.
- Record metadata is cached in `data/.record_<N>.json`. Delete it to re-query.

Default record 67727 is `TTJets_SingleLeptFromT_genMET-150`, UL16 NanoAODv9,
7.4 M events, 68 files. The four smallest are ~5.8 MB each, which is enough to
exercise the whole chain in under a minute.

---

## 2. `build_templates`

Reads NanoAOD, writes one 2-D template per (hadronFlavour, |eta| bin, pt bin)
cell into a single ROOT file.

```bash
./build_templates -c binning.conf -o templates.root data/*.root
```

| Option | Meaning |
|---|---|
| `-c, --config FILE` | Binning config. Default `binning.conf`. |
| `-o, --output FILE` | Output template file. Default `templates.root`. |
| *(positional)* | One or more NanoAOD `.root` inputs. |

Compile it by hand, or let `build_templates.sh` do it:

```bash
g++ -O2 -Wall build_templates.cpp $(root-config --cflags --libs) -o build_templates
```

### What you control, in `binning.conf`

| Key | Meaning |
|---|---|
| `PtBins` | Cell edges in pt. Any ascending, non-uniform list. |
| `AbsEtaBins` | Cell edges in \|eta\|. |
| `Flavors` | `Jet_hadronFlavour` values to build separate sets for: `5 4 0`. |
| `ScoreBinsB`, `ScoreBinsCvL` | Uniform score-axis resolution. |
| `ScoreEdgesB`, `ScoreEdgesCvL` | Explicit score edges; override the counts above. Must run 0 to 1. |
| `JetPtMin`, `JetAbsEtaMax` | Jet selection before binning. |
| `MinCellEntries` | Cells below this are flagged in the report. |

The pt/\|eta\| edges and flavour list are **written into the output file**, and
the Delphes module reads them back from there. Nothing restates the binning in
a card, so the two cannot drift apart. Change the binning, rebuild, done — no
card edit.

`binning_test.conf` is a deliberately coarse grid for smoke-testing on a few
files. Not for physics.

### Reading the report

```
    grid          : 3 flavour x 5 |eta| x 18 pt = 270 cells
    resolution    : 50 (B, variable) x 50 (CvL, variable) = 2500 bins/cell
==> jets seen 36381, filled 25677
    rejected: 10 outside grid, 0 unlisted flavour, 0 invalid score
    cells empty            : 4
    cells under 500        : 259
    sparsest cell          : tmpl_f5_eta4_pt17 (0 entries)
```

**Empty cells are the thing to watch.** `PseudoDeepFlavScore` refuses to load a
template file containing any, because `TH2::GetRandom2` hands back `(0,0)` for
an empty histogram — a legal-looking score that is silently wrong. Fix by
widening bins or adding input files.

### Why the templates are 2-D

DeepJet's grouped probabilities satisfy `B + C + L = 1`, with `CvB = C/(C+B)`
and `CvL = C/(C+L)`, so `(B, CvB, CvL)` has only **two** degrees of freedom.
Storing the joint `(B, CvL)` and deriving

```
C = CvL * (1 - B)        CvB = C / (C + B)
```

reproduces all three exactly. Verified against record 67727: derived CvB
matches stored CvB to `4.3e-4`, i.e. NanoAOD's float storage precision.

`(B, CvL)` is the pair to store. Going the other way — sampling `(B, CvB)` and
deriving CvL — divides by `(1 - CvB)` and loses all precision for high-B jets,
which is the region b-tagging exists to describe.

Sampling B and CvL from two separate 1-D templates would destroy their
correlation, which is large: `rho(B, CvL)` runs 0.69 to 0.84 by flavour, and
independent draws return `rho ~ 0`.

### Choosing the score resolution

`CvB = C/(C+B)` is nonlinear and most sensitive where B is small — exactly
where light and c jets sit — so uniform score bins bias the derived CvB for the
flavours a c-tagger cares about. Measured on 4 open-data files,
`|CvB_truth - CvB_derived|`:

| flavour | 50x50 uniform | 50x50 variable | 100x100 uniform |
|---------|---------------|----------------|-----------------|
| b       | 0.0010        | 0.0006         | 0.0011          |
| c       | 0.0153        | **0.0005**     | 0.0002          |
| light   | 0.0034        | **0.0001**     | 0.0026          |

`binning.conf` ships tanh-warped variable axes by default: 100x100 accuracy at
a quarter of the bins. That matters because bins-per-cell trades directly
against the statistics available in each cell.

### Checking your work

```bash
./validate_templates -t templates.root -f 5 data/*.root
```

`-f` selects `Jet_hadronFlavour` (5, 4 or 0); `-n` sets the sample count.
Per flavour it verifies the `B + C + L = 1` identity directly on NanoAOD, picks
the densest cell, and compares truth against joint sampling and against
independent marginals. `build_templates.sh` runs this automatically.

---

## 3. `PseudoDeepFlavScore`

The Delphes module. Add it to a card **after `JetFlavorAssociation`** (it needs
`jet->Flavor` populated) and before `TreeWriter`, or just use the ready card:

```bash
./DelphesHepMC2 cards/belphes_card_CMS_allflavors.tcl out.root in.hepmc
```

Card block:

```tcl
module PseudoDeepFlavScore PseudoDeepFlavScore {
  set JetInputArray JetEnergyScale/jets
  set TemplateFile  "belphes-templates/templates.root"
  set ClampPt       1
  set DefaultFlavor 0
  set RandomSeed    0
}
```

| Parameter | Meaning |
|---|---|
| `JetInputArray` | Jet collection to score. Default `JetEnergyScale/jets`. |
| `TemplateFile` | Output of `build_templates`. Binning is read from here. |
| `ClampPt` | `1`: jets above the top pt edge use the edge bin. `0`: leave unscored. |
| `DefaultFlavor` | Template used for jets that are neither b (5) nor c (4). Default `0`. |
| `RandomSeed` | Seeds a **private** `TRandom3`. `0` picks a random seed. |

### Output

Three branches on `Jet`:

```
Jet.Jet_btagDeepFlavB      DeepJet b+bb+lepb discriminant
Jet.Jet_btagDeepFlavCvB    DeepJet c vs b+bb+lepb
Jet.Jet_btagDeepFlavCvL    DeepJet c vs uds+g
```

Example:

```cpp
TTree *t = (TTree*) f->Get("Delphes");
t->Draw("Jet.Jet_btagDeepFlavB", "Jet.Flavor == 5 && Jet.PT > 30");
t->Draw("Jet.Jet_btagDeepFlavCvL:Jet.Jet_btagDeepFlavCvB", "Jet.Flavor == 4", "COLZ");
```

CMS c-tag working points are a 2-D region in the (CvL, CvB) plane — that plane
is meaningful here precisely because the two are drawn from one joint template
rather than sampled independently.

### Behaviour worth knowing

- **`-1` means "not computed."** It is what jets outside the grid carry, and
  also what `GenJet` and `FatJet` carry, since the module never touches those
  arrays. `0` is a legal score, so it can no longer double as "missing." (This
  changed on this branch — `Jet_btagDeepFlavB` previously defaulted to `0`.)
- **`ClampPt 1` is the sensible default.** A 2 TeV jet is real and taggable;
  the templates just stop. `|eta|` outside the grid is left unscored either
  way, since that is genuinely outside tracker acceptance.
- **Empty templates are fatal at `Init`**, not a silent `(0,0)` per jet.
- **Missing `JetFlavorAssociation` is reported.** Without it every jet has
  `Flavor == 0` and everything is sampled as light; the module says so at the
  end of the run instead of failing silently.
- **The flavour mapping is not one-to-one, by necessity.** Templates are keyed
  on CMS `Jet_hadronFlavour` (`5 / 4 / 0`, light = 0). Delphes' `jet->Flavor`
  is the highest parton PDG code in the cone, so light jets appear as 1, 2, 3,
  21 *or* 0. The module maps `5 -> 5`, `4 -> 4`, everything else to
  `DefaultFlavor`. A 150-event Z->bb run gives:

  ```
  Flavor==0 : 68    Flavor==3 :  6    Flavor==21 : 28
  Flavor==1 : 13    Flavor==4 :  3
  Flavor==2 : 17    Flavor==5 : 143
  ```

  so 132 of 278 jets reach the light template through the fallback. Intended,
  but the two conventions define "light" differently, which is a real (small)
  systematic rather than a bug.
- **The RNG is private.** Seeding it does not perturb the smearing and
  efficiency draws every other Delphes module makes off `gRandom`.

### Rebuilding after a change

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ~/CMSSW_15_0_0/src && cmsenv
cd /isilon/export/home/jhuan166/belphes
make -j"$(nproc)" HAS_PYTHIA8=1 VERBOSE=1 &> build_delphes_$(date +%Y%m%d).log
```

`Jet` is at `ClassDef(Jet, 7)` on this branch. Bump it again if you add members.

---

## Verified end to end

Built under `CMSSW_15_0_0` (`make -j HAS_PYTHIA8=1`, exit 0, no new warnings)
and run on 150 Z->bb events from `belphes-examples/pp_z_bb.hepmc`:

- `Jet ClassVersion = 7`; all three members visible to ROOT I/O on both `Jet`
  and `Candidate`.
- 262 jets, 182 scored, **80 unscored -- every one of them `|eta| > 2.5`**,
  i.e. the acceptance cut and nothing else.
- Derived CvB against `C/(C+B)` with `C = CvL*(1-B)`, measured on the output
  tree: mean `-1.1e-9`, RMS `1.5e-8`. Exact to float precision, as it must be.
- `GenJet` carries `-1` throughout, confirming "not computed" stays
  distinguishable from a genuine score of 0.

Reproduce with:

```bash
cd belphes-templates
./build_templates.sh --fetch 4 --smallest
./build_templates -c binning_test.conf -o templates.root data/*.root
cd ..
./DelphesHepMC2 cards/belphes_card_CMS_allflavors.tcl out.root in.hepmc
```

`RandomSeed 0` picks a random seed, so counts shift slightly run to run. Set
it non-zero for reproducibility.

---

## Relationship to `PseudoBTagScore`

The old module and `cards/belphes_card_CMS.tcl` are untouched and still work.
`PseudoDeepFlavScore` supersedes it: truth flavour is split b/c/light instead
of b/non-b, and B and CvL are drawn together instead of B alone. Use one or the
other in a given card, not both — they write the same
`Jet_btagDeepFlavB` field.

## Caveats on the default sample

Record 67727 is single-lepton ttbar with a `genMET > 150 GeV` filter:

- Every b-jet comes from `t -> Wb` — hard, central, well separated. None from
  gluon splitting `g -> bb`, where real tagger performance is markedly worse.
- Roughly half of hadronic W decays are `W -> cs`, so the light-jet sample is
  unusually c-rich for a generic process.

Neither is a bug. Both are reasons to check whether this sample matches the
process you are simulating before trusting absolute efficiencies.
