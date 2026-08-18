# belphes templates

Builds the DeepJet templates that `PseudoDeepFlavScore` samples from, in three
steps: fetch CMS Open Data, turn it into histograms, sample them inside Delphes.

```bash
cmsenv                                    # or any environment with root-config
./build_templates.sh --fetch 4 --smallest # ~23 MB, end-to-end in under a minute
./build_templates.sh --fetch 68           # the full 17.6 GB record
```

## Step 1 — `fetch_opendata.sh`

Queries the Open Data API for a record, resolves the NanoAOD file list, and
pulls files over XRootD.

```bash
./fetch_opendata.sh --list --limit 5 --smallest
./fetch_opendata.sh --limit 4 --smallest    # smallest files first, for testing
./fetch_opendata.sh --all                   # 68 files, 17.6 GB
```

Downloads land in `data/`. Re-running skips files already present at the
published size, so an interrupted fetch resumes. It refuses to run without
`--limit` or `--all`, since the default record is O(10) GB.

Default record is **67727** — `TTJets_SingleLeptFromT_genMET-150`, UL16
NanoAODv9, 7.4 M events. `--record N` points it elsewhere.

## Step 2 — `build_templates`

Reads NanoAOD, writes one 2-D template per (hadronFlavour, |eta| bin, pt bin)
cell. Everything configurable lives in `binning.conf`: pt edges, |eta| edges,
flavour list, score-axis resolution, jet selection.

```bash
./build_templates -c binning.conf -o templates.root data/*.root
```

The pt/|eta| edges and flavour list are **written into the output file**, and
the Delphes module reads them back from there. Nothing restates the binning in
a card, so the two cannot drift apart.

### Why the templates are 2-D and not 3-D

DeepJet's grouped probabilities satisfy `B + C + L = 1`, with

```
B   = probb + probbb + problepb
CvB = C / (C + B)
CvL = C / (C + L)
```

so `(B, CvB, CvL)` has only **two** degrees of freedom. Storing the joint
`(B, CvL)` and deriving

```
C = CvL * (1 - B)      CvB = C / (C + B)
```

reproduces all three exactly. Measured against record 67727, the derived CvB
agrees with the stored one to `4.3e-4` — NanoAOD's float storage precision.

`(B, CvL)` is the right pair to store. The other direction — sampling
`(B, CvB)` and deriving CvL — divides by `(1 - CvB)` and loses all precision
for high-B jets, which is the region b-tagging exists to describe.

### Why joint and not two 1-D templates

`rho(B, CvL)` is large — 0.69 to 0.84 depending on flavour. Sampling the two
marginals independently returns `rho ~ 0`, destroying it completely. See the
closure test output.

### Choosing the score resolution

`CvB = C/(C+B)` is nonlinear and most sensitive where B is small, which is
where light and c jets sit. Uniform score bins therefore bias the derived CvB
for exactly the flavours a c-tagger cares about. Measured on 4 open-data files
(`|CvB_truth - CvB_derived|`):

| flavour | 50x50 uniform | 50x50 variable | 100x100 uniform |
|---------|---------------|----------------|-----------------|
| b       | 0.0010        | 0.0006         | 0.0011          |
| c       | 0.0153        | **0.0005**     | 0.0002          |
| light   | 0.0034        | **0.0001**     | 0.0026          |

`binning.conf` ships the variable (tanh-warped) axes by default: they reach
100x100 accuracy at a quarter of the bins, and bins-per-cell trades directly
against the statistics available in each cell.

## Step 3 — `PseudoDeepFlavScore`

The Delphes module in `modules/`. Add it to a card after
`JetFlavorAssociation`, or use `cards/belphes_card_CMS_allflavors.tcl`:

```tcl
module PseudoDeepFlavScore PseudoDeepFlavScore {
  set JetInputArray JetEnergyScale/jets
  set TemplateFile  "belphes-templates/templates.root"
  set ClampPt       1
  set DefaultFlavor 0
  set RandomSeed    0
}
```

It writes three branches on `Jet`: `Jet_btagDeepFlavB`,
`Jet_btagDeepFlavCvB`, `Jet_btagDeepFlavCvL`.

Behaviour worth knowing:

- **Flavour** comes from Delphes' `jet->Flavor`. Anything that is not b (5) or
  c (4) uses the `DefaultFlavor` template, as `BTagging` does. If
  `JetFlavorAssociation` is missing from the ExecutionPath the module says so
  at the end of the run rather than silently tagging everything as light.
- **`-1` means "not computed"**, and is now the value carried by jets the
  module never sees (`GenJet`, `FatJet`) as well as jets outside the grid. `0`
  is a legal score, so it can no longer stand in for "missing".
- **`ClampPt 1`** scores jets above the top pt edge from the edge bin. A 2 TeV
  jet is real and taggable; the templates just stop. `|eta|` outside the grid
  is left unscored either way, since that is genuinely outside acceptance.
- **Empty templates are a fatal error at Init**, not a silent `(0,0)` per jet.
- The module uses a **private `TRandom3`**. Seeding it does not perturb the
  smearing and efficiency draws every other Delphes module makes off `gRandom`.

## `validate_templates`

Closure test, run automatically by `build_templates.sh`:

```bash
./validate_templates -t templates.root -f 5 data/*.root
```

Per flavour it checks the `B + C + L = 1` identity directly on NanoAOD, picks
the densest cell, and compares truth against joint sampling and against
independent marginals.

## Caveats on the default record

Record 67727 is single-lepton ttbar with a `genMET > 150 GeV` filter. Two
consequences for anything built from it:

- Every b-jet comes from `t -> Wb` — hard, central, well separated. None come
  from gluon splitting `g -> bb`, where the b's are collimated and real tagger
  performance is markedly worse.
- Roughly half of hadronic W decays are `W -> cs`, so the light-jet sample is
  unusually c-rich for a generic process. Splitting c out into its own template
  (which this pipeline does) matters more here than it would elsewhere.

Neither is a bug; both are reasons to think about whether this sample matches
the process you are simulating.
