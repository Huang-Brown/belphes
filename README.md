[![CI](https://github.com/delphes/delphes/actions/workflows/ci.yml/badge.svg)](https://github.com/delphes/delphes/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/21390046.svg)](https://zenodo.org/badge/latestdoi/21390046)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/delphes.svg)](https://anaconda.org/conda-forge/delphes)

# Belphes

**Author**: Jason Huang (jiashu_huang@brown.edu), Department of Physics, Brown University


🅱️elphes is a branch of Delphes—a C++ framework performing a fast multipurpose detector response simulation (see https://delphes.github.io). For each jet object, it samples
a pseudo b-tagging score.

## Quick start with Belphes

This branch assumes that you are using root and other softwares from some CMSSW.
After you have downloaded the files, go to a SCRAM-based area (CMSSW_XX_X_X/src)
and run 
```
scram tool info fmt
```

For example, for CMSSW_15_0_0, the system produces the following:
```
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Name : fmt
Version : 10.2.1-e35fd1db5eb3abc8ac0452e8ee427196
Revision : 1
++++++++++++++++++++

FMT_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fmt/10.2.1-e35fd1db5eb3abc8ac0452e8ee427196
INCLUDE=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fmt/10.2.1-e35fd1db5eb3abc8ac0452e8ee427196/include
LIB=fmt
LIBDIR=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fmt/10.2.1-e35fd1db5eb3abc8ac0452e8ee427196/lib
ROOT_INCLUDE_PATH=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fmt/10.2.1-e35fd1db5eb3abc8ac0452e8ee427196/include
```

Then, go into Makefile. Under the first comment, paste what follows after `INCLUDE`
 and `LIBDIR`:
```
CXXFLAGS += -I{$INCLUDE}
LDFLAGS += -L{$LIBDIR} -lfmt
```

Also note that Makefile has been modified
```
CXXFLAGS += -std=c++20 -I$(subst :, -I,$(CMSSW_FWLITE_INCLUDE_PATH)) ## allow c++20
```
so it may be compiled using a new C++ standard.

Now you can `make`. Do *NOT* do `./configure`, as it will reset Makefile. 

Finally, we can run Delphes:

```
./DelphesHepMC3
```

Command line parameters:

```
./DelphesHepMC3 config_file output_file [input_file(s)]
  config_file - configuration file in Tcl format
  output_file - output file in ROOT format,
  input_file(s) - input file(s) in HepMC format,
  with no input_file, or when input_file is -, read standard input.
```

Let's simulate some Z->ee events:

```
wget http://cp3.irmp.ucl.ac.be/downloads/z_ee.hep.gz
gunzip z_ee.hep.gz
./DelphesSTDHEP cards/delphes_card_CMS.tcl delphes_output.root z_ee.hep
```

or

```
curl -s http://cp3.irmp.ucl.ac.be/downloads/z_ee.hep.gz | gunzip | ./DelphesSTDHEP cards/delphes_card_CMS.tcl delphes_output.root
```

For the original documentation of Delphes, please visit https://delphes.github.io/workbook

<!-- ## Configure Delphes on lxplus.cern.ch

```
git clone https://github.com/delphes/delphes Delphes

cd Delphes

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh

make
``` -->

<!-- 

# Simple analysis using TTree::Draw

Now we can start [ROOT](https://root.cern) and look at the data stored in the output ROOT file.

Start ROOT and load Delphes shared library:

```
root -l
gSystem->Load("libDelphes");
```

Open ROOT file and do some basic analysis using Draw or TBrowser:

```
TFile *f = TFile::Open("delphes_output.root");
f->Get("Delphes")->Draw("Electron.PT");
TBrowser browser;
```

Notes:

- `Delphes` is the tree name. It can be learned e.g. from TBrowser.
- `Electron` is the branch name.
- `PT` is a variable (leaf) of this branch.

Complete description of all branches can be found in [doc/RootTreeDescription.html](doc/RootTreeDescription.html).

This information is also available at [this link](https://delphes.github.io/workbook/root-tree-description).

# Macro-based analysis

Analysis macro consists of histogram booking, event loop (histogram filling), histogram display.

Start ROOT and load Delphes shared library:

```
root -l
gSystem->Load("libDelphes");
```

Basic analysis macro:

```
{
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("delphes_output.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  // Book histograms
  TH1 *histElectronPT = new TH1F("electron pt", "electron P_{T}", 50, 0.0, 100.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // If event contains at least 1 electron
    if(branchElectron->GetEntries() > 0)
    {
      // Take first electron
      Electron *electron = (Electron*) branchElectron->At(0);

      // Plot electron transverse momentum
      histElectronPT->Fill(electron->PT);

      // Print electron transverse momentum
      cout << electron->PT << endl;
    }

  }

  // Show resulting histograms
  histElectronPT->Draw();
}
```

# More advanced macro-based analysis

The `examples` directory contains ROOT macros [Example1.C](examples/Example1.C), [Example2.C](examples/Example2.C) and [Example3.C](examples/Example3.C).

Here are the commands to run these ROOT macros:

```
root -l
.X examples/Example1.C("delphes_output.root");
```

or

```
root -l examples/Example1.C'("delphes_output.root")'
``` -->
