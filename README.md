## scharm to charm cutflow

simple cutflow for scharm to charm analysis

### what this is for
 
This code was designed to do several things: 
 - Serve as a reference: I've tried to keep `src/cutflow.cxx` as simple and readable as possible. Most of the cuts are hard-coded, and I've tried to use functions only when they make things more clear. 
 - Compile anywhere: Assuming you have ROOT and SUSYTools installed, and that `ROOTCOREDIR` is set, it should build with `make`. 
 - Test common tools: right now just the c-tagging SF code, which consists of three files: 
  - `src/CtagCalibration.cxx` and `include/CtagCalibration.hh`: a wrapper for `CalibrationDataInterface` (the official b-tagging scalefactor package), which applies the "medium" JetFitterCharm scale factors and hides some of the ugliness of the underlying package. 
  - `include/ctag_defs.hh`: lightweight miscellaneous c-tag definitions. 
 
### what this isn't for 
 
 This isn't a framework. It's not supposed to be used to fill histograms, create ntuples, optimize selections, lookup cross sections, or set limits. It's not supposed to have any dependencies that SUSYTools doesn't have. 

### simple instructions

If you have SUSYTools installed and `ROOTCOREDIR` is set, you should be able to type 

    make 

which will produce an executable called `cutflow`. Running 

    ./cutflow root_file1.root root_file2.root ...

will run the cutflow over these root files and print a summary. 

The executable expects several extra files in the run directory (it's probably easiest to sotflink them to the appropriate files): 

 - `grl.xml`: good runs list (only used with data)
 - `cdi.root`: b-tagging calibration file. For JetFitterCharm it should be `2013-Winter-rel17.2.1.4_MC12-83.root`

### what are all these other files?

Aside from the familiar `makefile` to build, and the `cutflow.cxx` file itself there are a few extra files kicking around: 
- `map_libs.sh` is used by `makefile` to link to SUSYTools
- `CtagCalibration` is a wrapper for the c-tagging calibration data interface, hopefully a bit simpler to use. 
- `CutCounter` is basically a map which keeps track of the order in which cuts are counted (so we can dump them in order)
- `LinkDef.hh` is used to tell ROOT how to read vectors from D3PDs
- `SmartChain` is derived from `TChain`, but attempts to fix a few of the more egregious design flaws: 

 - squash the normal 
```cxx
jet_something = 0; 
chain->SetBranchStatus("jet_something", 1); 
chain->SetBranchAddress("jet_something", &jet_something); 
```
song and dance into 
```cxx
chain->SetBranch("jet_something", &jet_something); 
```
because it's easier to write one line than three (and harder to make mistakes). 
 - Throw exceptions when branches are missing
 - Throw exceptions when files are corrupted

- `SusyBuffer` holds the information that's read out of the D3PD. It's something like what comes out of `MakeClass`, but with less garbage code. 
- `ctag_defs.hh` contains general definitions used by `CtagCalibration`. 
