## scharm to charm cutflow

simple cutflow for scharm to charm analysis

### what this is for
 
This code was designed to do several things: 
 - Serve as a reference: I've tried to keep `src/cutflow.cxx` as simple and readable as possible. Most of the cuts are hard-coded, and I've tried to use functions only when they make things more clear. 
 - Compile anywhere: Assuming you have ROOT and SUSYTools installed, and that `ROOTCOREDIR` is set, it should build with `make`. 
 - Test common tools: right now just the c-tagging SF code, which consists of three files: 
  - `src/CtagCalibration.cxx` and `include/CtagCalibration.hh`: a wrapper for `CalibrationDataInterface`, which applies the "medium" JetFitterCharm scale factors. 
  - `include/ctag_defs.hh`: lightweight miscellaneous c-tag definitions. 
 
### what this isn't for 
 
 This isn't a framework. It's not supposed to be used to fill histograms, create ntuples, optomize selections, lookup cross sections, or set limits. It's not supposed to have any dependencies that SUSYTools doesn't have. 
