# stopSUSY
SUSY project with Alexis S2019-F2019

### Run2
v1: Working with Run 2 NanoAOD data.  
Workflow instructions to make control plots after cloning this repo. _NOTE: All scripts have commenting at top with more detailed instructions on use. Also, in general, for the commands listed below, changing "all" to "test" will typically run a smaller number of jobs locally rather than submitting them to condor._  
1. (Very fast) Given datasets listed in bkgd_datasets, sig_datasets, and data_fileRedirector, creates lists of all the file names for those datasets. This step isn't technically necessary as these ntuple list files have already been included in the repo, but to recreate these files, execute:
```
bash getDASNtuples.sh
```
2. (Very slow) Produce ntuples. This will submit a condor job for each ntuple for each subprocess listed in bkgd_fileRedirector, for each channel, applying very loose cuts. 
- Change these directories before running if necessary:
  - `myReferenceDataDir`, `myDataDir` in `makeNtuple.py`
  - `baseDir`, `condorLogDir`, `x509Adr` in createCondorsubNtuplingSubprocess.sh
- Adjust which channels and processes to submit in `makeAllNtuples.py`:
```
bash makeAllNtuples.py all
```
3. (Fast) hadd the ntuples from step 2 so that there is only 1 ntuple for each subprocess. This will submit a condor job for each subprocess listed in bkgd_fileRedirector, for each channel.
- Change these directories before running if necessary:
  - `myDataDir` in `haddSubprocess.sh`
  - `baseDir` and `condorLogDir` in `createCondorsubHadd.sh`
- Adjust which channels and processes to submit in `haddAllNtuples.sh`. Then execute:
```
bash haddAllNtuples.sh all
```
  
_The next few steps can all be executed using the same command, but different sections of `plotAllVars.sh` should be commented out:_
```
bash plotAllVars.sh all save
```
4. (Slow) Plot with QCD from MC (aka round 1 of plotting): `plot1D_qcdMC.py`. For each of the ABCD regions, it will produce a root file containing the canvases for all the control variables which uses MC for QCD as well as each of the individual histograms. It will also save 2 files with the cutflow information; the .txt file is human readable and the .hdf file is for use when drawing the cutflow.
- Change these directories before running if necessary:
  - `myDataDir`, `statsDir`, `imgDir` in `plot1D_qcdMC.py`
- To submit all of these jobs to condor, uncomment section 1A "Normal 1d plots (QCD MC)" and adjust cuts, regions, and channels in `plotAllVars.sh` (note that all 4 regions for a particular cut + channel need to have been completed before step 5 can be execute for that cut + channel). Then execute: `bash plotAllVars.sh all save`
- You can also run this on 1 specific channel, last cut, and region by executing
```
python plot1D_qcdMC.py all save [mumu/muel/elel] [lastcut] [A/B/C/D/any]
```
5. (Very fast) Plot with QCD estimated from data using the ABCD method (aka round 2 of plotting): `plot1D_qcdData.py`. This will produce the same structure of root file as in step 4, but will only plot in the signal (B) region.
A
- Change these directories before running if necessary:
  - `imgDir` in `plot1D_qcdData.py`
  - `baseDir` and `condorLogDir` in `createCondorsubPlotting.sh`
- This step is run locally. Uncomment section 1B and adjust cuts and channels in `plotAllVars.sh` as needed. Then execute: `bash plotAllVars.sh all save`
- You can also run this on 1 specific channel and last cut by executing
```
python plot1D_qcdData.py all save [mumu/muel/elel] [lastcut]
```
The plots in the root files outputted by steps 4 and 5 can be viewed directly from the root files (all the canvases are saved in the root files), or:
6. (Very fast) To save the plots as actual pngs from the root files outputted in steps 4 and 5, uncomment sections 3A and/or 3B in `plotAllVars.sh`, then execute `bash plotAllVars.sh all save`
- Change these directories before running if necessary:
  - `imgDir` in `getPlots.py`
- You can also run this on 1 specific channel and last cut by executing
```
python getPlots.py all save [mumu/muel/elel] [lastcut] [qcdMC/qcdData] [A/B/C/D/ if for qcdMC]
```
7. (Very fast) To draw the cutflows and piecharts (for baseline and nJet<4 cuts) from the cutflow stats outputted in step 4, uncomment section 4 in `plotAllVars.sh`, then execute `bash plotAllVars.sh all save` 
- Change these directories before running if necessary:
  - `statsDir`, `cutflowPlotsDir` in `plotCutflows.py`
- You can also run this on 1 specific channel and last cut by executing
```
python plotCutflows.py all save [mumu/muel/elel] [A/B/C/D]
```

### HLLHC
v3CutSequence:
Handling stacking bkgd from multiple different sources, condor submit files are now automatically created, cutflows first created as a text file before an additional script plots the cutflows and pie charts for bkgd, actually FINALLY fixing correctly scaling everything using genweights.
Played around with machine learning.

v2CutSequence:
Make ttbar bkgd ntuple separately from sig bkgd, started adding parameters into plotting/ntupling scripts to streamline condor submission process.

v1CutSequence:
Separate plotting into an ntuple maker with extra cuts (makes 1 ntuple with all the trees from sig and ttbar bkgd), started using stopSelection to select leptons. Able to plot cutflows.

v0prep:
Basic plotting some var directly from signal and ttbar bkgd raw files (no cuts or selection)
Preprocess root: rewriting Alexis's code into python, similar to what I did later on but with more steps and from a different format of ntuple trees (e.g. needs b-tag upgrade/downgrade stuff)
