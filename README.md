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
- Adjust which channels and processes to submit in `makeAllNtuples.py`:
```
bash makeAllNtuples.py all
```
3. (Fast) hadd the ntuples from step 2 so that there is only 1 ntuple for each subprocess. This will submit a condor job for each subprocess listed in bkgd_fileRedirector, for each channel.
- Adjust which channels and processes to submit in `haddAllNtuples.sh`. Then execute:
```
bash haddAllNtuples.sh all
```
  
The next few steps can all be executed using the same command, but different sections of `plotAllVars.sh` should be commented out:
```
bash plotAllVars.sh all save
```

4. (Slow) Plot with QCD from MC (aka round 1 of plotting): `plot1D_qcdMC.py`. For each of the ABCD regions, it will produce a root file containing the canvases for all the control variables which uses MC for QCD as well as each of the individual histograms.  
- To submit all of these jobs to condor, uncomment section 1A "Normal 1d plots (QCD MC)" and adjust cuts, regions, and channels in `plotAllVars.sh` (note that all regions for a particular cut + channel need to have been completed before step 5 can be execute for that cut + channel). Then execute: `bash plotAllVars.sh all save`
5. (Very fast) Plot with QCD estimated from data using the ABCD method (aka round 2 of plotting): `plot1D_qcdData.py`. This will produce the same structure of root file as in step 4, but will only plot in the signal (B) region.
A
- This step is run locally. Uncomment section 1B and adjust cuts and channels in `plotAllVars.sh` as needed. Then execute: `bash plotAllVars.sh all save`
6. (Very fast) To save the plots as actual pngs from the root files outputted in steps 4 and 5, uncomment sections 4A and/or 4B in `plotAllVars.sh`, then execute `bash plotAllVars.sh all save`

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
