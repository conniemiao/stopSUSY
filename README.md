# stopSUSY
SUSY project with Alexis S2019-F2019

### Run2
v4CutSequence:
Working with Run 2 NanoAOD data.

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
