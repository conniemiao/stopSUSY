#!/usr/bin/env python

from root_pandas import read_root
import pandas as pd
import sys
import os

branches = ['met_pt','met_phi', 'njets','nbtag','jet_ht', 'lep1_pt', 'lep2_pt', \
        'lep1_mt', 'lep2_mt', 'mt_tot', 'mt_sum', 'm_eff']

# filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/stopCut_03Sig_elel.root"
# filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/stopCut_all_Bkgd_DY0Jets_MLL-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.root"
# filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/stopCut_all_Bkgd_DY1Jets_MLL-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.root"
filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/stopCut_all_Bkgd_TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.root"

# df = read_root(filename,"tSig0", columns=branches)
# df = read_root(filename,"tSig1", columns=branches)
df = read_root(filename,"tBkgd", columns=branches)

# df.to_hdf("/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/sig0_elel.hdf",key="tSig0")
# df.to_hdf("/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/sig1_elel.hdf",key="tSig1")
# df.to_hdf("/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/all_Bkgd_DY0Jets_MLL-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.hdf",key="tBkgd")
df.to_hdf("/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/stopCut_all_Bkgd_TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.hdf",key="tBkgd")

