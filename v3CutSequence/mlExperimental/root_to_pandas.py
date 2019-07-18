from root_pandas import read_root
import pandas as pd
import sys
import os

branches = ['met_pt','met_phi', 'njets','nbtag','jet_ht', 'lep1_pt', 'lep2_pt', \
        'lep1_mt', 'lep2_mt', 'mt_tot', 'mt_sum', 'm_eff']

# filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/stopCut_03Sig_elel.root"
# filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/stopCut_all_Bkgd_ST_tch_14TeV_top_incl-powheg-pythia8-madspin_elel.root"
# filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/stopCut_all_Bkgd_ST_tch_14TeV_antitop_incl-powheg-pythia8-madspin_elel.root"
filename = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/stopCut_all_Bkgd_DY1Jets_MLL-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.root"

# df = read_root(filename,"tSig0", columns=branches)
# df = read_root(filename,"tSig1", columns=branches)
df = read_root(filename,"tBkgd", columns=branches)

# df.to_hdf("pandas_files/sig0.hdf",key="tSig0")
# df.to_hdf("pandas_files/sig1.hdf",key="tSig1")
# df.to_hdf("pandas_files/all_Bkgd_ST_tch_14TeV_top_incl-powheg-pythia8-madspin_elel.hdf",key="tBkgd")
df.to_hdf("/afs/cern.ch/work/c/cmiao/private/myDataSusy/pandas_files/all_Bkgd_DY1Jets_MLL-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.hdf",key="tBkgd")

