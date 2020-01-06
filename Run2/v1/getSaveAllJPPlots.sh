# all the plots I need for the JP

echo "===================== ElEl ====================="
python getPlots.py "QCDMC_plot1D_all_ElEl_nJet<4_A.root" mt2
python getPlots.py "QCDMC_plot1D_all_ElEl_nJet<4_B.root" mt2
python getPlots.py "QCDMC_plot1D_all_ElEl_nJet<4_C.root" mt2
python getPlots.py "QCDMC_plot1D_all_ElEl_nJet<4_D.root" mt2
./plot1D_qcdData.py all save elel "nJet<4"
python getPlots.py "QCDData_plot1D_all_ElEl_nJet<4.root" mt2
python getPlots.py "fakeRegions_plot1D_all_ElEl_allCuts_sr.root" mt2
python getPlots.py "fakeRegions_plot1D_all_ElEl_allCuts_cr1a.root" mt2
python getPlots.py "fakeRegions_plot1D_all_ElEl_allCuts_cr1b.root" mt2
python getPlots.py "fakeRegions_plot1D_all_ElEl_allCuts_cr3.root" mt2 MET_pt mll
./plot1D_scaleDY_SR.py all save elel
python getPlots.py "fakeRegionsDYScaled_plot1D_all_ElEl_allCuts.root" mt2

echo "===================== MuEl ====================="
python getPlots.py "QCDMC_plot1D_all_MuEl_nJet<4_A.root" mt2
python getPlots.py "QCDMC_plot1D_all_MuEl_nJet<4_B.root" mt2
python getPlots.py "QCDMC_plot1D_all_MuEl_nJet<4_C.root" mt2
python getPlots.py "QCDMC_plot1D_all_MuEl_nJet<4_D.root" mt2
./plot1D_qcdData.py all save muel "nJet<4"
python getPlots.py "QCDData_plot1D_all_MuEl_nJet<4.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuEl_allCuts_sr.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuEl_allCuts_cr1a.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuEl_allCuts_cr1b.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuEl_allCuts_cr3.root" mt2 MET_pt mll
./plot1D_scaleDY_SR.py all save muel
python getPlots.py "fakeRegionsDYScaled_plot1D_all_MuEl_allCuts.root" mt2

echo "===================== MuMu ====================="
python getPlots.py "QCDMC_plot1D_all_MuMu_nJet<4_A.root" mt2
python getPlots.py "QCDMC_plot1D_all_MuMu_nJet<4_B.root" mt2
python getPlots.py "QCDMC_plot1D_all_MuMu_nJet<4_C.root" mt2
python getPlots.py "QCDMC_plot1D_all_MuMu_nJet<4_D.root" mt2
./plot1D_qcdData.py all save mumu "nJet<4"
python getPlots.py "QCDData_plot1D_all_MuMu_nJet<4.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuMu_allCuts_sr.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuMu_allCuts_cr1a.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuMu_allCuts_cr1b.root" mt2
python getPlots.py "fakeRegions_plot1D_all_MuMu_allCuts_cr3.root" mt2 MET_pt mll
./plot1D_scaleDY_SR.py all save mumu
python getPlots.py "fakeRegionsDYScaled_plot1D_all_MuMu_allCuts.root" mt2

echo
echo "To copy to local, execute scp -p cmiao@lxplus.cern.ch:/afs/cern.ch/user/c/cmiao/private/CMSSW_9_4_9/s2019_SUSY/plots/Run2/v1/plot1D/*.png ~/Desktop/stopSUSYPlots/Run2/v1/plot1D/"
