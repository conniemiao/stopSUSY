# Creates a plotSusyBkgd+Sig condorsub file for both cuts to nocut and cuts to
# njets<4 and for each dilepton channel and submits the condorsub file to condor.

variables=("lep1_pt" "lep1_eta" "lep1_relIso" "lep2_pt" "njets" "jet_pt" \
    "jet_eta" "nbtag" "dR_lep1_jet" "lep1_mt" "lep2_mt" "met_pt" "jet_ht" \
    "mt_sum" "mt_tot" "m_eff" "jet_ht_div_sqrt_met" "mt_tot_div_sqrt_met" \
    "m_eff_div_sqrt_met")
cuts=("njets<4" "nocut")

echo "Normal 1d plots:"
# Args to plotSusyBkgd+Sig.py: test mode, display, same flav, mu pref, last cut, 
# plotVar1, plotVar2, ...
for cut in "${cuts[@]}"
do
    ./createCondorsubPlotting.sh plotSusyBkgd+Sig.py 0 0 0 0 "$cut" \
        "${variables[@]}"
    condor_submit condorsub_plotting
    ./createCondorsubPlotting.sh plotSusyBkgd+Sig.py 0 0 1 0 "$cut" \
        "${variables[@]}"
    condor_submit condorsub_plotting
    ./createCondorsubPlotting.sh plotSusyBkgd+Sig.py 0 0 1 1 "$cut" \
        "${variables[@]}"
    condor_submit condorsub_plotting
done

echo
echo "Cutflows:"
# Args to plotCutflows.py: test mode, display, same flav, mu pref
# ./createCondorsubPlotting.sh plotCutflows.py 0 0 0 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh plotCutflows.py 0 0 1 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh plotCutflows.py 0 0 1 1
# condor_submit condorsub_plotting


echo
echo "2d plots:"
# Args to plot2DMetPtGraphs.py: test mode, display, same flav, mu pref, last cut,
# bkgd process, plotVarX, plotVarY
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 0 0 "nocut" jet_ht met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 0 0 "nocut" jet_ht lep1_mt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 0 0 "nocut" jet_ht mt_tot
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 0 0 "nocut" m_eff met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 0 0 "nocut" mt_sum met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 0 0 "nocut" mt_tot met_pt
condor_submit condorsub_plotting

./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 0 "nocut" jet_ht mt_tot
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 0 "nocut" jet_ht met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 0 "nocut" jet_ht lep1_mt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 0 "nocut" m_eff met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 0 "nocut" mt_sum met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 0 "nocut" mt_tot met_pt
condor_submit condorsub_plotting

./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 1 "nocut" jet_ht mt_tot
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 1 "nocut" jet_ht met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 1 "nocut" jet_ht lep1_mt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 1 "nocut" m_eff met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 1 "nocut" mt_sum met_pt
condor_submit condorsub_plotting
./createCondorsubPlotting.sh plot2DMetPtGraphs.py 0 0 1 1 "nocut" mt_tot met_pt
condor_submit condorsub_plotting
