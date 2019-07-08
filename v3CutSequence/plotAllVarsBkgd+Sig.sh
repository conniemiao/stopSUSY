# Creates a plotSusyBkgd+Sig condorsub file for both cuts to nocut and cuts to
# njets<4 and for each dilepton channel and submits the condorsub file to condor.
# Args to plotSusyBkgd+Sig.py: test mode, display, same flav, mu pref, last cut, 
# plotVar1, plotVar2, ...

variables=("lep1_pt" "lep1_eta" "lep1_relIso" "lep2_pt" "njets" "jet_pt" \
    "jet_eta" "nbtag" "dR_lep1_jet" "lep1_mt" "lep2_mt" "met_pt" "jet_ht" \
    "mt_sum" "mt_tot" "m_eff" "jet_ht_div_sqrt_met" "mt_tot_div_sqrt_met" \
    "m_eff_div_sqrt_met")
cuts=("njets<4" "nocut")

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

# ./createCondorsubPlotting.sh plotCutflows.py 0 0 0 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh plotCutflows.py 0 0 1 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh plotCutflows.py 0 0 1 1
# condor_submit condorsub_plotting


# python plot2DMetPtGraphs.py 0 0 0 0 "nocut" lep1_pt met_pt
# echo
# echo
# python plot2DMetPtGraphs.py 0 0 1 0 "nocut" lep1_pt met_pt
# echo                                          
# echo                                          
# python plot2DMetPtGraphs.py 0 0 1 1 "nocut" lep1_pt met_pt
# echo
# echo
# python plot2DMetPtGraphs.py 0 0 0 0 "nocut" lep1_eta met_pt
# echo
# echo
# python plot2DMetPtGraphs.py 0 0 1 0 "nocut" lep1_eta met_pt
# echo
# echo
# python plot2DMetPtGraphs.py 0 0 1 1 "nocut" lep1_eta met_pt
# echo
# echo
# python plot2DMetPtGraphs.py 0 0 0 0 "nocut" lep1_mt met_pt
# echo
# echo
# python plot2DMetPtGraphs.py 0 0 1 0 "nocut" lep1_mt met_pt
# echo
# echo
# python plot2DMetPtGraphs.py 0 0 1 1 "nocut" lep1_mt met_pt
