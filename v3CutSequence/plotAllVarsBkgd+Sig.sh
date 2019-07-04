# plots all the variables in the ntuples and the cutflows for all the variables
# and for all the possible combinations of settings

variables=("lep1_pt" "lep1_eta" "lep1_relIso" "lep2_pt" "njets" "jet_pt" \
    "jet_eta" "nbtag" "dR_lep1_jet" "lep1_mt" "lep2_mt" "met_pt" "jet_ht" \
    "mt_sum" "mt_tot" "m_eff" "jet_ht_div_sqrt_met" "mt_tot_div_sqrt_met" \
    "m_eff_div_sqrt_met")

# plotSusyBkgd+Sig.py: testMode, displayMode, findingSameFlavor, muPreference, 
# lastcut, plotVar1, plotVar2, ...

#python plotSusyBkgd+Sig.py 0 0 0 0 "njets<4" "${variables[@]}"
#echo
#echo
#python plotSusyBkgd+Sig.py 0 0 1 0 "njets<4" "${variables[@]}"
#echo
#echo
python plotSusyBkgd+Sig.py 0 0 1 1 "njets<4" "${variables[@]}"

#python plotSusyBkgd+Sig.py 0 0 0 0 "nocut" "${variables[@]}"
#echo
#echo
# python plotSusyBkgd+Sig.py 0 0 1 0 "nocut" "${variables[@]}"
#echo
#echo
# python plotSusyBkgd+Sig.py 0 0 1 1 "nocut" "${variables[@]}"


#python plotCutflows.py 0 0 0 0
#echo
#echo
#python plotCutflows.py 0 0 1 0
#echo
#echo
python plotCutflows.py 0 0 1 1


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
