# Creates a plotSusyBkgd+Sig condorsub file for both cuts to nocut and cuts to
# njets<4 and for each dilepton channel and submits the condorsub file to condor.

variables=("lep1_pt" "lep1_eta" "lep1_relIso" "lep2_pt" "njets" "jet_pt" \
    "jet_eta" "nbtag" "dR_lep1_jet" "lep1_mt" "lep2_mt" "met_pt" "jet_ht" \
    "mt_sum" "mt_tot" "m_eff" "jet_ht_div_sqrt_met" "mt_tot_div_sqrt_met" \
    "m_eff_div_sqrt_met")

cuts=("njets<4" "nocut")

# processes=("W-Jets" "Drell-Yan" "Diboson" "Single-Top" "TT+X")
processes=("Single-Top")

#--------------------------------------------------------------------------------#

# Args to plotSusyBkgd+Sig.py: test mode, display, same flav, mu pref, last cut, 
# plotVar1, plotVar2, ...
echo "Normal 1d plots:"
# for cut in "${cuts[@]}"
# do
#     ./createCondorsubPlotting.sh plotSusyBkgd+Sig.py 0 0 0 0 "$cut" \
#         "${variables[@]}"
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plotSusyBkgd+Sig.py 0 0 1 0 "$cut" \
#         "${variables[@]}"
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plotSusyBkgd+Sig.py 0 0 1 1 "$cut" \
#         "${variables[@]}"
#     condor_submit condorsub_plotting
# done

#--------------------------------------------------------------------------------#

# Args to generateCutflows.py: test mode, experimental, same flav, mu pref
echo
echo "Generate cutflow stats:"
# ./createCondorsubPlotting.sh generateCutflows.py 0 0 0 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh generateCutflows.py 0 0 1 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh generateCutflows.py 0 0 1 1
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh generateCutflows.py 0 1 0 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh generateCutflows.py 0 1 1 0
# condor_submit condorsub_plotting
# ./createCondorsubPlotting.sh generateCutflows.py 0 1 1 1
# condor_submit condorsub_plotting

#--------------------------------------------------------------------------------#

# Args to generateCutflows.py: test mode, experimental, same flav, mu pref
echo
echo "Plot cutflows and pie charts:"
python plotCutflows.py 0 0 0 0
echo
python plotCutflows.py 0 0 1 0
echo
python plotCutflows.py 0 0 1 1
echo
python plotCutflows.py 0 1 0 0
echo
python plotCutflows.py 0 1 1 0
echo
python plotCutflows.py 0 1 1 1
echo

#--------------------------------------------------------------------------------#

# Args to plot2D.py: test mode, display, same flav, mu pref, last cut,
# bkgd process, plotVarX, plotVarY
echo
echo "2d plots:"
# for process in "${processes[@]}"
# do
#     ./createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process jet_ht met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process jet_ht lep1_mt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process jet_ht mt_tot
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process m_eff met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process mt_sum met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process mt_tot met_pt
#     condor_submit condorsub_plotting
#     
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process jet_ht mt_tot
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process jet_ht met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process jet_ht lep1_mt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process m_eff met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process mt_sum met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process mt_tot met_pt
#     condor_submit condorsub_plotting
#     
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process jet_ht mt_tot
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process jet_ht met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process jet_ht lep1_mt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process m_eff met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process mt_sum met_pt
#     condor_submit condorsub_plotting
#     ./createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process mt_tot met_pt
#     condor_submit condorsub_plotting
# done

#--------------------------------------------------------------------------------#

