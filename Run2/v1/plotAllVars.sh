# Creates a plotting condorsub file for each possible plotting type, channel, last
# cut, etc., and submits the condorsub file to condor.

# declare -a channels=("mumu" "muel" "elel")
declare -a channels=("elel")

#--------------------------------------------------------------------------------#

cuts=("nocut" "njets<4")
# cuts=("nbtag<2" "MET>80")

for channel in "${channels[@]}"
do
    #--------------------------------------------------------------------------------#

    # Args to plot1D.py: testMode {test, all}, displayMode {show, save}, channel 
    # {mumu, elel, muel}, lastcut
    echo "Normal 1d plots:"
    for cut in "${cuts[@]}"
    do
        bash createCondorsubPlotting.sh plot1D.py all save $channel "$cut"
        condor_submit condorsub_plotting
    done

    #--------------------------------------------------------------------------------#

    # Args to generateCutflows.py: test mode, experimental, same flav, mu pref
    echo
    echo "Generate cutflow stats:"
    
    # bash createCondorsubPlotting.sh generateCutflows.py all $channel
    # condor_submit condorsub_plotting

    #--------------------------------------------------------------------------------#

done

# Args to generateCutflows.py: test mode, experimental, same flav, mu pref
echo
echo "Plot cutflows and pie charts:"

# python plotCutflows.py 0 0 0 0
# echo
# python plotCutflows.py 0 0 1 0
# echo
# python plotCutflows.py 0 0 1 1
# echo
# python plotCutflows.py 0 1 0 0
# echo
# python plotCutflows.py 0 1 1 0
# echo
# python plotCutflows.py 0 1 1 1
# echo

#--------------------------------------------------------------------------------#

# Args to plot2D.py: test mode, display, same flav, mu pref, last cut,
# bkgd process, plotVarX, plotVarY
echo
echo "2d plots:"
# processes=("W-Jets" "Drell-Yan" "Diboson" "Single-Top" "TT+X")
processes=("Single-Top")

# for process in "${processes[@]}"
# do
#     bash createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process jet_ht met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process jet_ht lep1_mt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process jet_ht mt_tot
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process m_eff met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process mt_sum met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 0 0 "nocut" $process mt_tot met_pt
#     condor_submit condorsub_plotting
#     
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process jet_ht mt_tot
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process jet_ht met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process jet_ht lep1_mt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process m_eff met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process mt_sum met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 0 "nocut" $process mt_tot met_pt
#     condor_submit condorsub_plotting
#     
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process jet_ht mt_tot
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process jet_ht met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process jet_ht lep1_mt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process m_eff met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process mt_sum met_pt
#     condor_submit condorsub_plotting
#     bash createCondorsubPlotting.sh plot2D.py 0 0 1 1 "nocut" $process mt_tot met_pt
#     condor_submit condorsub_plotting
# done

#--------------------------------------------------------------------------------#

