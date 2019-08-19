# NOTE: requires 2 args, testMode {test, all} and displayMode {show, save} 
# Creates a plotting condorsub file for each possible plotting type, channel, last
# cut, etc., and submits the condorsub file to condor.

testMode=$1
if [[ "$testMode" == "test" ]]; then 
    channels=("elel")
    bkgdProcesses=("Diboson")
    cuts=("nocut")
    plotVars2D=("lep1_pt" "MET_pt""Jet_ht" "mt_tot")
    regions=("any")
elif [[ "$testMode" == "all" ]]; then
    channels=("mumu" "muel" "elel")
    bkgdProcesses=("TTBar" "TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top")
    cuts=("nocut" "njets<4")
    plotVars2D=("lep1_pt" "lep2_pt" "lep1_mt" "lep2_mt" "MET_pt" "lep1_eta" \
        "lep2_eta" "Jet_ht" "mt_tot" "mt_sum" "m_eff")
    regions=("A" "B" "C" "D")
else
    echo "need {test, all} as 1st arg to plotAllVars.sh"
    exit 1
fi

displayMode=$2
if [ "$displayMode" != "show" ] && [ "$displayMode" != "save" ]; then 
    echo "need {show, save} as 2nd arg to plotAllVars.sh"
    exit 1
fi

#--------------------------------------------------------------------------------#

for channel in "${channels[@]}"
do
    #--------------------------------------------------------------------------------#

    # Args to plot1D_qcdMC.py: testMode {test, all}, displayMode {show, save}, channel
    # {mumu, elel, muel}, lastcut, region {A,B,C,D}
    echo "Normal 1d plots (QCD MC):"
    for cut in "${cuts[@]}"
    do
        for region in "${regions[@]}"
        do
            bash createCondorsubPlotting.sh plot1D_qcdMC.py $testMode $displayMode \
                "$channel" "$cut" "$region"
            if [[ "$testMode" == "all" ]]; then 
                condor_submit condorsub_plotting
            else
                ./plot1D_qcdMC.py $testMode $displayMode "$channel" "$cut" "$region"
            fi
            echo
        done
    done

    #--------------------------------------------------------------------------------#

    # Args to generateCutflows.py: testMode {test, all}, channel {mumu, elel, muel}
    echo
    echo "Generate cutflow stats:"
#     bash createCondorsubPlotting.sh generateCutflows.py $testMode "$channel"
#     if [[ "$testMode" == "all" ]]; then 
#         condor_submit condorsub_plotting
#     else
#         ./generateCutflows.py $testMode "$channel"
#     fi

    #--------------------------------------------------------------------------------#

    # Args to plot2D.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}, lastcut, process, plotVarX, plotVarY
    echo
    echo "2d plots:"
#     for process in "${bkgdProcesses[@]}"
#     do
#         for plotVarX in "${plotVars2D[@]}"
#         do
#             for plotVarY in "${plotVars2D[@]}"
#             do
#                 if [ "$plotVarX" == "$plotVarY" ]; then
#                     continue
#                 fi
#                 bash createCondorsubPlotting.sh plot2D.py $testMode $displayMode \
#                     "$channel" "$cut" "$process" "$plotVarX" "$plotVarY"
#                 if [[ "$testMode" == "all" ]]; then 
#                     condor_submit condorsub_plotting
#                 else
#                     ./plot2D.py $testMode $displayMode "$channel" "$cut" "$process" \
#                         "$plotVarX" "$plotVarY"
#                 fi
# 
#             done
#         done
#     done

    #--------------------------------------------------------------------------------#

done

#--------------------------------------------------------------------------------#

echo
echo "Get 1d plots, cutflows, and pie charts:"
# for channel in "${channels[@]}"
# do
#     # Args to getPlots.py: testMode {test, all}, displayMode {show, save}, 
#     # channel {mumu, elel, muel}, lastcut
#     for cut in "${cuts[@]}"
#     do
#         for region in "${regions[@]}"
#         do
#             python getPlots.py $testMode $displayMode "$channel" "$cut" "$region"
#         done
#     done
# 
#     # Args to plotCutflows.py: testMode {test, all}, displayMode {show, save}, 
#     # channel {mumu, elel, muel}
#     python plotCutflows.py $testMode $displayMode "$channel"
# done

#--------------------------------------------------------------------------------#
