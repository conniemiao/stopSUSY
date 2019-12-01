# NOTE: requires 2 args, testMode {test, all} and displayMode {show, save} 
# Creates a plotting condorsub file for each possible plotting type, channel, last
# cut, etc., and submits the condorsub file to condor.

testMode=$1
if [[ "$testMode" == "test" ]]; then 
    channels=("elel")
    bkgdProcesses=("Diboson")
    cuts=("baseline")
    plotVars2D=("lep1_pt" "MET_pt""Jet_ht" "mt_tot")
    regions=("A")
elif [[ "$testMode" == "all" ]]; then
    # channels=("mumu" "muel" "elel")
    channels=("muel" "elel")
    bkgdProcesses=("TTBar" "TT+X" "Diboson" "W-Jets" "Drell-Yan" "Single-Top" "QCD")
    cuts=("baseline" "nJet<4")
    plotVars2D=("lep1_pt" "lep2_pt" "lep1_mt" "lep2_mt" "MET_pt" "lep1_eta" \
        "lep2_eta" "Jet_ht" "mt_tot" "mt_sum" "m_eff")
    # regions=("A" "B" "C" "D")
    regions=("any")
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
    echo "------------------------------- $channel -------------------------------"

    #--------------------------------------------------------------------------------#
    # SECTION 1A

    # Args to plot1D_qcdMC.py: testMode {test, all}, displayMode {show, save}, channel
    # {mumu, elel, muel}, lastcut, region {A,B,C,D}
    echo "------------------ Normal 1d plots (QCD MC) ------------------"
    for cut in "${cuts[@]}"
    do
        for region in "${regions[@]}"
        do
            bash createCondorsubPlotting.sh plot1D_qcdMC.py $testMode $displayMode \
                $channel $cut $region
            if [[ "$testMode" == "all" ]]; then 
                condor_submit condorsub_plotting
            else
                ./plot1D_qcdMC.py $testMode $displayMode $channel $cut $region
            fi
            # ./plot1D_qcdMC.py $testMode $displayMode $channel $cut $region
            echo
        done
    done

    #--------------------------------------------------------------------------------#
    # SECTION 1B

    # Args to plot1D_qcdData.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}, lastcut
    echo "------------------ Normal 1d plots (QCD Data) ------------------"
#     for cut in "${cuts[@]}"
#     do
#         # faster to just run locally than to submit to condor
#         ./plot1D_qcdData.py $testMode $displayMode $channel $cut
#         echo
#     done

    #--------------------------------------------------------------------------------#
    # SECTION 2

    # Args to plot2D.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}, lastcut, process, plotVarX, plotVarY
    echo
    echo "------------------ 2d plots ------------------"
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
#                     $channel $cut $process $plotVarX $plotVarY
#                 if [[ "$testMode" == "all" ]]; then 
#                     condor_submit condorsub_plotting
#                 else
#                     ./plot2D.py $testMode $displayMode $channel $cut $process \
#                         $plotVarX $plotVarY
#                 fi
# 
#             done
#         done
#     done

    #--------------------------------------------------------------------------------#

    echo
done

#--------------------------------------------------------------------------------#

echo
for channel in "${channels[@]}"
do
    echo "------------------------------- $channel -------------------------------"

    #--------------------------------------------------------------------------------#
    # SECTION 3A

    # Args to getPlots.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}, lastcut
    echo
    echo "------------------ Get 1d plots from qcd MC ------------------"
#     for cut in "${cuts[@]}"
#     do
#         for region in "${regions[@]}"
#         do
#             python getPlots.py $testMode $displayMode $channel $cut qcdMC $region
#         done
#     done

    #--------------------------------------------------------------------------------#
    # SECTION 3B

    # Args to getPlots.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}, lastcut
    echo
    echo "------------------ Get 1d plots from qcd data ------------------"
#     for cut in "${cuts[@]}"
#     do
#         python getPlots.py $testMode $displayMode $channel $cut qcdData
#     done

    #--------------------------------------------------------------------------------#
    # SECTION 4

    # Args to plotCutflows.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}
    echo
    echo "------------------ Get cutflows and piecharts ------------------"
#     for region in "${regions[@]}"
#     do
#         python plotCutflows.py $testMode $displayMode $channel $region
#     done

    #--------------------------------------------------------------------------------#

    echo
done

#--------------------------------------------------------------------------------#
