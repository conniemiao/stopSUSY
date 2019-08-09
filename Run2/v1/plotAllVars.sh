# Creates a plotting condorsub file for each possible plotting type, channel, last
# cut, etc., and submits the condorsub file to condor.

# declare -a channels=("mumu" "muel" "elel")
declare -a channels=("elel")

#--------------------------------------------------------------------------------#

cuts=("nocut" "njets<4")
# cuts=("nbtag<2" "MET>80")

# processes=("W-Jets" "Drell-Yan" "Diboson" "Single-Top" "TT+X")
processes=("Single-Top")

plotVars2D=("lep1_pt" "lep2_pt" "lep1_mt" "lep2_mt" "MET_pt" "lep1_eta" "lep2_eta" \
    "Jet_ht" "mt_tot" "mt_sum" "m_eff")

for channel in "${channels[@]}"
do
    #--------------------------------------------------------------------------------#

    # Args to plot1D.py: testMode {test, all}, displayMode {show, save}, channel 
    # {mumu, elel, muel}, lastcut
    echo "Normal 1d plots:"
    for cut in "${cuts[@]}"
    do
        bash createCondorsubPlotting.sh plot1D.py all save "$channel" "$cut"
        # condor_submit condorsub_plotting
    done

    #--------------------------------------------------------------------------------#

    # Args to generateCutflows.py: testMode {test, all}, channel {mumu, elel, muel}
    echo
    echo "Generate cutflow stats:"
    
    bash createCondorsubPlotting.sh generateCutflows.py all "$channel"
    # condor_submit condorsub_plotting

    #--------------------------------------------------------------------------------#

    # Args to plot2D.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}, lastcut, process, plotVarX, plotVarY
    echo
    echo "2d plots:"

    for process in "${processes[@]}"
    do
        for plotVarX in "${plotVars2D[@]}"
        do
            for plotVarY in "${plotVars2D[@]}"
            do
                if [ "$plotVarX" == "$plotVarY" ]; then
                    continue
                fi
                bash createCondorsubPlotting.sh plot2D.py all save "$channel" "$cut" \
                    "$process" "$plotVarX" "$plotVarY"
                # condor_submit condorsub_plotting
            done
        done
    done

    #--------------------------------------------------------------------------------#

done

#--------------------------------------------------------------------------------#

echo
echo "Get 1d plots, cutflows, and pie charts:"
for channel in "${channels[@]}"
do
    # Args to getPlots.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}, lastcut
    for cut in "${cuts[@]}"
    do
        python getPlots.py all save "$channel" "$cut"
    done

    # Args to plotCutflows.py: testMode {test, all}, displayMode {show, save}, 
    # channel {mumu, elel, muel}
    python plotCutflows.py all save "$channel"
done

#--------------------------------------------------------------------------------#
