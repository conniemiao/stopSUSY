# plots all the variables in the ntuples and the cutflows for all the variables
# and for all the possible combinations of settings

# variables=("lep1_pt" "lep1_eta" "lep1_relIso" "lep2_pt" \
#     "njets" "jet_pt" "jet_eta" "nbtag" "deltaR_lep1_jet" \
#     "mtlep1" "mtlep2" "met_pt")
variables=("deltaR_lep1_jet")

python plotSusyBkgd+Sig.py 0 0 0 0 "${variables[@]}"
echo
python plotSusyBkgd+Sig.py 0 1 0 0 "${variables[@]}"
echo
python plotSusyBkgd+Sig.py 0 0 1 0 "${variables[@]}"
echo
python plotSusyBkgd+Sig.py 0 1 1 0 "${variables[@]}"
echo
python plotSusyBkgd+Sig.py 0 0 1 1 "${variables[@]}"
echo
python plotSusyBkgd+Sig.py 0 1 1 1 "${variables[@]}"
echo
# python plotCutflowsBkgd+Sig.py 0 1 0 0
# echo
# python plotCutflowsBkgd+Sig.py 0 1 1 0
# echo
# python plotCutflowsBkgd+Sig.py 0 1 1 1
