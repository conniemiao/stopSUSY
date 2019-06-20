# calls plotSusyBkgd+Sig.py for all the variables and for all the possible
# combinations of settings

variables=("lep1_pt" "lep1_eta" "lep1_relIso" "lep2_pt" "lep2_eta" "lep2_relIso" \
    "njets" "jet_pt" "jet_eta" "nbtag" "deltaR_lep1_jet" "deltaR_lep2_jet" \
    "mtlep1" "mtlep2" "met_pt" "met_phi")

for var in "$variables"
do
    python plotSusyBkgd+Sig.py 0 0 0 0 "$var"
    python plotSusyBkgd+Sig.py 0 1 0 0 "$var"
    python plotSusyBkgd+Sig.py 0 0 1 0 "$var"
    python plotSusyBkgd+Sig.py 0 1 1 0 "$var"
    python plotSusyBkgd+Sig.py 0 0 1 1 "$var"
    python plotSusyBkgd+Sig.py 0 1 1 1 "$var"
done
