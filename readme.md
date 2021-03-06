#ReadMe

## SelEvents.cc

`.L SelEvents`

`SelEvents("NtupleName")`

produces a reduced ntuple (NtupleName_sel.root) with the main old branches plus:
* corrMet: corrected t1pfmet 
* commonSel: it's == 1 if (nEle < 2 && nMuons == 0 && metFilters == 1 && DeltaPhi(gg,Met) > 2.1 && maxDeltaPhi(jet,Met) < 2.7 && minDeltaPhi(jet,Met) > 0.5)


## OptScan.C
`.x OptScan.C(<Zmass>,<A0mass>)`

(takes as input the reduced ntuples from SelEvents.cc)


4D scan on pt1/mgg, pt2/mgg, corrMet and ptgg variables, on the events passing the common selection

the range of the scan for each variable needs to be changed manually in the code through the initial values and the required number of steps 


## Significance.C
`.x Significance.C(<Zmass>,<A0mass>,<pt1_cut>,<pt2_cut>,<met_cut>,<ptgg_cut>)`

(using reduced ntuples from SelEvents.cc)

S events, B events and significance applying the cuts

## Plots2DSignificance.C
`.x Plots2DSignificance.C(<Zmass>,<A0mass>,<pt1_cut>,<pt2_cut>,<met_cut>,<ptgg_cut>)`

(the input files for this macro are the original ntuples, *not the reduced ones*)

produces 2D plots of the significance of the threesholds; the common selection (nEle < 2, nMuons == 0, metFilters == 1, DeltaPhi(gg,Met) > 2.1, maxDeltaPhi(jet,Met) < 2.7, minDeltaPhi(jet,Met) > 0.5) is applied, the MET is not corrected
