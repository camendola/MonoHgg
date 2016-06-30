#ReadMe

## SelEvents.cc

`.L SelEvents`

`SelEvents("NtupleName")`

produces a reduced ntuple (NtupleName_sel.root) with the main old branches plus:
* corrMet (corrected t1pfmet)
* commonSel (== 1 if nEle < 2, nMuons == 0, metFilters == 1, DeltaPhi(gg,Met) > 2.1, maxDeltaPhi(jet,Met) < 2.7, minDeltaPhi(jet,Met) > 0.5)


##OptScan.C
`.x OptScan(<Zmass>,<A0mass>)`

(using reduced ntuples from SelEvents.cc)
4D scan on pt1/mgg, pt2/mgg, corrMet and ptgg
