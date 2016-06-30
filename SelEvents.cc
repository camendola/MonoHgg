#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;

void SelEvents(const char* filename) {

  cout << "File " << filename << endl;  
  
  TFile *fileOrig = 0;
  TTree *treeOrig = 0;
  
  fileOrig = TFile::Open(Form("%s.root",filename));
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("DiPhotonTree");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }
  
  
  fileOrig->cd();
  if (!treeOrig) {
    cout << "Tree DiPhotonTree not existing !" << endl; 
    return;    
  }
  
  int nentriesOrig = treeOrig->GetEntries();  
  
  TFile *fileNew = TFile::Open(Form("%s_sel.root",filename),"recreate");
  fileNew->ls();  
  fileNew->cd(); 
  
  
  //  vector<TTree*> trees; 
  //  trees.push_back(treeNew);
  
  // original tree leaves
  Int_t           sampleID;
  Int_t           run;
  Int_t           event;
  Int_t           lumi;
  Int_t           nvtx;
  Float_t         weight; 
  Float_t         pu_weight;
  Float_t	  t1pfmet;
  Float_t	  t1pfmetPhi;
  Float_t	  t1pfmetSumEt;
  Float_t         ptgg;
  Float_t         mgg;
  Float_t         pt1;
  Float_t         eta1;
  Float_t         phi1;
  Float_t         ptOverM1;
  Float_t         pt2;
  Float_t         eta2;
  Float_t         phi2;
  Float_t         ptOverM2;
  Float_t         mggNominal;
  Float_t         mggGen;
  Float_t         ptJetLead;
  Float_t         etaJetLead;
  Float_t         phiJetLead;
  Float_t         massJetLead;
  Float_t         ptJetSubLead;
  Float_t         etaJetSubLead;
  Float_t         phiJetSubLead;
  Float_t         massJetSubLead;
  Float_t         ptJet3;
  Float_t         etaJet3;
  Float_t         phiJet3;
  Float_t         massJet3;
  Float_t         ptJet4;
  Float_t         etaJet4;
  Float_t         phiJet4;
  Float_t         massJet4;
  Int_t           nEle;   
  Int_t           nMuons;    
  Int_t           eleveto1;    
  Int_t           eleveto2;    
  Int_t           genmatch1;    
  Int_t           genmatch2;    
  Int_t           metF_GV; 
  Int_t           metF_HBHENoise; 
  Int_t           metF_HBHENoiseIso;    
  Int_t           metF_CSC;   
  Int_t           metF_eeBadSC;   
  Int_t           metF_MuonBadTrack;   
  Int_t           metF_HadronTrackRes;   

  // List of branches - original tree
  TBranch        *b_sampleID; 
  TBranch        *b_run; 
  TBranch        *b_event;
  TBranch        *b_lumi;
  TBranch        *b_nvtx;
  TBranch        *b_weight;
  TBranch        *b_pu_weight;
  TBranch	 *b_t1pfmet;
  TBranch	 *b_t1pfmetPhi;
  TBranch	 *b_t1pfmetSumEt;
  TBranch        *b_ptgg;
  TBranch        *b_mgg; 
  TBranch        *b_pt1; 
  TBranch        *b_eta1; 
  TBranch        *b_phi1; 
  TBranch        *b_ptOverM1;
  TBranch        *b_pt2; 
  TBranch        *b_eta2; 
  TBranch        *b_phi2; 
  TBranch        *b_ptOverM2;
  TBranch        *b_mggNominal;
  TBranch        *b_mggGen;
  TBranch        *b_ptJetLead;   //!
  TBranch        *b_etaJetLead;   //!
  TBranch        *b_phiJetLead;   //!
  TBranch        *b_massJetLead;   //!
  TBranch        *b_ptJetSubLead;   //!
  TBranch        *b_etaJetSubLead;   //!
  TBranch        *b_phiJetSubLead;   //!
  TBranch        *b_massJetSubLead;   //!
  TBranch        *b_ptJet3;   //!
  TBranch        *b_etaJet3;   //!
  TBranch        *b_phiJet3;   //!
  TBranch        *b_massJet3;   //!
  TBranch        *b_ptJet4;   //!
  TBranch        *b_etaJet4;   //!
  TBranch        *b_phiJet4;   //!
  TBranch        *b_massJet4;   //!
  TBranch        *b_nEle;   //!  
  TBranch        *b_nMuons;   //!  
  TBranch        *b_eleveto1;   //!  
  TBranch        *b_eleveto2;   //!  
  TBranch        *b_genmatch1;   //!  
  TBranch        *b_genmatch2;   //!  
  TBranch        *b_metF_GV;   //!  
  TBranch        *b_metF_HBHENoise;   //!
  TBranch        *b_metF_HBHENoiseIso;   //! 
  TBranch        *b_metF_CSC;   //!  
  TBranch        *b_metF_eeBadSC;   //! 
  TBranch        *b_metF_MuonBadTrack;   //!  
  TBranch        *b_metF_HadronTrackRes;   //!   
  

  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
  treeOrig->SetBranchAddress("run", &run, &b_run);
  treeOrig->SetBranchAddress("event", &event, &b_event);
  treeOrig->SetBranchAddress("lumi", &lumi, &b_lumi);
  treeOrig->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  treeOrig->SetBranchAddress("weight", &weight, &b_weight);
  treeOrig->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  treeOrig->SetBranchAddress("t1pfmet", &t1pfmet, &b_t1pfmet);
  treeOrig->SetBranchAddress("t1pfmetPhi", &t1pfmetPhi, &b_t1pfmetPhi);
  treeOrig->SetBranchAddress("t1pfmetSumEt", &t1pfmetSumEt, &b_t1pfmetSumEt);
  treeOrig->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
  treeOrig->SetBranchAddress("mgg", &mgg, &b_mgg);
  treeOrig->SetBranchAddress("pt1", &pt1, &b_pt1);
  treeOrig->SetBranchAddress("eta1", &eta1, &b_eta1);
  treeOrig->SetBranchAddress("phi1", &phi1, &b_phi1);
  treeOrig->SetBranchAddress("ptOverM1", &ptOverM1, &b_ptOverM1);
  treeOrig->SetBranchAddress("pt2", &pt2, &b_pt2);
  treeOrig->SetBranchAddress("eta2", &eta2, &b_eta2);
  treeOrig->SetBranchAddress("phi2", &phi2, &b_phi2);
  treeOrig->SetBranchAddress("ptOverM2", &ptOverM2, &b_ptOverM2);
  treeOrig->SetBranchAddress("mggNominal", &mggNominal, &b_mggNominal);
  treeOrig->SetBranchAddress("mggGen", &mggGen, &b_mggGen);
  treeOrig->SetBranchAddress("ptJetLead", &ptJetLead, &b_ptJetLead);
  treeOrig->SetBranchAddress("etaJetLead", &etaJetLead, &b_etaJetLead);
  treeOrig->SetBranchAddress("phiJetLead", &phiJetLead, &b_phiJetLead);
  treeOrig->SetBranchAddress("massJetLead", &massJetLead, &b_massJetLead);
  treeOrig->SetBranchAddress("ptJetSubLead", &ptJetSubLead, &b_ptJetSubLead);
  treeOrig->SetBranchAddress("etaJetSubLead", &etaJetSubLead, &b_etaJetSubLead);
  treeOrig->SetBranchAddress("phiJetSubLead", &phiJetSubLead, &b_phiJetSubLead);
  treeOrig->SetBranchAddress("massJetSubLead", &massJetSubLead, &b_massJetSubLead);
  treeOrig->SetBranchAddress("ptJet3", &ptJet3, &b_ptJet3);
  treeOrig->SetBranchAddress("etaJet3", &etaJet3, &b_etaJet3);
  treeOrig->SetBranchAddress("phiJet3", &phiJet3, &b_phiJet3);
  treeOrig->SetBranchAddress("massJet3", &massJet3, &b_massJet3);
  treeOrig->SetBranchAddress("ptJet4", &ptJet4, &b_ptJet4);
  treeOrig->SetBranchAddress("etaJet4", &etaJet4, &b_etaJet4);
  treeOrig->SetBranchAddress("phiJet4", &phiJet4, &b_phiJet4);
  treeOrig->SetBranchAddress("massJet4", &massJet4, &b_massJet4);
  treeOrig->SetBranchAddress("nEle", &nEle, &b_nEle);      
  treeOrig->SetBranchAddress("nMuons", &nMuons, &b_nMuons);  
  treeOrig->SetBranchAddress("eleveto1", &eleveto1, &b_eleveto1);  
  treeOrig->SetBranchAddress("eleveto2", &eleveto2, &b_eleveto2);  
  treeOrig->SetBranchAddress("genmatch1", &genmatch1, &b_genmatch1);  
  treeOrig->SetBranchAddress("genmatch2", &genmatch2, &b_genmatch2);  
  treeOrig->SetBranchAddress("metF_GV", &metF_GV, &b_metF_GV);    
  treeOrig->SetBranchAddress("metF_HBHENoise", &metF_HBHENoise, &b_metF_HBHENoise);   
  treeOrig->SetBranchAddress("metF_HBHENoiseIso", &metF_HBHENoiseIso, &b_metF_HBHENoiseIso);  
  treeOrig->SetBranchAddress("metF_CSC", &metF_CSC, &b_metF_CSC); 
  treeOrig->SetBranchAddress("metF_eeBadSC", &metF_eeBadSC, &b_metF_eeBadSC);      
  treeOrig->SetBranchAddress("metF_MuonBadTrack", &metF_MuonBadTrack, &b_metF_MuonBadTrack); 
  treeOrig->SetBranchAddress("metF_HadronTrackRes", &metF_HadronTrackRes, &b_metF_HadronTrackRes); 
  
 
 // new variables to be added
  Int_t commonSel;
  Float_t corrMet;
  Float_t deltaPhiMaxJet;
  Float_t deltaPhiMinJet;
  Float_t deltaPhiMetGG;
  Int_t metFilters;
 
  TTree *Tnew = new TTree("DiPhotonTree","tree with 2 photon selection");
  Tnew->SetAutoSave(-99999999999);     
  Tnew->SetAutoFlush(-99999999999);   
  Tnew->Branch("commonSel", &commonSel, "commonSel/I");
  Tnew->Branch("corrMet", &corrMet, "corrMet/F");
  Tnew->Branch("deltaPhiMaxJet", &deltaPhiMaxJet, "deltaPhiMaxJet/F");
  Tnew->Branch("deltaPhiMinJet", &deltaPhiMinJet, "deltaPhiMinJet/F");
  Tnew->Branch("deltaPhiMetGG", &deltaPhiMetGG, "deltaPhiMetGG/F");
  Tnew->Branch("metFilters",&metFilters,"metFilters/I");
  
  // Copy branches
  Tnew->Branch("run", &run, "run/I");
  Tnew->Branch("event", &event, "event/I");
  Tnew->Branch("lumi", &lumi, "lumi/I");
  Tnew->Branch("nvtx", &nvtx, "nvtx/I");
  Tnew->Branch("weight", &weight, "weight/F");
  Tnew->Branch("pu_weight", &pu_weight, "pu_weight/F");
  Tnew->Branch("t1pfmet", &t1pfmet, "t1pfmet/F");
  Tnew->Branch("t1pfmetPhi", &t1pfmetPhi, "t1pfmetPhi/F");
  Tnew->Branch("t1pfmetSumEt", &t1pfmetSumEt, "t1pfmetSumEt/F");
  Tnew->Branch("ptgg", &ptgg, "ptgg/F");
  Tnew->Branch("mgg", &mgg, "mgg/F");
  Tnew->Branch("pt1", &pt1, "pt1/F");
  Tnew->Branch("ptOverM1", &ptOverM1, "ptOverM1/F");
  Tnew->Branch("pt2", &pt2, "pt2/F");
  Tnew->Branch("ptOverM2", &ptOverM2, "ptOverM2/F");
  Tnew->Branch("nEle", &nEle, "nEle/I");      
  Tnew->Branch("nMuons", &nMuons, "nMuons/I");  
  Tnew->Branch("eleveto1", &eleveto1, "eleveto1/I");  
  Tnew->Branch("eleveto2", &eleveto2, "eleveto2/I");  
  Tnew->Branch("genmatch1", &genmatch1, "genmatch1/I");  
  Tnew->Branch("genmatch2", &genmatch2, "genmatch2/I");  
 
  

  
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0){
      std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    }
    treeOrig->GetEntry(i);
    
    // met corrections
    float corrMetPx, corrMetPy, corrMetE;
    if (sampleID>0 && sampleID<10000) { //MC
      corrMetPx = t1pfmet*cos(t1pfmetPhi) -(-2.32806 + -0.0178001 * t1pfmetSumEt);
      corrMetPy = t1pfmet*sin(t1pfmetPhi) -(0.454459 + 0.00144021 * t1pfmetSumEt);
    
    } else { //Data   
      corrMetPx = t1pfmet*cos(t1pfmetPhi) -(-2.35978 + -0.0135418 * t1pfmetSumEt);
      corrMetPy = t1pfmet*sin(t1pfmetPhi) -(1.01529 + 0.0198305 * t1pfmetSumEt);    
  
    }
    corrMetE = sqrt(corrMetPx*corrMetPx + corrMetPy*corrMetPy);
    TLorentzVector correctedMetTLV;
    correctedMetTLV.SetPxPyPzE(corrMetPx,corrMetPy,0,corrMetE);
    corrMet = correctedMetTLV.Pt();  
    // deltaPhi jets - met    
    float dPhiMetJ1 = -1000.;
    float dPhiMetJ2 = -1000.;
    float dPhiMetJ3 = -1000.;
    float dPhiMetJ4 = -1000.;
    if (ptJetLead>50) {
      TLorentzVector tlvJet1;
      tlvJet1.SetPtEtaPhiM(ptJetLead,etaJetLead,phiJetLead,massJetLead);
      dPhiMetJ1 = correctedMetTLV.DeltaPhi(tlvJet1);
    }
    if (ptJetSubLead>50) {
      TLorentzVector tlvJet2;
      tlvJet2.SetPtEtaPhiM(ptJetSubLead,etaJetSubLead,phiJetSubLead,massJetSubLead);
      dPhiMetJ2 = correctedMetTLV.DeltaPhi(tlvJet2);
    }
    if (ptJet3>50) {
      TLorentzVector tlvJet3;
      tlvJet3.SetPtEtaPhiM(ptJet3,etaJet3,phiJet3,massJet3);
      dPhiMetJ3 = correctedMetTLV.DeltaPhi(tlvJet3);
    }
    if (ptJet4>50) {
      TLorentzVector tlvJet4;
      tlvJet4.SetPtEtaPhiM(ptJet4,etaJet4,phiJet4,massJet4);
      dPhiMetJ4 = correctedMetTLV.DeltaPhi(tlvJet4);
    }
    
    float maxDPhiMetJ = 0.;
    if (dPhiMetJ1>=-4 && fabs(dPhiMetJ1)>maxDPhiMetJ) maxDPhiMetJ = fabs(dPhiMetJ1);
    if (dPhiMetJ2>=-4 && fabs(dPhiMetJ2)>maxDPhiMetJ) maxDPhiMetJ = fabs(dPhiMetJ2);
    if (dPhiMetJ3>=-4 && fabs(dPhiMetJ3)>maxDPhiMetJ) maxDPhiMetJ = fabs(dPhiMetJ3);
    if (dPhiMetJ4>=-4 && fabs(dPhiMetJ4)>maxDPhiMetJ) maxDPhiMetJ = fabs(dPhiMetJ4);
    deltaPhiMaxJet = maxDPhiMetJ;
    
    float minDPhiMetJ = 500000.;
    if (dPhiMetJ1>=-4 && fabs(dPhiMetJ1)<minDPhiMetJ) minDPhiMetJ = fabs(dPhiMetJ1);
    if (dPhiMetJ2>=-4 && fabs(dPhiMetJ2)<minDPhiMetJ) minDPhiMetJ = fabs(dPhiMetJ2);
    if (dPhiMetJ3>=-4 && fabs(dPhiMetJ3)<minDPhiMetJ) minDPhiMetJ = fabs(dPhiMetJ3);
    if (dPhiMetJ4>=-4 && fabs(dPhiMetJ4)<minDPhiMetJ) minDPhiMetJ = fabs(dPhiMetJ4);
    deltaPhiMinJet = minDPhiMetJ;
    
    //check
    /*if (i%10000 == 0){
      std::cout << "Jet1:"<<endl;
      std::cout << "ptJetLead:"<<ptJetLead<<"; phiJetLead: "<<phiJetLead<<"; dPhiMetJ1 = "<<dPhiMetJ1<<endl;
      std::cout << "Jet2:"<<endl;
      std::cout << "ptJetSubLead:"<<ptJetSubLead<<"; phiJetSubLead: "<<phiJetSubLead<<"; dPhiMetJ2 = "<<dPhiMetJ2<<endl;
      std::cout << "Jet3:"<<endl;
      std::cout << "ptJet3:"<<ptJet3<<"; phiJet3: "<<phiJet3<<"; dPhiMetJ3 = "<<dPhiMetJ3<<endl;
      std::cout << "Jet4:"<<endl;
      std::cout << "ptJet4:"<<ptJet4<<"; phiJet4: "<<phiJet4<<"; dPhiMetJ4 = "<<dPhiMetJ4<<endl;
      std::cout << "minDPhiMetJ = "<<minDPhiMetJ<<"; maxDPhiMetJ = "<<maxDPhiMetJ<<endl;
      }*/

    // deltaPhi gg - met
    TLorentzVector g1TLV, g2TLV, ggTLV;
    g1TLV.SetPtEtaPhiM(pt1,eta1,phi1,0);
    g2TLV.SetPtEtaPhiM(pt2,eta2,phi2,0);
    ggTLV = g1TLV + g2TLV; 
    deltaPhiMetGG = fabs(ggTLV.DeltaPhi(correctedMetTLV));

    // met filters - data only
    metFilters = 1;
    if (sampleID>=10000) metFilters = metF_GV && metF_HBHENoise && metF_HBHENoiseIso && metF_CSC && metF_eeBadSC && metF_MuonBadTrack && metF_HadronTrackRes;
    
    // check if full selection is passed - including met filters on data
    commonSel = ( (nEle<2) && (nMuons==0) && (maxDPhiMetJ<2.7) && (minDPhiMetJ>0.5) && (deltaPhiMetGG>2.1) && metFilters==1 );

    Tnew->Fill();
  }

  // new info
  fileNew->cd();
  Tnew->Write();
  fileNew->Close();

  fileOrig->cd();
  fileOrig->Close();
  
}




