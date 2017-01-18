#include <TH1.h>
#include <THStack.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TText.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "FDelta.C"
#include "CMS_lumi.C"
double MAX(double ptJetLead, double phiJetLead, double ptJetSubLead, double phiJetSubLead,double ptJet3, double phiJet3,double ptJet4, double phiJet4,double t1pfmetPhi){
  double Max = 0.;
  double ptJetCut = 50.;
  double deltaLead = -10;
  if(ptJetLead>ptJetCut){
deltaLead = phiJetLead - t1pfmetPhi;

  if(deltaLead>3.14){
    deltaLead-=6.28;
  }
  if(deltaLead<-3.14){
    deltaLead+=6.28;
  }
  deltaLead =fabs(deltaLead);
  }
  double deltaSub = -10;
  if(ptJetSubLead>ptJetCut){
    deltaSub = phiJetSubLead - t1pfmetPhi;
    if(deltaSub>3.14){
      deltaSub-=6.28;
  }
    if(deltaSub<-3.14){
      deltaSub+=6.28;
    }
    deltaSub =fabs(deltaSub);
  }
  double delta3 = -10;
if(ptJet3>ptJetCut){
  delta3 = phiJet3 - t1pfmetPhi;
  if(delta3>3.14){
    delta3-=6.28;
  }
  if(delta3<-3.14){
    delta3+=6.28;
  }
  delta3 =fabs(delta3);
 }
 double delta4 = -10;
 if(ptJet4>ptJetCut){
   delta4= phiJet4 - t1pfmetPhi;
  if(delta4>3.14){
    delta4-=6.28;
  }
  if(delta4<-3.14){
    delta4+=6.28;
  }
  delta4 =fabs(delta4);
 }
  
 float maxDPhiMetJ = 0.;
 if (deltaLead>=-4 && fabs(deltaLead)>maxDPhiMetJ) maxDPhiMetJ = fabs(deltaLead);
 if (deltaSub>=-4 && fabs(deltaSub)>maxDPhiMetJ) maxDPhiMetJ = fabs(deltaSub);
 if (delta3>=-4 && fabs(delta3)>maxDPhiMetJ) maxDPhiMetJ = fabs(delta3);
 if (delta4>=-4 && fabs(delta4)>maxDPhiMetJ) maxDPhiMetJ = fabs(delta4);
 Max = maxDPhiMetJ;
  return Max;

}



double myMIN(double ptJetLead, double phiJetLead, double ptJetSubLead, double phiJetSubLead,double ptJet3, double phiJet3,double  ptJet4, double phiJet4,double t1pfmetPhi){
  double Min =3.14;
  double ptJetCut = 50.;

  double deltaLead = -10;
  if(ptJetLead>ptJetCut){
    deltaLead = phiJetLead - t1pfmetPhi;

    if(deltaLead>3.14){
      deltaLead-=6.28;
  }
    if(deltaLead<-3.14){
      deltaLead+=6.28;
    }
    deltaLead =fabs(deltaLead);
  }
  double deltaSub = -10;
  if(ptJetSubLead>ptJetCut){
    deltaSub = phiJetSubLead - t1pfmetPhi;
    if(deltaSub>3.14){
      deltaSub-=6.28;
    }
    if(deltaSub<-3.14){
      deltaSub+=6.28;
    }
    deltaSub =fabs(deltaSub);
  }
  double delta3 = -10;
  if(ptJet3>ptJetCut){
    delta3 = phiJet3 - t1pfmetPhi;
    if(delta3>3.14){
      delta3-=6.28;
    }
    if(delta3<-3.14){
      delta3+=6.28;
    }
    delta3 =fabs(delta3);
  }
  double delta4 = -10;
  if(ptJet4>ptJetCut){
    delta4= phiJet4 - t1pfmetPhi;
    if(delta4>3.14){
      delta4-=6.28;
    }
    if(delta4<-3.14){
      delta4+=6.28;
    }
    delta4 =fabs(delta4);
  }
  
  float minDPhiMetJ = 500000.;
  if (deltaLead>=-4 && fabs(deltaLead)<minDPhiMetJ) minDPhiMetJ = fabs(deltaLead);
  if (deltaSub>=-4 && fabs(deltaSub)<minDPhiMetJ) minDPhiMetJ = fabs(deltaSub);
  if (delta3>=-4 && fabs(delta3)<minDPhiMetJ) minDPhiMetJ = fabs(delta3);
  if (delta4>=-4 && fabs(delta4)<minDPhiMetJ) minDPhiMetJ = fabs(delta4);
  Min = minDPhiMetJ;
  return Min;
  
}



void Plots2DSignificance(int sigMass, double cut_pt1IN, double cut_pt2IN, double cut_metIN,double cut_ptggIN){

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);



  /*
  TFile *sig = TFile::Open(Form("./50ns_betaV4/NewWeightDMHtoGG_M%d.root",sigMass),"READ");
  TFile *bkg1 =  TFile::Open("./50ns_betaV4/NewWeightDiPhoton.root","READ"); 
  TFile *bkg2 =  TFile::Open("./50ns_betaV4/NewWeightDYJetsToLL.root","READ");
  TFile *bkg3 =  TFile::Open("./50ns_betaV4/NewWeightGJets.root","READ");  
  TFile *bkg4 =  TFile::Open("./50ns_betaV4/NewWeightGluGluHToGG.root","READ");
  TFile *bkg5 =  TFile::Open("./50ns_betaV4/NewWeightQCD.root","READ");
  TFile *bkg6 =  TFile::Open("./50ns_betaV4/NewWeightVH.root","READ");
  */  
  TFile *data =  TFile::Open("./25ns_2246inv_v3/DoubleEG.root","READ"); 
  TFile *sig = TFile::Open(Form("./25ns_2246inv_v3/2HDM_mZP%d.root",sigMass),"READ");
  if(sig == NULL) cout<<"sig fails"<<endl;
  TFile *bkg1 =  TFile::Open("./25ns_2246inv_v3/DiPhoton.root","READ"); 
  if(bkg1 == NULL) cout<<"bkg1 fails"<<endl;
  TFile *bkg2 =  TFile::Open("./25ns_2246inv_v3/DYJetsToLL.root","READ");
  if(bkg2 == NULL) cout<<"bkg2 fails"<<endl;
  TFile *bkg3 =  TFile::Open("./25ns_2246inv_v3/GJets.root","READ");  
  if(bkg3 == NULL) cout<<"bkg3 fails"<<endl;
  TFile *bkg4 =  TFile::Open("./25ns_2246inv_v3/GluGluHToGG.root","READ");
  if(bkg4 == NULL) cout<<"bkg4 fails"<<endl;
  TFile *bkg5 =  TFile::Open("./25ns_2246inv_v3/QCD.root","READ");
  if(bkg5 == NULL) cout<<"bkg5 fails"<<endl;
  TFile *bkg6 =  TFile::Open("./25ns_2246inv_v3/VH.root","READ");
  if(bkg6 == NULL) cout<<"bkg6 fails"<<endl;
  TFile *bkg7 =  TFile::Open("./25ns_2246inv_v3/ttHJetToGG.root","READ");
  TFile *bkg8 =  TFile::Open("./25ns_2246inv_v3/VBFHToGG.root","READ");
  TFile *bkg9 =  TFile::Open("./25ns_2246inv_v3/TGJets.root","READ");
  TFile *bkg10 =  TFile::Open("./25ns_2246inv_v3/TTGJets.root","READ");
  TFile *bkg11 =  TFile::Open("./25ns_2246inv_v3/WGToLNuG.root","READ");
  TFile *bkg12 =  TFile::Open("./25ns_2246inv_v3/ZGTo2LG.root","READ");


 
  TTree *tree_data = (TTree*) data->Get("DiPhotonTree"); 
  TTree *tree_sig = (TTree*) sig->Get("DiPhotonTree"); 
  TTree *tree_bkg1 = (TTree*) bkg1->Get("DiPhotonTree");
  TTree *tree_bkg2 = (TTree*) bkg2->Get("DiPhotonTree");
  TTree *tree_bkg3 = (TTree*) bkg3->Get("DiPhotonTree");
  TTree *tree_bkg4 = (TTree*) bkg4->Get("DiPhotonTree");
  TTree *tree_bkg5 = (TTree*) bkg5->Get("DiPhotonTree");
  TTree *tree_bkg6 = (TTree*) bkg6->Get("DiPhotonTree");
  TTree *tree_bkg7 = (TTree*) bkg7->Get("DiPhotonTree");
  TTree *tree_bkg8 = (TTree*) bkg8->Get("DiPhotonTree");
  TTree *tree_bkg9 = (TTree*) bkg9->Get("DiPhotonTree");
  TTree *tree_bkg10 = (TTree*) bkg10->Get("DiPhotonTree");
  TTree *tree_bkg11 = (TTree*) bkg11->Get("DiPhotonTree");
  TTree *tree_bkg12 = (TTree*) bkg12->Get("DiPhotonTree");


  TCut Cut_pt1;
  TCut Cut_pt2;
  TCut Cut_met;
  TCut Cut_ptgg;
  TCut mggmax = "mgg<130";
  TCut mggmin = "mgg>120";
  TCut eveto1 = "eleveto1 == 1";
  TCut eveto2 = "eleveto2 == 1";
  TCut eveto = eveto1 && eveto2;
  TCut genmatch = "((genmatch1==1 && genmatch2==0)||(genmatch1==0 && genmatch2==1)||(genmatch1==0 && genmatch2==0))";  
  
  TCut DPHIcut = "FDelta(pt1,eta1,phi1,0.,pt2,eta2,phi2,0.,t1pfmetPhi)>2.1";
  TCut maxDPHIJetcut = "MAX(ptJetLead,phiJetLead,ptJetSubLead,phiJetSubLead,ptJet3,phiJet3,ptJet4,phiJet4,t1pfmetPhi)<2.7";
  TCut minDPHIJetcut = "myMIN(ptJetLead,phiJetLead,ptJetSubLead,phiJetSubLead,ptJet3,phiJet3,ptJet4,phiJet4,t1pfmetPhi)>0.5";
  TCut nMuonCut = "nMuons == 0";
  TCut nEleCut = "nEle <= 1";  
  
  TCanvas *canvas = new TCanvas("canvas","",600,400);   
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  double cut_pt1 = 0;
  double cut_pt2 = 0;
  double cut_met = 0;
  double cut_ptgg = 0;
  
  
  double Stot_pt1 = 0;
  double Stot_pt2 = 0;
  double Stot_met = 0;
  double Stot_ptgg = 0;
  double S = 0;
  double F1 = 0;
  double F2 = 0;
  double F3 = 0;
  double F4 = 0;
  double F5 = 0;
  double F6 = 0;
  double F7 = 0;
  double F8 = 0;
  double F9 = 0;
  double F10 = 0;
  double F11 = 0;
  double F12 = 0;
  double B = 0;
  double Eff = 0;
  
  double Sig3 = 0;
  double Sig5 = 0;
  int NstepsPt1 = 20; 
  int NstepsPt2 = 20; 
  int NstepsMET = 20; 
  int NstepsPtgg = 20; 
  double temp = 0;
  
  double firstCutPt1 = cut_pt1IN  - 0.05*NstepsPt1/2; 
  double firstCutPt2 = 1./4.;
  double firstCutMET = cut_metIN - 5*NstepsMET/2;
  double firstCutPtgg = cut_ptggIN - 5*NstepsPtgg/2;
  
  double crossSec = 1.; ////
  double N = 0.;
  
  
  
  cout<<"############Segnale M"<<sigMass<<endl;
  cout<<"120 < mgg < 130"<<endl;  
  ios_base::fmtflags oldflags = std::cout.flags();
  TH2F *hpt1VSpt2 = new TH2F("hpt1VSpt2","",NstepsPt1,firstCutPt1,firstCutPt1+0.05*NstepsPt1,NstepsPt2,firstCutPt2,firstCutPt2+0.05*NstepsPt2); 
  TH2F *hpt1VSmet = new TH2F("hpt1VSmet","",NstepsPt1,firstCutPt1,firstCutPt1+0.05*NstepsPt1,NstepsMET,firstCutMET,firstCutMET+5*NstepsMET);
  TH2F *hptggVSmet = new TH2F("hptggVSmet","",NstepsPtgg,firstCutPtgg,firstCutPtgg+5*NstepsPtgg,NstepsMET,firstCutMET,firstCutMET+5*NstepsMET);
  TH2F *hmetVSpt2 = new TH2F("metVSpt2","",NstepsMET,firstCutMET,firstCutMET+5*NstepsMET,NstepsPt2,firstCutPt2,firstCutPt2+0.05*NstepsPt2);
  TH2F *hptggVSpt2 = new TH2F("ptggVSpt2","",NstepsPtgg,firstCutPtgg,firstCutPtgg+5*NstepsPtgg,NstepsPt2,firstCutPt2,firstCutPt2+0.05*NstepsPt2);
  TH2F *hpt1VSptgg = new TH2F("hpt1VSptgg","",NstepsPt1,firstCutPt1,firstCutPt1+0.05*NstepsPt1,NstepsPtgg,firstCutPtgg,firstCutPtgg+5*NstepsPtgg);
  
  tree_sig->Draw("(pt1)>>pt1_tot1(30,0,1000)","weight"*(mggmin && mggmax && eveto));
  TH1F *pt1_tot1 =(TH1F*)gPad->GetPrimitive("pt1_tot1");
  pt1_tot1->Scale(crossSec); 
  Stot_pt1 = pt1_tot1->Integral();
  
  
  Cut_met = Form("t1pfmet>%lf",cut_metIN);
  Cut_ptgg = Form("ptgg>%lf",cut_ptggIN);
  for(i = 0; i<NstepsPt1;i++){
    cut_pt1 = firstCutPt1 + 0.05*i +0.001;
    Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1);   
    for(j = 0; j<NstepsPt2;j++){
      cut_pt2 = firstCutPt2 + 0.05*j +0.001;
      Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2);
      tree_data->Draw("(pt1)>>pt1_data(30,0,1000)","(mgg<115||mgg>135)"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
      TH1F *pt1_data =(TH1F*)gPad->GetPrimitive("pt1_data");
      N = pt1_data->GetEntries();
      Sig5= 0.;
      if(N>=0){
	
	
	tree_sig->Draw("(pt1)>>pt1_1(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto &&  maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_1 =(TH1F*)gPad->GetPrimitive("pt1_1"); 
	pt1_1->Scale(crossSec);      
	S = pt1_1->Integral();   
	tree_bkg1->Draw("(pt1)>>pt1_bkg1(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg1 =(TH1F*)gPad->GetPrimitive("pt1_bkg1");
	
	F1 = pt1_bkg1->Integral();   
	tree_bkg2->Draw("(pt1)>>pt1_bkg2(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg2 =(TH1F*)gPad->GetPrimitive("pt1_bkg2");
	
	F2 = pt1_bkg2->Integral();
	if(F2<0) F2 = 0;          
	tree_bkg3->Draw("(pt1)>>pt1_bkg3(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg3 =(TH1F*)gPad->GetPrimitive("pt1_bkg3");
	
	F3 = pt1_bkg3->Integral();  
	tree_bkg4->Draw("(pt1)>>pt1_bkg4(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg4 =(TH1F*)gPad->GetPrimitive("pt1_bkg4");
	
	F4 = pt1_bkg4->Integral();    
	tree_bkg5->Draw("(pt1)>>pt1_bkg5(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg5 =(TH1F*)gPad->GetPrimitive("pt1_bkg5");
	
	F5 = pt1_bkg5->Integral();    
	tree_bkg6->Draw("(pt1)>>pt1_bkg6(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut)); 
	TH1F *pt1_bkg6 =(TH1F*)gPad->GetPrimitive("pt1_bkg6");
	F6 = pt1_bkg6->Integral();  
	tree_bkg7->Draw("(pt1)>>pt1_bkg7(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg7 =(TH1F*)gPad->GetPrimitive("pt1_bkg7");
	F7 = pt1_bkg7->Integral();
	tree_bkg8->Draw("(pt1)>>pt1_bkg8(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                            
	TH1F *pt1_bkg8 =(TH1F*)gPad->GetPrimitive("pt1_bkg8");
	F8 = pt1_bkg8->Integral();
	tree_bkg9->Draw("(pt1)>>pt1_bkg9(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg9 =(TH1F*)gPad->GetPrimitive("pt1_bkg9");
	F9 = pt1_bkg9->Integral();
	if(F9<0) F9=0;
	tree_bkg10->Draw("(pt1)>>pt1_bkg10(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg10 =(TH1F*)gPad->GetPrimitive("pt1_bkg10");
	F10 = pt1_bkg10->Integral();
	tree_bkg11->Draw("(pt1)>>pt1_bkg11(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                       
	TH1F *pt1_bkg11 =(TH1F*)gPad->GetPrimitive("pt1_bkg11");
	F11 = pt1_bkg11->Integral();
	tree_bkg12->Draw("(pt1)>>pt1_bkg12(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg12 =(TH1F*)gPad->GetPrimitive("pt1_bkg12");
	F12 = pt1_bkg12->Integral();
	
	
	
	B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12;
	Eff = S/Stot_pt1;
	Sig3 = Eff/(1.5+sqrt(B));
	Sig5 = Eff/(2.5+sqrt(B));
      }
      hpt1VSpt2->Fill(cut_pt1,cut_pt2,Sig5);
    }
  }
  
  /////////////////////////////////////////////////////////
  Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2IN);
  Cut_ptgg = Form("ptgg>%lf",cut_ptggIN);
  for(i = 0; i<NstepsPt1;i++){
    cut_pt1 = firstCutPt1 + 0.05*i +0.001;
    Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1);
    for(j = 0; j<NstepsMET;j++){
      cut_met = firstCutMET + 5*j;
      Cut_met = Form("t1pfmet>%lf",cut_met);
      tree_data->Draw("(pt1)>>pt1_data(30,0,1000)","(mgg<115||mgg>135)"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
      TH1F *pt1_data =(TH1F*)gPad->GetPrimitive("pt1_data");
      N = pt1_data->GetEntries();
      Sig5= 0.;
      if(N>=0){
	
	
	tree_sig->Draw("(pt1)>>pt1_1(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_1 =(TH1F*)gPad->GetPrimitive("pt1_1"); 
	pt1_1->Scale(crossSec);   
	S = pt1_1->Integral();   
	tree_bkg1->Draw("(pt1)>>pt1_bkg1(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg1 =(TH1F*)gPad->GetPrimitive("pt1_bkg1");
	
	F1 = pt1_bkg1->Integral();   
	tree_bkg2->Draw("(pt1)>>pt1_bkg2(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg2 =(TH1F*)gPad->GetPrimitive("pt1_bkg2");
	
	F2 = pt1_bkg2->Integral();
	if(F2<0) F2 = 0;     
	tree_bkg3->Draw("(pt1)>>pt1_bkg3(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg3 =(TH1F*)gPad->GetPrimitive("pt1_bkg3");
	
	F3 = pt1_bkg3->Integral();  
	tree_bkg4->Draw("(pt1)>>pt1_bkg4(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg4 =(TH1F*)gPad->GetPrimitive("pt1_bkg4");
	
	F4 = pt1_bkg4->Integral();    
	tree_bkg5->Draw("(pt1)>>pt1_bkg5(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg5 =(TH1F*)gPad->GetPrimitive("pt1_bkg5");
	
	F5 = pt1_bkg5->Integral();   
 	tree_bkg6->Draw("(pt1)>>pt1_bkg6(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut)); 
	TH1F *pt1_bkg6 =(TH1F*)gPad->GetPrimitive("pt1_bkg6");
	F6 = pt1_bkg6->Integral();  
	tree_bkg7->Draw("(pt1)>>pt1_bkg7(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_met && Cut_pt2 && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut)); 
	TH1F *pt1_bkg7 =(TH1F*)gPad->GetPrimitive("pt1_bkg7");
	F7 = pt1_bkg7->Integral();
	tree_bkg8->Draw("(pt1)>>pt1_bkg8(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                            
	TH1F *pt1_bkg8 =(TH1F*)gPad->GetPrimitive("pt1_bkg8");
	F8 = pt1_bkg8->Integral();
	tree_bkg9->Draw("(pt1)>>pt1_bkg9(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg9 =(TH1F*)gPad->GetPrimitive("pt1_bkg9");
	F9 = pt1_bkg9->Integral();
	if(F9<0) F9 = 0;
	tree_bkg10->Draw("(pt1)>>pt1_bkg10(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg10 =(TH1F*)gPad->GetPrimitive("pt1_bkg10");
	F10 = pt1_bkg10->Integral();
	tree_bkg11->Draw("(pt1)>>pt1_bkg11(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                       
	TH1F *pt1_bkg11 =(TH1F*)gPad->GetPrimitive("pt1_bkg11");
	F11 = pt1_bkg11->Integral();
	tree_bkg12->Draw("(pt1)>>pt1_bkg12(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg12 =(TH1F*)gPad->GetPrimitive("pt1_bkg12");
	F12 = pt1_bkg12->Integral();
	

	
	B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12;
	
	
	

	
	Eff = S/Stot_pt1;
	Sig3 = Eff/(1.5+sqrt(B));
	Sig5 = Eff/(2.5+sqrt(B));
      }
	hpt1VSmet->Fill(cut_pt1,cut_met,Sig5);
    }
  }
  
  

  /////////////////////////////////////////////////////////
  Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1IN);
  Cut_ptgg = Form("ptgg>%lf",cut_ptggIN);
  for(i = 0; i<NstepsPt2;i++){
    cut_pt2 = firstCutPt2 + 0.05*i+0.001;
    Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2);
    for(j = 0; j<NstepsMET;j++){
      cut_met = firstCutMET + 5*j;
	Cut_met = Form("t1pfmet>%lf",cut_met);
	tree_data->Draw("(pt1)>>pt1_data(30,0,1000)","(mgg<115||mgg>135)"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_data =(TH1F*)gPad->GetPrimitive("pt1_data");
	N = pt1_data->GetEntries();
	Sig5= 0.;
	if(N>=0){
	  tree_sig->Draw("(pt1)>>pt1_1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	  TH1F *pt1_1 =(TH1F*)gPad->GetPrimitive("pt1_1"); 
	  pt1_1->Scale(crossSec);   	
	  S = pt1_1->Integral();   
	  tree_bkg1->Draw("(pt1)>>pt1_bkg1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	  TH1F *pt1_bkg1 =(TH1F*)gPad->GetPrimitive("pt1_bkg1");
	  
	  F1 = pt1_bkg1->Integral();   
	  tree_bkg2->Draw("(pt1)>>pt1_bkg2(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	  TH1F *pt1_bkg2 =(TH1F*)gPad->GetPrimitive("pt1_bkg2");
	
	  F2 = pt1_bkg2->Integral();
	  if(F2<0) F2 = 0;     
	  tree_bkg3->Draw("(pt1)>>pt1_bkg3(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	  TH1F *pt1_bkg3 =(TH1F*)gPad->GetPrimitive("pt1_bkg3");
	  
	  F3 = pt1_bkg3->Integral();  
	  tree_bkg4->Draw("(pt1)>>pt1_bkg4(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	  TH1F *pt1_bkg4 =(TH1F*)gPad->GetPrimitive("pt1_bkg4");
	  
	  F4 = pt1_bkg4->Integral();    
	  tree_bkg5->Draw("(pt1)>>pt1_bkg5(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	  TH1F *pt1_bkg5 =(TH1F*)gPad->GetPrimitive("pt1_bkg5");
	
	  F5 = pt1_bkg5->Integral();    
	  tree_bkg6->Draw("(pt1)>>pt1_bkg6(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut)); 
	  TH1F *pt1_bkg6 =(TH1F*)gPad->GetPrimitive("pt1_bkg6");
	  
	  F6 = pt1_bkg6->Integral();  
	  
	  tree_bkg7->Draw("(pt1)>>pt1_bkg7(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	  TH1F *pt1_bkg7 =(TH1F*)gPad->GetPrimitive("pt1_bkg7");
	  F7 = pt1_bkg7->Integral();
	  tree_bkg8->Draw("(pt1)>>pt1_bkg8(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                            
	  TH1F *pt1_bkg8 =(TH1F*)gPad->GetPrimitive("pt1_bkg8");
	  F8 = pt1_bkg8->Integral();
	  tree_bkg9->Draw("(pt1)>>pt1_bkg9(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	  TH1F *pt1_bkg9 =(TH1F*)gPad->GetPrimitive("pt1_bkg9");
	  F9 = pt1_bkg9->Integral();
	  if(F9<0) F9=0;
	  tree_bkg10->Draw("(pt1)>>pt1_bkg10(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	  TH1F *pt1_bkg10 =(TH1F*)gPad->GetPrimitive("pt1_bkg10");
	  F10 = pt1_bkg10->Integral();
	  tree_bkg11->Draw("(pt1)>>pt1_bkg11(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                       
	  TH1F *pt1_bkg11 =(TH1F*)gPad->GetPrimitive("pt1_bkg11");
	  F11 = pt1_bkg11->Integral();
	  tree_bkg12->Draw("(pt1)>>pt1_bkg12(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	  TH1F *pt1_bkg12 =(TH1F*)gPad->GetPrimitive("pt1_bkg12");
	  F12 = pt1_bkg12->Integral();
	  
	
	  
	  B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12;
	  
	  
	  Eff = S/Stot_pt1;
	  Sig3 = Eff/(1.5+sqrt(B));
	  Sig5 = Eff/(2.5+sqrt(B));
	}
	hmetVSpt2->Fill(cut_met,cut_pt2,Sig5);
    }
  }
  
  
  
  
  /////////////////////////////////////////////////////////
  Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2IN);
  Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1IN);
  for(i = 0; i<NstepsPtgg;i++){
    cut_ptgg = firstCutPtgg + 5*i;
    Cut_ptgg = Form("ptgg>%lf",cut_ptgg);
    for(j = 0; j<NstepsMET;j++){
      cut_met = firstCutMET + 5*j;
      Cut_met = Form("t1pfmet>%lf",cut_met);
      tree_data->Draw("(pt1)>>pt1_data(30,0,1000)","(mgg<115||mgg>135)"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
      TH1F *pt1_data =(TH1F*)gPad->GetPrimitive("pt1_data");
      N = pt1_data->GetEntries();
      Sig5= 0.;
      if(N>=0){
	tree_sig->Draw("(pt1)>>pt1_1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_1 =(TH1F*)gPad->GetPrimitive("pt1_1"); 
	pt1_1->Scale(crossSec);   	
	S = pt1_1->Integral();   
	tree_bkg1->Draw("(pt1)>>pt1_bkg1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg1 =(TH1F*)gPad->GetPrimitive("pt1_bkg1");
	
	F1 = pt1_bkg1->Integral();   
	tree_bkg2->Draw("(pt1)>>pt1_bkg2(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg2 =(TH1F*)gPad->GetPrimitive("pt1_bkg2");
	
	F2 = pt1_bkg2->Integral();
	if(F2<0) F2 = 0;     
	tree_bkg3->Draw("(pt1)>>pt1_bkg3(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg3 =(TH1F*)gPad->GetPrimitive("pt1_bkg3");
	
	F3 = pt1_bkg3->Integral();  
	tree_bkg4->Draw("(pt1)>>pt1_bkg4(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg4 =(TH1F*)gPad->GetPrimitive("pt1_bkg4");
	
	F4 = pt1_bkg4->Integral();    
	tree_bkg5->Draw("(pt1)>>pt1_bkg5(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg5 =(TH1F*)gPad->GetPrimitive("pt1_bkg5");
	
	F5 = pt1_bkg5->Integral();    
	tree_bkg6->Draw("(pt1)>>pt1_bkg6(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut)); 
	TH1F *pt1_bkg6 =(TH1F*)gPad->GetPrimitive("pt1_bkg6");
	
	F6 = pt1_bkg6->Integral();  
	
	tree_bkg7->Draw("(pt1)>>pt1_bkg7(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg7 =(TH1F*)gPad->GetPrimitive("pt1_bkg7");
	F7 = pt1_bkg7->Integral();
	tree_bkg8->Draw("(pt1)>>pt1_bkg8(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                            
	TH1F *pt1_bkg8 =(TH1F*)gPad->GetPrimitive("pt1_bkg8");
	F8 = pt1_bkg8->Integral();
	tree_bkg9->Draw("(pt1)>>pt1_bkg9(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg9 =(TH1F*)gPad->GetPrimitive("pt1_bkg9");
	F9 = pt1_bkg9->Integral();
	if(F9<0) F9=0;
	tree_bkg10->Draw("(pt1)>>pt1_bkg10(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg10 =(TH1F*)gPad->GetPrimitive("pt1_bkg10");
	F10 = pt1_bkg10->Integral();
	tree_bkg11->Draw("(pt1)>>pt1_bkg11(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                       
	TH1F *pt1_bkg11 =(TH1F*)gPad->GetPrimitive("pt1_bkg11");
	F11 = pt1_bkg11->Integral();
	tree_bkg12->Draw("(pt1)>>pt1_bkg12(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg12 =(TH1F*)gPad->GetPrimitive("pt1_bkg12");
	F12 = pt1_bkg12->Integral();
	
	
	
	B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12;
	
	
	Eff = S/Stot_pt1;
	Sig3 = Eff/(1.5+sqrt(B));
	Sig5 = Eff/(2.5+sqrt(B));
      }
      hptggVSmet->Fill(cut_ptgg,cut_met,Sig5);
    }
  }
  
  
  /////////////////////////////////////////////////////////
  Cut_met = Form("t1pfmet>%lf",cut_metIN);
  Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1IN);
  for(i = 0; i<NstepsPtgg;i++){
    cut_ptgg = firstCutPtgg + 5*i;
    Cut_ptgg = Form("ptgg>%lf",cut_ptgg);
    for(j = 0; j<NstepsPt2;j++){
      cut_pt2 = firstCutPt2 + 0.05*j+0.001;
      Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2);
      tree_data->Draw("(pt1)>>pt1_data(30,0,1000)","(mgg<115||mgg>135)"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
      TH1F *pt1_data =(TH1F*)gPad->GetPrimitive("pt1_data");
      N = pt1_data->GetEntries();
      Sig5= 0.;
      if(N>=0){
	tree_sig->Draw("(pt1)>>pt1_1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_1 =(TH1F*)gPad->GetPrimitive("pt1_1"); 
	pt1_1->Scale(crossSec);   	
	S = pt1_1->Integral();   
	tree_bkg1->Draw("(pt1)>>pt1_bkg1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg1 =(TH1F*)gPad->GetPrimitive("pt1_bkg1");
	
	F1 = pt1_bkg1->Integral();   
	tree_bkg2->Draw("(pt1)>>pt1_bkg2(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg2 =(TH1F*)gPad->GetPrimitive("pt1_bkg2");
	
	F2 = pt1_bkg2->Integral();
	if(F2<0) F2 = 0;     
	tree_bkg3->Draw("(pt1)>>pt1_bkg3(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg3 =(TH1F*)gPad->GetPrimitive("pt1_bkg3");
	
	F3 = pt1_bkg3->Integral();  
	tree_bkg4->Draw("(pt1)>>pt1_bkg4(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg4 =(TH1F*)gPad->GetPrimitive("pt1_bkg4");
	
	F4 = pt1_bkg4->Integral();    
	tree_bkg5->Draw("(pt1)>>pt1_bkg5(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg5 =(TH1F*)gPad->GetPrimitive("pt1_bkg5");
	
	F5 = pt1_bkg5->Integral();    
	tree_bkg6->Draw("(pt1)>>pt1_bkg6(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut)); 
	TH1F *pt1_bkg6 =(TH1F*)gPad->GetPrimitive("pt1_bkg6");

	F6 = pt1_bkg6->Integral();  
	
	tree_bkg7->Draw("(pt1)>>pt1_bkg7(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg7 =(TH1F*)gPad->GetPrimitive("pt1_bkg7");
	F7 = pt1_bkg7->Integral();
	tree_bkg8->Draw("(pt1)>>pt1_bkg8(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                            
	TH1F *pt1_bkg8 =(TH1F*)gPad->GetPrimitive("pt1_bkg8");
	F8 = pt1_bkg8->Integral();
	tree_bkg9->Draw("(pt1)>>pt1_bkg9(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg9 =(TH1F*)gPad->GetPrimitive("pt1_bkg9");
	F9 = pt1_bkg9->Integral();
	if(F9<0) F9=0;
	tree_bkg10->Draw("(pt1)>>pt1_bkg10(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg10 =(TH1F*)gPad->GetPrimitive("pt1_bkg10");
	F10 = pt1_bkg10->Integral();
	tree_bkg11->Draw("(pt1)>>pt1_bkg11(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                       
	TH1F *pt1_bkg11 =(TH1F*)gPad->GetPrimitive("pt1_bkg11");
	F11 = pt1_bkg11->Integral();
	tree_bkg12->Draw("(pt1)>>pt1_bkg12(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg12 =(TH1F*)gPad->GetPrimitive("pt1_bkg12");
	F12 = pt1_bkg12->Integral();
	
      
      
	B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12;
	
	
	Eff = S/Stot_pt1;
	Sig3 = Eff/(1.5+sqrt(B));
	Sig5 = Eff/(2.5+sqrt(B));
	}
	hptggVSpt2->Fill(cut_ptgg,cut_pt2,Sig5);
    }
  }
  
  /////////////////////////////////////////////////////////
  Cut_met = Form("t1pfmet>%lf",cut_metIN);
  Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2IN);
  for(i = 0; i<NstepsPtgg;i++){
    cut_ptgg = firstCutPtgg + 5*i;
    Cut_ptgg = Form("ptgg>%lf",cut_ptgg);
    for(j = 0; j<NstepsPt1;j++){
      cut_pt1 = firstCutPt1 + 0.05*j+0.001;
      Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1);
      tree_data->Draw("(pt1)>>pt1_data(30,0,1000)","(mgg<115||mgg>135)"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
      TH1F *pt1_data =(TH1F*)gPad->GetPrimitive("pt1_data");
      N = pt1_data->GetEntries();
      Sig5= 0.;
      if(N>=0){
	tree_sig->Draw("(pt1)>>pt1_1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_1 =(TH1F*)gPad->GetPrimitive("pt1_1"); 
	pt1_1->Scale(crossSec);   	
	S = pt1_1->Integral();   
	tree_bkg1->Draw("(pt1)>>pt1_bkg1(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg1 =(TH1F*)gPad->GetPrimitive("pt1_bkg1");
	
	F1 = pt1_bkg1->Integral();   
	tree_bkg2->Draw("(pt1)>>pt1_bkg2(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg2 =(TH1F*)gPad->GetPrimitive("pt1_bkg2");
	
	F2 = pt1_bkg2->Integral();
	if(F2<0) F2 = 0;     
	tree_bkg3->Draw("(pt1)>>pt1_bkg3(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg3 =(TH1F*)gPad->GetPrimitive("pt1_bkg3");

	F3 = pt1_bkg3->Integral();  
	tree_bkg4->Draw("(pt1)>>pt1_bkg4(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut));
	TH1F *pt1_bkg4 =(TH1F*)gPad->GetPrimitive("pt1_bkg4");
	
	F4 = pt1_bkg4->Integral();    
	tree_bkg5->Draw("(pt1)>>pt1_bkg5(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut && genmatch));
	TH1F *pt1_bkg5 =(TH1F*)gPad->GetPrimitive("pt1_bkg5");
	
	F5 = pt1_bkg5->Integral();    
	tree_bkg6->Draw("(pt1)>>pt1_bkg6(30,0,1000)","weight"*(Cut_pt2 && Cut_met && Cut_pt1 && Cut_ptgg && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut &&nMuonCut &&nEleCut)); 
	TH1F *pt1_bkg6 =(TH1F*)gPad->GetPrimitive("pt1_bkg6");

	F6 = pt1_bkg6->Integral();  
	
	tree_bkg7->Draw("(pt1)>>pt1_bkg7(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg7 =(TH1F*)gPad->GetPrimitive("pt1_bkg7");
	F7 = pt1_bkg7->Integral();
	tree_bkg8->Draw("(pt1)>>pt1_bkg8(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                            
	TH1F *pt1_bkg8 =(TH1F*)gPad->GetPrimitive("pt1_bkg8");
	F8 = pt1_bkg8->Integral();
	tree_bkg9->Draw("(pt1)>>pt1_bkg9(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                             
	TH1F *pt1_bkg9 =(TH1F*)gPad->GetPrimitive("pt1_bkg9");
	F9 = pt1_bkg9->Integral();
	if(F9<0) F9=0;
	tree_bkg10->Draw("(pt1)>>pt1_bkg10(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg10 =(TH1F*)gPad->GetPrimitive("pt1_bkg10");
	F10 = pt1_bkg10->Integral();
	tree_bkg11->Draw("(pt1)>>pt1_bkg11(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut && minDPHIJetcut && nMuonCut && nEleCut ));                                       
	TH1F *pt1_bkg11 =(TH1F*)gPad->GetPrimitive("pt1_bkg11");
	F11 = pt1_bkg11->Integral();
	tree_bkg12->Draw("(pt1)>>pt1_bkg12(30,0,1000)","weight"*(Cut_pt1 && Cut_ptgg && Cut_pt2 && Cut_met && mggmin && mggmax && eveto && maxDPHIJetcut  && DPHIcut &&minDPHIJetcut && nMuonCut && nEleCut ));                                           
	TH1F *pt1_bkg12 =(TH1F*)gPad->GetPrimitive("pt1_bkg12");
	F12 = pt1_bkg12->Integral();
	
      
      
	B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12;
	
	
	Eff = S/Stot_pt1;
	Sig3 = Eff/(1.5+sqrt(B));
	Sig5 = Eff/(2.5+sqrt(B));
	}
	hpt1VSptgg->Fill(cut_pt1,cut_ptgg,Sig5);
      }
  }
  
   
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
  int iPos=0; 
  
  TCanvas *fancycanvas1 = new TCanvas("fancycanvas1","",550,500);
  
  fancycanvas1->SetRightMargin(0.17);   
  hpt1VSpt2->GetXaxis()->SetTitle("p_{T1}/m_{#gamma#gamma}");  
  hpt1VSpt2->GetYaxis()->SetTitle("p_{T2}/m_{#gamma#gamma}");  
  hpt1VSpt2->GetZaxis()->SetTitle("Sign(5#sigma)");  
  
  hpt1VSpt2->GetZaxis()->SetTitleOffset(1.6);  
  hpt1VSpt2->GetYaxis()->SetTitleOffset(1.2);  
  hpt1VSpt2->Draw("colz");
  // CMS_lumi( fancycanvas1,true,iPos,true );
  gStyle->SetOptStat(0); 
  
  
  TCanvas *fancycanvas2 = new TCanvas("fancycanvas2","",550,500);
  
  
  fancycanvas2->SetRightMargin(0.17);   
  hpt1VSmet->GetXaxis()->SetTitle("p_{T1}/m_{#gamma#gamma}");  
  hpt1VSmet->GetYaxis()->SetTitle("MET [GeV]");  
  hpt1VSmet->GetZaxis()->SetTitle("Sign(5#sigma)");  
  hpt1VSmet->GetZaxis()->SetTitleOffset(1.6);  
  hpt1VSmet->GetYaxis()->SetTitleOffset(1.2);  
  hpt1VSmet->Draw("colz");
  // CMS_lumi( fancycanvas2,true,iPos,true );
  gStyle->SetOptStat(0);  
  
  
  TCanvas *fancycanvas3 = new TCanvas("fancycanvas3","",550,500); 
  
  
  fancycanvas3->SetRightMargin(0.17);   
  hmetVSpt2->GetXaxis()->SetTitle("MET [GeV]");  
  hmetVSpt2->GetYaxis()->SetTitle("p_{T2}/m_{#gamma#gamma}");  
  hmetVSpt2->GetZaxis()->SetTitle("Sign(5#sigma)");  
  
  hmetVSpt2->GetZaxis()->SetTitleOffset(1.6);  
  hmetVSpt2->GetYaxis()->SetTitleOffset(1.2);  
  hmetVSpt2->Draw("colz");
  // CMS_lumi( fancycanvas3,true,iPos,true );
  
  gStyle->SetOptStat(0);  
  
  
  TCanvas *fancycanvas4 = new TCanvas("fancycanvas4","",550,500); 
  
  
  fancycanvas4->SetRightMargin(0.17);   
  hptggVSmet->GetXaxis()->SetTitle("p_{T#gamma#gamma} [GeV]");  
  hptggVSmet->GetYaxis()->SetTitle("MET [GeV]");  
  hptggVSmet->GetZaxis()->SetTitle("Sign(5#sigma)");  
  
  hptggVSmet->GetZaxis()->SetTitleOffset(1.6);  
  hptggVSmet->GetYaxis()->SetTitleOffset(1.2);  
  hptggVSmet->Draw("colz");
  // CMS_lumi( fancycanvas4,true,iPos,true );
  
  gStyle->SetOptStat(0);  
  
  
  TCanvas *fancycanvas5 = new TCanvas("fancycanvas5","",550,500); 
  
  
  fancycanvas5->SetRightMargin(0.17);   
  hptggVSpt2->GetXaxis()->SetTitle("p_{T#gamma#gamma} [GeV]");  
  hptggVSpt2->GetYaxis()->SetTitle("p_{T2}/m_{#gamma#gamma}");  
  hptggVSpt2->GetZaxis()->SetTitle("Sign(5#sigma)");  
  
  hptggVSpt2->GetZaxis()->SetTitleOffset(1.6);  
  hptggVSpt2->GetYaxis()->SetTitleOffset(1.2);  
  hptggVSpt2->Draw("colz");
  // CMS_lumi( fancycanvas5,true,iPos,true );
  
  gStyle->SetOptStat(0);  
  

  TCanvas *fancycanvas6 = new TCanvas("fancycanvas6","",550,500); 
  
  
  fancycanvas6->SetRightMargin(0.17);   
  hpt1VSptgg->GetXaxis()->SetTitle("p_{T1}/m_{#gamma#gamma}");  
  hpt1VSptgg->GetYaxis()->SetTitle("p_{T#gamma#gamma} [GeV]");  
  hpt1VSptgg->GetZaxis()->SetTitle("Sign(5#sigma)");  

  hpt1VSptgg->GetZaxis()->SetTitleOffset(1.6);  
  hpt1VSptgg->GetYaxis()->SetTitleOffset(1.2);  
  hpt1VSptgg->Draw("colz");
  //  CMS_lumi( fancycanvas6,true,iPos,true );

  gStyle->SetOptStat(0);  
  
  
  /* fancycanvas1->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPt1vsPt2_mass_DPHIandLepVeto_all4D.png",sigMass));
    fancycanvas2->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPt1vsMET_mass_DPHIandLepVeto_all4D.png",sigMass));
    fancycanvas3->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsMETvsPt2_mass_DPHIandLepVeto_all4D.png",sigMass));
    fancycanvas4->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPtggvsMET_mass_DPHIandLepVeto_all4D.png",sigMass));
    fancycanvas5->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPtggvsPt2_mass_DPHIandLepVeto_all4D.png",sigMass));
    fancycanvas6->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPt1vsPtgg_mass_DPHIandLepVeto_all4D.png",sigMass));*/
  fancycanvas1->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPt1vsPt2_mass_DPHIandLepVeto_all4D_nolabel.pdf",sigMass));
  fancycanvas2->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPt1vsMET_mass_DPHIandLepVeto_all4D_nolabel.pdf",sigMass));
  // fancycanvas3->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsMETvspt2_mass_DPHIandLepVeto_all4D.pdf",sigMass));
  fancycanvas4->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPtggvsMET_mass_DPHIandLepVeto_all4D_nolabel.pdf",sigMass));
  // fancycanvas5->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPtggvsPt2_mass_DPHIandLepVeto_all4D.pdf",sigMass));
  // fancycanvas6->SaveAs(Form("./25ns_2246inv_v3/plots/cuts_MassDep/M%d_CutsPt1vsPtgg_mass_DPHIandLepVeto_all4D.pdf",sigMass));
  
  

}

 
  


