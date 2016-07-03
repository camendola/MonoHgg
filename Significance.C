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



double errSig5(double errB, double errS, double errN, double S, double B, double N){
  double A = pow(errS,2)*pow(N*(2.5+sqrt(B)),-2);
  double BI = pow(errN,2)*pow(S,2)*pow((pow(N,2)*(2.5+sqrt(B))),-2);
  double C = pow(errB,2)*pow(((S/N)*0.5*pow((sqrt(B)+2.5),-2)),2)/B;
  double err = sqrt(A+BI+C); 
  return err;
}

double errSig(double errB, double errS, double S, double B){
  double A = pow(errS,2)*(pow((B+0.5*S),2))/pow((B+S),3);
  double BI = pow(errB,2)*pow((0.5*S),2)*pow((S+B),-3);
  double err = sqrt(A+BI); 
  return err;
}




void Significance(int sigMass, int A0mass,  double cut_pt1IN, double cut_pt2IN, double cut_metIN, double cut_ptggIN){



  TFile *sig = TFile::Open(Form("../25ns_76/25ns_2246inv_A0/2HDM_mZP%d_mA0%d_sel.root",sigMass,A0mass),"READ");
  TFile *bkg1 =  TFile::Open("../25ns_76/25ns_2246inv_76/DiPhoton_sel.root","READ"); 
  TFile *bkg2 =  TFile::Open("../25ns_76/25ns_2246inv_76/DYJetsToLL_sel.root","READ");
  TFile *bkg3 =  TFile::Open("../25ns_76/25ns_2246inv_76/GJets_sel.root","READ");  
  TFile *bkg4 =  TFile::Open("../25ns_76/25ns_2246inv_76/GluGluHToGG_sel.root","READ");
  TFile *bkg5 =  TFile::Open("../25ns_76/25ns_2246inv_76/QCD_sel.root","READ");
  TFile *bkg6 =  TFile::Open("../25ns_76/25ns_2246inv_76/VH_sel.root","READ");
  TFile *bkg7 =  TFile::Open("../25ns_76/25ns_2246inv_76/ttHJetToGG_sel.root","READ");
  TFile *bkg8 =  TFile::Open("../25ns_76/25ns_2246inv_76/VBFHToGG_sel.root","READ");
  TFile *bkg9 =  TFile::Open("../25ns_76/25ns_2246inv_76/TGJets_sel.root","READ");
  TFile *bkg10 =  TFile::Open("../25ns_76/25ns_2246inv_76/TTGJets_sel.root","READ");
  TFile *bkg11 =  TFile::Open("../25ns_76/25ns_2246inv_76/WGToLNuG_sel.root","READ");
  TFile *bkg12 =  TFile::Open("../25ns_76/25ns_2246inv_76/ZGTo2LG_sel.root","READ");
  TFile *bkg13 =  TFile::Open("../25ns_76/25ns_2246inv_76/ZZTo2L2Nu_sel.root","READ");
  TFile *bkg14 =  TFile::Open("../25ns_76/25ns_2246inv_76/TTGG_0Jets_sel.root","READ");


 

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
 TTree *tree_bkg13 = (TTree*) bkg13->Get("DiPhotonTree");
  TTree *tree_bkg14 = (TTree*) bkg14->Get("DiPhotonTree"); 

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

  
  TCanvas *canvas = new TCanvas("canvas","",600,400);   
  int i = 0;
  int j = 0;
  int k = 0;

  double cut_pt1 = 0;
  double cut_pt2 = 0;
  double cut_met = 0;

   
  
  double Stot= 0;

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
 double F13 = 0;
 double F14 = 0;
  double B = 0;
  double Eff = 0;
  
  double Sig3 = 0;
  double Sig5 = 0;
  


 
  


     
  
  
  cout<<"Signal mZP"<<sigMass<<" mA0"<<A0mass<<endl;
  cout<<"120 < mgg < 130"<<endl;  
  cout<<"*****"<<end;

  ios_base::fmtflags oldflags = std::cout.flags();
  
  
  tree_sig->Draw("(mgg)>>mgg_tot1(30,100,180)","weight*((mgg>120)&&(mgg<130))"*(eveto));
  TH1F *mgg_tot1 =(TH1F*)gPad->GetPrimitive("mgg_tot1");

  Stot = mgg_tot1->Integral();
  

  Cut_met = Form("corrMet>%lf",cut_metIN);
  Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2IN);
  Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1IN);
  Cut_ptgg = Form("ptgg>%lf",cut_ptggIN);

  
  tree_sig->Draw("(mgg)>>mgg_1(30,100,180)","weight*(commonSel==1)*((mgg<130)&&(mgg>120))"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && eveto));
  TH1F *mgg_1 =(TH1F*)gPad->GetPrimitive("mgg_1");
  tree_bkg1->Draw("(mgg)>>mgg_bkg1(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg1 =(TH1F*)gPad->GetPrimitive("mgg_bkg1");
  tree_bkg2->Draw("(mgg)>>mgg_bkg2(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto ));
  TH1F *mgg_bkg2 =(TH1F*)gPad->GetPrimitive("mgg_bkg2");
  tree_bkg3->Draw("(mgg)>>mgg_bkg3(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto  && genmatch));
  TH1F *mgg_bkg3 =(TH1F*)gPad->GetPrimitive("mgg_bkg3");
  tree_bkg4->Draw("(mgg)>>mgg_bkg4(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg4 =(TH1F*)gPad->GetPrimitive("mgg_bkg4");
  tree_bkg5->Draw("(mgg)>>mgg_bkg5(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto && genmatch));
  TH1F *mgg_bkg5 =(TH1F*)gPad->GetPrimitive("mgg_bkg5");
  tree_bkg6->Draw("(mgg)>>mgg_bkg6(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg6 =(TH1F*)gPad->GetPrimitive("mgg_bkg6");
  tree_bkg7->Draw("(mgg)>>mgg_bkg7(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg7 =(TH1F*)gPad->GetPrimitive("mgg_bkg7");
  tree_bkg8->Draw("(mgg)>>mgg_bkg8(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg8 =(TH1F*)gPad->GetPrimitive("mgg_bkg8");
  tree_bkg9->Draw("(mgg)>>mgg_bkg9(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg9 =(TH1F*)gPad->GetPrimitive("mgg_bkg9");
  tree_bkg10->Draw("(mgg)>>mgg_bkg10(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg10 =(TH1F*)gPad->GetPrimitive("mgg_bkg10");
  tree_bkg11->Draw("(mgg)>>mgg_bkg11(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto));
  TH1F *mgg_bkg11 =(TH1F*)gPad->GetPrimitive("mgg_bkg11");
  tree_bkg12->Draw("(mgg)>>mgg_bkg12(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto ));
  TH1F *mgg_bkg12 =(TH1F*)gPad->GetPrimitive("mgg_bkg12");
  tree_bkg13->Draw("(mgg)>>mgg_bkg13(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto ));
  TH1F *mgg_bkg13 =(TH1F*)gPad->GetPrimitive("mgg_bkg13");
  tree_bkg14->Draw("(mgg)>>mgg_bkg14(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_ptgg && Cut_pt2 && Cut_pt1 && mggmin && mggmax && eveto ));
  TH1F *mgg_bkg14 =(TH1F*)gPad->GetPrimitive("mgg_bkg14");

  S = mgg_1->Integral();
  F1 = mgg_bkg1->Integral();
  F2 = mgg_bkg2->Integral();
  F3 = mgg_bkg3->Integral();
  F4 = mgg_bkg4->Integral();
  F5 = mgg_bkg5->Integral();
  F6 = mgg_bkg6->Integral();
  F7 = mgg_bkg7->Integral();
  F8 = mgg_bkg8->Integral();
  F9 = mgg_bkg9->Integral();
  F10 = mgg_bkg10->Integral();
  F11 = mgg_bkg11->Integral();
  F12 = mgg_bkg12->Integral();
  F13 = mgg_bkg13->Integral();
  F14 = mgg_bkg14->Integral();
  
  if(F1 < 0) F1 = 0;  
  if(F2 < 0) F2 = 0;  
  if(F3 < 0) F3 = 0;  
  if(F4 < 0) F4 = 0;  
  if(F5 < 0) F5 = 0;  
  if(F6 < 0) F6 = 0;  
  if(F7 < 0) F7 = 0;  
  if(F8 < 0) F8 = 0;  
  if(F9 < 0) F9 = 0;  
  if(F10 < 0) F10 = 0;  
  if(F11 < 0) F11 = 0;  
  if(F12 < 0) F12 = 0;  
  if(F13 < 0) F13 = 0;  
  if(F14 < 0) F14 = 0;  
  
  
  double DeltaS = 0;
  double DeltaStot = 0;
  double DeltaF1 = 0;
  double DeltaF2 = 0;
  double DeltaF3 = 0;
  double DeltaF4 = 0;
  double DeltaF5 = 0;
  double DeltaF6 = 0;
  double DeltaF7 = 0;
  double DeltaF8 = 0;
  double DeltaF9 = 0;
  double DeltaF10 = 0;
  double DeltaF11 = 0;
  double DeltaF12 = 0;
  double DeltaF13 = 0;
  double DeltaF14 = 0;
  
  double DeltaSi= 0;
  double DeltaStoti= 0;
  double DeltaF1i= 0;
  double DeltaF2i= 0;
  double DeltaF3i = 0;
  double DeltaF4i = 0;
  double DeltaF5i = 0;
  double DeltaF6i = 0;
  double DeltaF7i = 0;
  double DeltaF8i = 0;
  double DeltaF9i = 0;
  double DeltaF10i = 0;
  double DeltaF11i = 0;
  double DeltaF12i = 0;
double DeltaF13i = 0;
  double DeltaF14i = 0;
  
  for(i=0; i<31; i++){
    DeltaSi=mgg_1->GetBinError(i);      
    DeltaStoti=mgg_tot1->GetBinError(i);      
    DeltaF1i=mgg_bkg1->GetBinError(i);      
    DeltaF2i= mgg_bkg2->GetBinError(i);
    DeltaF3i= mgg_bkg3->GetBinError(i);
    DeltaF4i= mgg_bkg4->GetBinError(i);      
    DeltaF5i= mgg_bkg5->GetBinError(i);      
    DeltaF6i= mgg_bkg6->GetBinError(i);      
    DeltaF7i= mgg_bkg7->GetBinError(i);      
    DeltaF8i= mgg_bkg8->GetBinError(i);      
    DeltaF9i= mgg_bkg9->GetBinError(i);      
    DeltaF10i= mgg_bkg10->GetBinError(i);      
    DeltaF11i= mgg_bkg11->GetBinError(i);      
    DeltaF12i= mgg_bkg12->GetBinError(i);      
    DeltaF13i= mgg_bkg13->GetBinError(i);      
    DeltaF14i= mgg_bkg14->GetBinError(i);      
    
    DeltaS +=  pow(DeltaSi,2);   
    DeltaF1 +=  pow(DeltaF1i,2);
    DeltaF2 +=  pow(DeltaF2i,2);
    DeltaF3 +=  pow(DeltaF3i,2);
    DeltaF4 +=  pow(DeltaF4i,2);
    DeltaF5 +=  pow(DeltaF5i,2);
    DeltaF6 +=  pow(DeltaF6i,2);
    DeltaF7 +=  pow(DeltaF7i,2);
    DeltaF8 +=  pow(DeltaF8i,2);
    DeltaF9 +=  pow(DeltaF9i,2);
    DeltaF10 +=  pow(DeltaF10i,2);
    DeltaF11 +=  pow(DeltaF11i,2);
    DeltaF12 +=  pow(DeltaF12i,2);
    DeltaF13 +=  pow(DeltaF13i,2);
    DeltaF14 +=  pow(DeltaF14i,2);
    
  }
  DeltaS  = sqrt(DeltaS);
  DeltaStot  = sqrt(DeltaStot);
  DeltaF1  = sqrt(DeltaF1);
  DeltaF2  = sqrt(DeltaF2);
  DeltaF3  = sqrt(DeltaF3);
  DeltaF4  = sqrt(DeltaF4);
  DeltaF5  = sqrt(DeltaF5);
  DeltaF6  = sqrt(DeltaF6);
  DeltaF7  = sqrt(DeltaF7);
  DeltaF8  = sqrt(DeltaF8);
  DeltaF9  = sqrt(DeltaF9);
  DeltaF10  = sqrt(DeltaF10);
  DeltaF11  = sqrt(DeltaF11);
  DeltaF12  = sqrt(DeltaF12);
  DeltaF13  = sqrt(DeltaF13);
  DeltaF14  = sqrt(DeltaF14);
  
  
  B = (F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12+F13+F14);
  double DeltaB = sqrt(pow(DeltaF1,2)+pow(DeltaF2,2)+pow(DeltaF3,2)+pow(DeltaF4,2)+pow(DeltaF5,2)+pow(DeltaF6,2)+pow(DeltaF7,2)+pow(DeltaF8,2)+pow(DeltaF9,2)+pow(DeltaF10,2)+pow(DeltaF11,2)+pow(DeltaF12,2)+pow(DeltaF13,2)+pow(DeltaF14,2));      
  DeltaB = sqrt(B); 
  
  
  
  Eff = S/Stot;
  Sig3 = Eff/(1.5+sqrt(B)); //punzi's significance 3sigma
  Sig5 = Eff/(2.5+sqrt(B)); //punzi's significance 5sigma
  double DeltaSig5 = errSig5(DeltaB,DeltaS,DeltaStot,S,B,Stot);     
  double Sign = S/sqrt(S+B); 
  double DeltaSig = errSig(DeltaB,DeltaS,S,B);
  cout<<"S = "<<S<<" +/- "<<DeltaS<<endl;
  cout<<"B = "<<B<<" +/- "<<DeltaB<<endl;
  cout<<"Sig5 = "<<Sig5<<" +/- "<<DeltaSig5<<endl;
  cout<<"Eff = "<<Eff<<endl;
  cout<<"Sign = "<<Sign<<" +/- "<<DeltaSig<<endl;
  cout<<"(sigmaB/B) = "<<1./sqrt(B)<<endl;  
     
}




