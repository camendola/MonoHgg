#include <TH1.h>
#include <THStack.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>

#include <iostream>
#include <fstream>
#include <iomanip>



void OptScan(int sigMass,int A0mass){

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
  
  TCut eveto1 = "eleveto1==1";
  TCut eveto2 = "eleveto2==1";
  TCut eveto = eveto1 && eveto2;
  TCut genmatch = "((genmatch1==1 && genmatch2==0)||(genmatch1==0 && genmatch2==1)||(genmatch1==0 && genmatch2==0))";  
 
  
  TCanvas *canvas = new TCanvas("canvas","",600,400);   
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  
  double  cut_pt1 = 0;
  double  cut_pt2 = 0;
  double  cut_met = 0;
  double  cut_ptgg = 0;  
  
  double Stot = 0;
 
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
  int NstepsPt1 = 10; 
  int NstepsPt2 = 10; 
  int NstepsMET = 10; 
  int NstepsPtgg =10; 
  double signMax = 0;
  double OldSig5 = 0;
  double optCut_pt1 = 0;
  double optCut_pt2 = 0;
  double optCut_met = 0;
  double optCut_ptgg = 0;
  double temp = 0;
  
  double firstCut[4];
  double crossSec = 1.;
  
  if(sigMass==600){
    firstCut[0]=0.4; 
    firstCut[1]=0.25; 
    firstCut[2]=50;
    firstCut[3]=50;
    
  } 
  
  if(sigMass==800){
    firstCut[0]=0.4; 
    firstCut[1]=0.25; 
    firstCut[2]=50;
    firstCut[3]=50;
  } 
  
  
  if(sigMass==1000){
    firstCut[0]=0.4; 
    firstCut[1]=0.25; 
    firstCut[2]=50;
    firstCut[3]=50;
  }
  
  if(sigMass==1200){
    firstCut[0]=0.4; 
    firstCut[1]=0.25; 
    firstCut[2]=50;
  firstCut[3]=50;

  }
  
  if(sigMass==1400){
    firstCut[0]=0.4; 
    firstCut[1]=0.25; 
    firstCut[2]=50;
 firstCut[3]=50;

  }


  if(sigMass==1700){
    firstCut[0]=1.3; 
    firstCut[1]=0.25; 
    firstCut[2]=280;
 firstCut[3]=300;

  }



  if(sigMass==2500){
    firstCut[0]=1.3; 
    firstCut[1]=0.25; 
    firstCut[2]=240;
 firstCut[3]=240;

  }
  
  
  cout<<"############Segnale M"<<sigMass<<endl;
  cout<<"120 < mgg < 130"<<endl;  
  ios_base::fmtflags oldflags = std::cout.flags();

  
  
  
  
  for(i = 0; i<NstepsPt1;i++){
    
    
    cut_pt1 = firstCut[0] + 0.05*i;
    Cut_pt1 = Form("pt1/mgg>%lf",cut_pt1);
    
    for(j = 0; j<NstepsPt2;j++){
      
      cut_pt2 = firstCut[1] + 0.05*j;
      Cut_pt2 = Form("pt2/mgg>%lf",cut_pt2);
      
      for(k = 0; k<NstepsMET;k++){
	cut_met = firstCut[2] + 5*k;
	Cut_met = Form("corrMet>%lf",cut_met);

	
	

  tree_sig->Draw("(mgg)>>mgg_tot(30,100,180)","weight"*(mggmin && mggmax && eveto));
      TH1F *mgg_tot =(TH1F*)gPad->GetPrimitive("mgg_tot"); 
      
      Stot = mgg_tot->Integral();  
      cout<<"pt1 > "<<cut_pt1<<"; pt2 > "<<cut_pt2<<"; met > "<<cut_met<<"; Stot ="<<Stot<<endl;
      cout<<"cut\tS\tDiPho\tDY\tGJets\tggH\tQCD\tVH\tB\tEff\tSig(3sigma)\tSig(5sigma)"<<endl;
      for(l = 0; l<NstepsPtgg;l++){
	cut_ptgg = firstCut[3] + 5*l;
	Cut_ptgg = Form("ptgg>%lf",cut_ptgg);

	Sig5 = 0.;

   
	tree_sig->Draw("(mgg)>>mgg_1(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg && mggmin && mggmax));
	TH1F *mgg_1 =(TH1F*)gPad->GetPrimitive("mgg_1"); 
	tree_bkg1->Draw("(mgg)>>mgg_bkg1(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
   TH1F *mgg_bkg1 =(TH1F*)gPad->GetPrimitive("mgg_bkg1");
   tree_bkg2->Draw("(mgg)>>mgg_bkg2(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
   TH1F *mgg_bkg2 =(TH1F*)gPad->GetPrimitive("mgg_bkg2");
   tree_bkg3->Draw("(mgg)>>mgg_bkg3(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax && genmatch));
	TH1F *mgg_bkg3 =(TH1F*)gPad->GetPrimitive("mgg_bkg3");
	tree_bkg4->Draw("(mgg)>>mgg_bkg4(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg4 =(TH1F*)gPad->GetPrimitive("mgg_bkg4");
	tree_bkg5->Draw("(mgg)>>mgg_bkg5(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax && genmatch));
	TH1F *mgg_bkg5 =(TH1F*)gPad->GetPrimitive("mgg_bkg5");
	tree_bkg6->Draw("(mgg)>>mgg_bkg6(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg6 =(TH1F*)gPad->GetPrimitive("mgg_bkg6");
	tree_bkg7->Draw("(mgg)>>mgg_bkg7(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg7 =(TH1F*)gPad->GetPrimitive("mgg_bkg7");
	tree_bkg8->Draw("(mgg)>>mgg_bkg8(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg8 =(TH1F*)gPad->GetPrimitive("mgg_bkg8");	
	tree_bkg9->Draw("(mgg)>>mgg_bkg9(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg9 =(TH1F*)gPad->GetPrimitive("mgg_bkg9");
	tree_bkg10->Draw("(mgg)>>mgg_bkg10(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg10 =(TH1F*)gPad->GetPrimitive("mgg_bkg10");
	tree_bkg11->Draw("(mgg)>>mgg_bkg11(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg11 =(TH1F*)gPad->GetPrimitive("mgg_bkg11");
	tree_bkg12->Draw("(mgg)>>mgg_bkg12(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg12 =(TH1F*)gPad->GetPrimitive("mgg_bkg12");
	tree_bkg13->Draw("(mgg)>>mgg_bkg13(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
	TH1F *mgg_bkg13 =(TH1F*)gPad->GetPrimitive("mgg_bkg13");
	tree_bkg14->Draw("(mgg)>>mgg_bkg14(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
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
	
	B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12+F13+F14;
	Eff = S/Stot;
	Sig3 = Eff/(1.5+sqrt(B));
	Sig5 = Eff/(2.5+sqrt(B));
	cout<<fixed<<setprecision(2);
	cout<<cut_ptgg<<"\t"<<S<<"\t"<<F1<<"\t"<<F2<<"\t"<<F3<<"\t"<<F4<<"\t"<<F5<<"\t"<<F6<<"\t"<<B<<"\t"<<Eff<<"\t";
	cout.flags(oldflags);
	cout.precision(7);
	cout<<"\t"<<Sig3<<"\t"<<Sig5<<endl;









	if(Sig5>OldSig5){
	  temp = Sig5;  
	}      
	if(temp > signMax){
	  signMax = Sig5;
	  optCut_pt1 = cut_pt1;
	  optCut_pt2 = cut_pt2;
	  optCut_met = cut_met; 
	  optCut_ptgg = cut_ptgg;
	}
	OldSig5 = Sig5;
      }    
      cout<<""<<endl;
    }
    cout<<""<<endl;
  }
   cout<<""<<endl;
  }

  
  canvas->cd();
  cout<<"##################"<<endl;
  cout<<"Optimal cuts:"<<endl;
  cout<<"pt1/mgg > "<<optCut_pt1<<endl;
  cout<<"pt2/mgg > "<<optCut_pt2<<endl;
  cout<<"MET > "<<optCut_met<<"GeV"<<endl;
 cout<<"ptgg > "<<optCut_ptgg<<endl;
  Cut_pt1 = Form("pt1/mgg>%lf",optCut_pt1);
  Cut_pt2 = Form("pt2/mgg>%lf",optCut_pt2);
  Cut_met = Form("corrMet>%lf",optCut_met);
Cut_ptgg = Form("ptgg>%lf",optCut_ptgg);
  cout<<"S\tStot\tDiPho\tDY\tGJets\tggH\tQCD\tVH\tB\tEff\tSig(3sigma)\tSig(5sigma)"<<endl;
  tree_sig->Draw("(mgg)>>Mgg_1(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_1 =(TH1F*)gPad->GetPrimitive("Mgg_1"); 
   tree_bkg1->Draw("(mgg)>>Mgg_bkg1(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg1 =(TH1F*)gPad->GetPrimitive("Mgg_bkg1");
  tree_bkg2->Draw("(mgg)>>Mgg_bkg2(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg2 =(TH1F*)gPad->GetPrimitive("Mgg_bkg2");
  tree_bkg3->Draw("(mgg)>>Mgg_bkg3(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax && genmatch));
  TH1F *Mgg_bkg3 =(TH1F*)gPad->GetPrimitive("Mgg_bkg3");
  tree_bkg4->Draw("(mgg)>>Mgg_bkg4(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg4 =(TH1F*)gPad->GetPrimitive("Mgg_bkg4");
  tree_bkg5->Draw("(mgg)>>Mgg_bkg5(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax && genmatch));
  TH1F *Mgg_bkg5 =(TH1F*)gPad->GetPrimitive("Mgg_bkg5");
  tree_bkg6->Draw("(mgg)>>Mgg_bkg6(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg6 =(TH1F*)gPad->GetPrimitive("Mgg_bkg6");
  tree_bkg7->Draw("(mgg)>>Mgg_bkg7(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg7 =(TH1F*)gPad->GetPrimitive("Mgg_bkg7");
  tree_bkg8->Draw("(mgg)>>Mgg_bkg8(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg8 =(TH1F*)gPad->GetPrimitive("Mgg_bkg8");	
  tree_bkg9->Draw("(mgg)>>Mgg_bkg9(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg9 =(TH1F*)gPad->GetPrimitive("Mgg_bkg9");
  tree_bkg10->Draw("(mgg)>>Mgg_bkg10(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg10 =(TH1F*)gPad->GetPrimitive("Mgg_bkg10");
  tree_bkg11->Draw("(mgg)>>Mgg_bkg11(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg11 =(TH1F*)gPad->GetPrimitive("Mgg_bkg11");
  tree_bkg12->Draw("(mgg)>>Mgg_bkg12(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg12 =(TH1F*)gPad->GetPrimitive("Mgg_bkg12");
  tree_bkg13->Draw("(mgg)>>Mgg_bkg13(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg13 =(TH1F*)gPad->GetPrimitive("Mgg_bkg13");
  tree_bkg14->Draw("(mgg)>>Mgg_bkg14(30,100,180)","weight*(commonSel==1)"*(Cut_met && Cut_pt2 && Cut_pt1 && Cut_ptgg  && mggmin && mggmax));
  TH1F *Mgg_bkg14 =(TH1F*)gPad->GetPrimitive("Mgg_bkg14");
  
  S = Mgg_1->Integral();
  F1 = Mgg_bkg1->Integral();  
  F2 = Mgg_bkg2->Integral();  
  F3 = Mgg_bkg3->Integral();  
  F4 = Mgg_bkg4->Integral();  
  F5 = Mgg_bkg5->Integral();  
  F6 = Mgg_bkg6->Integral();  
  F7 = Mgg_bkg7->Integral();  
  F8 = Mgg_bkg8->Integral();  
  F9 = Mgg_bkg9->Integral();  
  F10 = Mgg_bkg10->Integral();  
  F11 = Mgg_bkg11->Integral();  	
  F12 = Mgg_bkg12->Integral();  
  F13 = Mgg_bkg13->Integral();  	
  F14 = Mgg_bkg14->Integral();  
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
  B = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12+F13+F14;
  Eff = S/Stot;
  Sig3 = Eff/(1.5+sqrt(B));
  Sig5 = Eff/(2.5+sqrt(B));
  cout<<fixed<<setprecision(2);
  cout<<S<<"\t"<<Stot<<F1<<"\t"<<F2<<"\t"<<F3<<"\t"<<F4<<"\t"<<F5<<"\t"<<F6<<"\t"<<B<<"\t"<<Eff<<"\t";
  cout.flags(oldflags);
  cout.precision(7);
  cout<<"\t"<<Sig3<<"\t"<<Sig5<<endl;
 
  



 

  
}



 
  
