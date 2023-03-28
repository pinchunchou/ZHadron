#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <cmath>
#include <TCanvas.h>
#include <TLegend.h>
#include <TChain.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;
void style(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetTitleFont(43);
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetStatFont(43);

  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetOptStat(0); /*don't show statistics box*/
  gStyle->SetOptTitle(0); /*don't show histogram titles*/
  gStyle->SetTitleSize(48, "xyz");
  gStyle->SetTitleOffset(1, "xyz");
  gStyle->SetLabelSize(36, "xyz");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetLineScalePS(1.5);

  gROOT->ForceStyle();
}

void AnaSumPfE(){

  style();

  std::cout<<"Start running..."<<std::endl;

	TChain *tMC = new TChain("particleFlowAnalyser/pftree");
  tMC->Add("/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbSpring21MiniAOD-mva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/20221017_ZHadronInitialRunMCDY/221017_113040/0000/*.root");
  
  TChain *tHJ = new TChain("particleFlowAnalyser/pftree"); 
  tHJ->Add("/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/*.root");

  int entries_MC = tMC->GetEntries();
  int entries_HJ = tHJ->GetEntries();

  std::cout<<"entries_MC = "<<entries_MC<<", entries_HJ = "<<entries_HJ<<std::endl;

  vector<float>* pfE_MC = 0;
  vector<float>* pfEta_MC = 0;
  int nPF_MC;
  tMC->SetBranchAddress("pfE", &pfE_MC);
  tMC->SetBranchAddress("nPF", &nPF_MC);
  tMC->SetBranchAddress("pfEta",&pfEta_MC);

  std::cout<<"a"<<std::endl;

  vector<float>* pfE_HJ = 0;
  vector<float>* pfEta_HJ = 0;
  int nPF_HJ;
  tHJ->SetBranchAddress("pfE", &pfE_HJ);
  tHJ->SetBranchAddress("nPF",&nPF_HJ);
  tHJ->SetBranchAddress("pfEta",&pfEta_HJ);

  std::cout<<"b"<<std::endl;

  double *dataMC; double *dataHJ;
  dataMC = (double*) malloc(entries_MC*sizeof(double));
  dataHJ = (double*) malloc(entries_HJ*sizeof(double));
  //double dataMC[entries_MC];
  //std::cout<<"c"<<std::endl;
  //double dataHJ[entries_HJ];

  std::cout<<"Filling MC"<<std::endl;

  double sumPfE_tmp = 0;
  int count_no0_MC = 0, count_no0_HJ = 0;

  for (int i = 0; i < entries_MC; i++){
    //std::cout<<"i = "<<i<<std::endl;
    if (i%1000==0) cout <<i<<"/"<<entries_MC<<endl;
    tMC->GetEntry(i);
    //std::cout<<"d"<<std::endl;
    sumPfE_tmp = 0;
    //std::cout<<"nPF_MC = "<<nPF_MC<<std::endl;
    for(int ipf=0;ipf<nPF_MC;ipf++){
      //std::cout<<"ipf = "<<ipf<<std::endl;
      if(fabs((*pfEta_MC)[ipf])>3.&&fabs((*pfEta_MC)[ipf])<5.) sumPfE_tmp+=(*pfE_MC)[ipf];
    }
    //std::cout<<"Storing..."<<std::endl;
    dataMC[i] = sumPfE_tmp;
    if(sumPfE_tmp>0) count_no0_MC++;
    //std::cout<<"Stored."<<std::endl;
  }

  std::cout<<"Filling HJ"<<std::endl;

  for (int i = 0; i < entries_HJ; i++){
    if (i%1000==0) cout <<i<<"/"<<entries_HJ<<endl;
    tHJ->GetEntry(i);
    sumPfE_tmp = 0;
    for(int ipf=0;ipf<nPF_HJ;ipf++){
      if(fabs((*pfEta_HJ)[ipf])>3.&&fabs((*pfEta_HJ)[ipf])<5.) sumPfE_tmp+=(*pfE_HJ)[ipf];
    }
    dataHJ[i] = sumPfE_tmp;
    if(sumPfE_tmp>0) count_no0_HJ++;
  }

  std::cout<<"count_no0_MC = "<<count_no0_MC<<", count_no0_HJ = "<<count_no0_HJ<<std::endl;

  std::sort(dataMC, dataMC + entries_MC, std::greater<int>());
  std::sort(dataHJ, dataHJ + entries_HJ, std::greater<int>());

  //0-10,10-20,20-30,...,90-100.

  TH1D *hMC_PfE_avg = new TH1D("hMC_PfE_avg","",10,0,100);
  TH1D *hHJ_PfE_avg = new TH1D("hHJ_PfE_avg","",10,0,100);
  TH1D *hMH_PfE_avg = new TH1D("hMH_PfE_avg","",10,0,100);

  double sum_sum_sum_MC = 0;
  double sum_sum_sum_HJ = 0;

  std::cout<<"Filling hist!"<<std::endl;

  for(int i=0;i<10;i++){

  	std::cout<<"i = "<<i<<std::endl;

  	double sum_sum_MC = 0;
  	int count_MC = 0;

  	std::cout<<"Filling MC"<<std::endl;

  	for(int j=(int)(i*count_no0_MC/10);j<(int)((i+1)*count_no0_MC/10);j++){
  		sum_sum_MC+=dataMC[j];
  		sum_sum_sum_MC+=dataMC[j];
  		count_MC++;
  	}
  	sum_sum_MC/=count_MC;
    std::cout<<"sum_sum_MC = "<<sum_sum_MC<<std::endl;
  	hMC_PfE_avg->Fill(i*10+5,sum_sum_MC);

  	double sum_sum_HJ = 0;
  	int count_HJ = 0;

  	std::cout<<"Filling HJ"<<std::endl;

  	for(int j=(int)(i*count_no0_HJ/10);j<(int)((i+1)*count_no0_HJ/10);j++){
  		sum_sum_HJ+=dataHJ[j];
  		sum_sum_sum_HJ+=dataHJ[j];
  		count_HJ++;
  	}
  	sum_sum_HJ/=count_HJ;
    std::cout<<"sum_sum_HJ = "<<sum_sum_HJ<<std::endl;
  	hHJ_PfE_avg->Fill(i*10+5,sum_sum_HJ);
    hMH_PfE_avg->Fill(i*10+5,sum_sum_HJ-sum_sum_MC);
  }

  std::cout<<"Drawing..."<<std::endl;

  sum_sum_sum_MC/=count_no0_MC;
  sum_sum_sum_HJ/=count_no0_HJ;

  std::cout<<"<sum_sum_sum_MC> = "<<sum_sum_sum_MC<<std::endl;
  std::cout<<"<sum_sum_sum_HJ> = "<<sum_sum_sum_HJ<<std::endl;

  TCanvas *c1 = new TCanvas("c1","",1200,800);

  hMC_PfE_avg->SetMarkerStyle(20);
  hHJ_PfE_avg->SetMarkerStyle(21);

  hMC_PfE_avg->SetMarkerColor(TColor::GetColor("#377eb8"));
  hHJ_PfE_avg->SetMarkerColor(TColor::GetColor("#e41a1c"));

  hMC_PfE_avg->SetMinimum(0);
  hHJ_PfE_avg->SetMinimum(0);

  hMC_PfE_avg->SetXTitle("Percentage");
  hHJ_PfE_avg->SetXTitle("Percentage");

  hMC_PfE_avg->SetYTitle("Average of HF energy");
  hHJ_PfE_avg->SetYTitle("Average of HF energy");

  hMC_PfE_avg->Draw("hist p");
  hHJ_PfE_avg->Draw("hist p same");

  TLegend leg(0.58,0.78,0.9,0.9);
  leg.AddEntry(hMC_PfE_avg ,"Z+hydjet MC","p");
  leg.AddEntry(hHJ_PfE_avg ,"Hydjet only MC","p");
  leg.SetFillColorAlpha(kWhite,0);
  leg.SetLineColor(kBlack);
  leg.SetLineWidth(1);
  leg.Draw();

  c1->SaveAs("figs/SumPfE.png"); 

  hMC_PfE_avg->SetMinimum(1);
  hHJ_PfE_avg->SetMinimum(1);
  c1->SetLogy(1);

  c1->SaveAs("figs/SumPfE_log.png"); 

  c1->Clear();
  c1->SetLogy(0);

  hMH_PfE_avg->SetMarkerStyle(20);
  hMH_PfE_avg->SetMarkerColor(kBlack);
  hMH_PfE_avg->SetMinimum(0);

  hMH_PfE_avg->SetXTitle("Percentage");
  hMH_PfE_avg->SetYTitle("Difference of HF energy average");

  hMH_PfE_avg->Draw("hist p");

  c1->SaveAs("figs/SumPfE_diff.png"); 

  free(dataMC);
  free(dataHJ);
}