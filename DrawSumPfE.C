#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCut.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <cmath>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TChain.h>
#include <TLine.h>
#include <iostream>
#include <algorithm>
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

void DrawSumPfE(){
	TChain *tMC = new TChain("t");
   	tMC->Add("/eos/cms/store/group/phys_heavyions_ops/pchou/MC/*.root?#t");
   
   	TChain *tHJ = new TChain("t");
   	tHJ->Add("/eos/cms/store/group/phys_heavyions_ops/pchou/Hydjet/output_MC0000.root?#t");

   	int entries_MC = tMC->GetEntries();
   	int entries_HJ = tHJ->GetEntries();

   	double dataMC[entries_MC];
   	double dataHJ[entries_HJ];

   	double sumPfE_MC;
    tMC->SetBranchAddress("sumPfE", &sumPfE_MC);

    double sumPfE_HJ;
    tHJ->SetBranchAddress("sumPfE", &sumPfE_HJ);


   	for(int i=0;i<entries_MC;i++){
   		tMC->GetEntry(i);
   		dataMC[i] = sumPfE_MC;
   	}
   	for(int i=0;i<entries_HJ;i++){
   		tHJ->GetEntry(i);
   		dataHJ[i] = sumPfE_HJ;
   	}

   	std::sort(dataMC, dataMC + entries_MC, std::greater<int>());
   	std::sort(dataHJ, dataHJ + entries_HJ, std::greater<int>());

   	//0-10,10-20,20-30,...,90-100.

   	TH1D *hMC_PfE_avg = new TH1D("hMC_PfE_avg","",10,0,100);
   	TH1D *hHJ_PfE_avg = new TH1D("hHJ_PfE_avg","",10,0,100);

   	double sum_sum_sum_MC = 0;
   	double sum_sum_sum_HJ = 0;

   	for(int i=0;i<10;i++){

   		double sum_sum_MC = 0;
   		int count_MC = 0;
   		for(int j=(int)(i*entries_MC/10);j<(int)((i+1)*entries_MC/10);j++){
   			sum_sum_MC+=dataMC[j];
   			sum_sum_sum_MC+=dataMC[j];
   			count_MC++;
   		}
   		sum_sum_MC/=count_MC;
   		hMC_PfE_avg->Fill(i*10+5,sum_sum_MC);

   		double sum_sum_HJ = 0;
   		int count_HJ = 0;

   		for(int j=(int)(i*entries_HJ/10);j<(int)((i+1)*entries_HJ/10);j++){
   			sum_sum_HJ+=dataHJ[j];
   			sum_sum_sum_HJ+=dataHJ[j];
   			count_HJ++;
   		}
   		sum_sum_HJ/=count_HJ;
   		hHJ_PfE_avg->Fill(i*10+5,sum_sum_HJ);
   	}

   	sum_sum_sum_MC/=entries_MC;
   	sum_sum_sum_HJ/=entries_HJ;

   	std::cout<<"sum_sum_sum_MC = "<<sum_sum_sum_MC<<std::endl;
   	std::cout<<"sum_sum_sum_HJ = "<<sum_sum_sum_HJ<<std::endl;

   	TCanvas *c1 = new TCanvas("c1","",800,800);

   	hMC_PfE_avg->SetMarkerStyle(24);
   	hHJ_PfE_avg->SetMarkerStyle(25);

   	hMC_PfE_avg->SetMarkerColor(TColor::GetColor("#377eb8"));
   	hHJ_PfE_avg->SetMarkerColor(TColor::GetColor("#e41a1c"));

   	hMC_PfE_avg->SetMinimum(0);
   	hHJ_PfE_avg->SetMinimum(0);

   	hMC_PfE_avg->SetXTitle("Percentage");
   	hHJ_PfE_avg->SetXTitle("Percentage");

   	hMC_PfE_avg->SetYTitle("Average of HF energy");
   	hHJ_PfE_avg->SetYTitle("Average of HF energy");

   	hMC_PfE_avg->Draw("p");
   	hHJ_PfE_avg->Draw("p same");

   	TLegend leg(0.18,0.68,0.58,0.9);
   	leg.AddEntry(hMC_PfE_avg ,"Z+hydjet MC","p");
   	leg.AddEntry(hHJ_PfE_avg ,"Hydjet only MC","p");
   	leg.SetFillColorAlpha(kWhite,0);
   	leg.SetLineColor(kBlack);
   	leg.SetLineWidth(1);
   	leg.Draw();


   	c1->SaveAs("figs/SumPfE.png"); 




}