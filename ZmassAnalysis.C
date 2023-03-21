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

const char *typeofdata = "doubleMu";
const char *typeofdatatext = "Double muon";

//const char *typeofdata = "HardProbes";
//const char *typeofdatatext = "hard probes";

//const char *typeofdata = "MinBias7";
//const char *typeofdatatext = "minimum bias 7";


//const float hf_diff[5] = {24.7924, 484.422, 1342.21, 3039.47, 999999};
//const int cent_diff[5] = {90, 50, 30, 10, 0};

const float hf_diff[5] = {999999, 3039.47, 1342.21, 484.422, 24.7924};
const float cent_diff[5] = {0, 10, 30, 50, 90};

Double_t bwfun(Double_t *x, Double_t *par) {

   //Fit parameters:
   
   //par[0]=mean
   //par[1]=width
   //par[2]=area
   
   //par[3]=width of the convoluted Gaussian function
   
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 1000.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::BreitWigner(xx,par[0], par[1]);
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
//         sum += fland * ROOT::Math::crystalball_function(x[0], par[3], par[4], 1, xx);
         xx = xupp - (i-.5) * step;
         fland = TMath::BreitWigner(xx,par[0], par[1]);
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
//         sum += fland * ROOT::Math::crystalball_function(x[0], par[3], par[4], 1, xx);
      }

      return (par[2] * step * sum * invsq2pi +par[4]*exp(par[6]+par[5]*x[0]));
 //     return (par[2] * step * sum * invsq2pi / par[3]+par[4]+par[5]*x[0]+par[6]*x[0]*x[0]+par[7]*x[0]*x[0]*x[0]);
}

//TH1D* ZmassAnalysis_single(double ptL=0,double ptH=2000,int centL=0,int centH=200)
TH1D* ZmassAnalysis_single(double ptL=0,double ptH=2000,int centL=0,int centH=4)
{

   int binsize = 40;
   if(centL==3&&centH==4) binsize = 30;
   if(ptL>55) binsize = 25;

   TH1D *hData = new TH1D("hData","",binsize,81.2,101.2);
   TH1D *hMC = new TH1D("hMC","",binsize,81.2,101.2);

   TCanvas *c = new TCanvas("c","",800,800);

   TH1D *hData_eta = new TH1D("hData_eta","",binsize,-8,8);
   TH1D *hMC_eta = new TH1D("hMC_eta","",binsize,-8,8);
   
   TH1D *hData_phi = new TH1D("hData_phi","",binsize,-3.1415926,3.1415926);
   TH1D *hMC_phi = new TH1D("hMC_phi","",binsize,-3.1415926,3.1415926);
     
   TChain *tMC = new TChain("t");
   tMC->Add("/eos/cms/store/group/phys_heavyions_ops/pchou/MC/*.root?#t");
   
   TChain *tData = new TChain("t");
   tData->Add(Form("/eos/cms/store/group/phys_heavyions_ops/pchou/%s/*.root?#t",typeofdata));

   TH1D *h = new TH1D("h","",100,0,1);

   //TFile *infData = new TFile("output_doubleMu_221107.root");
   //TFile *infMC = new TFile("outputMC1.root");

   //TFile *infData = new TFile("/eos/cms/store/group/phys_heavyions_ops/pchou/output_doubleMu50.root");
   //TFile *infData = new TFile("/eos/cms/store/group/phys_heavyions_ops/pchou/output_doubleMu0000.root");
   //TFile *infMC = new TFile("/eos/cms/store/group/phys_heavyions_ops/pchou/outputMC.root");

   //TTree *tData = (TTree*)infData->Get("t");
   //TTree *tMC = (TTree*)infMC->Get("t");
   
   TH1D *hData_pt = new TH1D("hData_pt","",binsize,0,200);
   TH1D *hMC_pt = new TH1D("hMC_pt","",binsize,0,200);

   TH1D *hData_muPt1 = new TH1D("hData_muPt1","",40,0,120);
   TH1D *hMC_muPt1 = new TH1D("hMC_muPt1","",40,0,120);

   TH1D *hData_muPt2 = new TH1D("hData_muPt2","",40,0,120);
   TH1D *hMC_muPt2 = new TH1D("hMC_muPt2","",40,0,120);

   TH1D *hMC_muPtRatio1 = new TH1D("hMC_muPtRatio1","",40,0.8,1.2);
   TH1D *hMC_muPtRatio2 = new TH1D("hMC_muPtRatio2","",40,0.8,1.2);

   TH1D *hMC_genMuPt1 = new TH1D("hMC_genMuPt1","",40,0,120);
   TH1D *hMC_genMuPt2 = new TH1D("hMC_genMuPt2","",40,0,120);

   TH2D *hData_muPt12 = new TH2D("hData_muPt12","",40,0,120,40,0,120);
   TH2D *hMC_muPt12 = new TH2D("hMC_muPt12","",40,0,120,40,0,120);
   TH2D *hMC_genMuPt12 = new TH2D("hMC_genMuPt12","",40,0,120,40,0,120);

   TH2D *hData_muPt12_dr0 = new TH2D("hData_muPt12_dr0","",40,0,120,40,0,120);
   TH2D *hMC_muPt12_dr0 = new TH2D("hMC_muPt12_dr0","",40,0,120,40,0,120);
   TH2D *hMC_genMuPt12_dr0 = new TH2D("hMC_genMuPt12_dr0","",40,0,120,40,0,120);

   TH2D *hData_zMDR = new TH2D("hData_zMDR","",40,0,5,40,0,150);
   TH2D *hMC_zMDR = new TH2D("hMC_zMDR","",40,0,5,40,0,150);
   TH2D *hMC_genzMDR = new TH2D("hMC_genzMDR","",40,0,5,40,0,150);

   TH2D *hData_zPtM = new TH2D("hData_zPtM","",40,0,150,40,0,200);
   TH2D *hMC_zPtM = new TH2D("hMC_zPtM","",40,0,150,40,0,200);
   TH2D *hMC_genzPtM = new TH2D("hMC_genzPtM","",40,0,150,40,0,200);

   TH1D *hData_muDeta = new TH1D("hData_muDeta","",40,-6.284,6.284);
   TH1D *hMC_muDeta = new TH1D("hMC_muDeta","",40,-6.284,6.284);
   TH1D *hMC_genMuDeta = new TH1D("hMC_genMuDeta","",40,-6.284,6.284);
   TH1D *hData_muDphi = new TH1D("hData_muDphi","",40,-3.1415926,3.1415926);
   TH1D *hMC_muDphi = new TH1D("hMC_muDphi","",40,-3.1415926,3.1415926);
   TH1D *hMC_genMuDphi = new TH1D("hMC_genMuDphi","",40,-3.1415926,3.1415926);

   TH1D *hData_muDR = new TH1D("hData_muDR","",40,0,7.025);
   TH1D *hMC_muDR = new TH1D("hMC_muDR","",40,0,7.025);
   TH1D *hMC_genMuDR = new TH1D("hMC_genMuDR","",40,0,7.025);
   TH1D *hData_muDphiS = new TH1D("hData_muDphiS","",40,-1,1);
   TH1D *hMC_muDphiS = new TH1D("hMC_muDphiS","",40,-1,1);
   TH1D *hMC_genMuDphiS = new TH1D("hMC_genMuDphiS","",40,-1,1);

   TCut zMassRange = "zMass>60";
   TCut zPtRange = Form("zPt>%f&&zPt<%f",ptL,ptH);
   TCut muPt1Range = Form("max(genMuPt1,genMuPt2)>%f&&max(genMuPt1,genMuPt2)<%f",ptL,ptH);
   TCut muPt2Range = Form("min(genMuPt1,genMuPt2)>%f&&min(genMuPt1,genMuPt2)<%f",ptL,ptH);
   TCut hiHFRange = Form("hiHF<=%.4f&&hiHF>%.4f",hf_diff[centL],hf_diff[centH]);
   
   tData->Draw("zMass>>hData",zMassRange&&zPtRange&&hiHFRange);
   hData->Sumw2();

   tMC->Draw("zMass>>hMC",zMassRange&&zPtRange&&hiHFRange);
   hMC->Sumw2();

   tData->Draw("zEta>>hData_eta",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("zEta>>hMC_eta",zMassRange&&zPtRange&&hiHFRange);

   tData->Draw("zPhi>>hData_phi",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("zPhi>>hMC_phi",zMassRange&&zPtRange&&hiHFRange);

   tData->Draw("zPt>>hData_pt",hiHFRange);
   tMC->Draw("zPt>>hMC_pt",hiHFRange);

   tData->Draw("muPt1>>hData_muPt1",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("muPt1>>hMC_muPt1",zMassRange&&zPtRange&&hiHFRange);

   tData->Draw("muPt2>>hData_muPt2",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("muPt2>>hMC_muPt2",zMassRange&&zPtRange&&hiHFRange);

   tMC->Draw("max(genMuPt1,genMuPt2)>>hMC_genMuPt1",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("min(genMuPt1,genMuPt2)>>hMC_genMuPt2",zMassRange&&zPtRange&&hiHFRange);
   
   tMC->Draw("muPt1/max(genMuPt1,genMuPt2)>>hMC_muPtRatio1",zMassRange&&muPt1Range&&hiHFRange);
   tMC->Draw("muPt2/min(genMuPt1,genMuPt2)>>hMC_muPtRatio2",zMassRange&&muPt2Range&&hiHFRange);

   tData->Draw("muPt1:muPt2>>hData_muPt12",zPtRange&&hiHFRange);
   tMC->Draw("muPt1:muPt2>>hMC_muPt12",zPtRange&&hiHFRange);
   tMC->Draw("max(genMuPt1,genMuPt2):min(genMuPt1,genMuPt2)>>hMC_genMuPt12",zPtRange&&hiHFRange);

   tData->Draw("muPt1:muPt2>>hData_muPt12_dr0","muDR<1"&&zPtRange&&hiHFRange);
   tMC->Draw("muPt1:muPt2>>hMC_muPt12_dr0","muDR<1"&&zPtRange&&hiHFRange);
   tMC->Draw("max(genMuPt1,genMuPt2):min(genMuPt1,genMuPt2)>>hMC_genMuPt12_dr0","muDR<1"&&zPtRange&&hiHFRange);

   tData->Draw("zMass:muDR>>hData_zMDR",zPtRange&&hiHFRange);
   tMC->Draw("zMass:muDR>>hMC_zMDR",zPtRange&&hiHFRange);
   tMC->Draw("genZMass:genMuDR>>hMC_genzMDR",zPtRange&&hiHFRange);

   tData->Draw("zPt:zMass>>hData_zPtM",hiHFRange);
   tMC->Draw("zPt:zMass>>hMC_zPtM",hiHFRange);
   tMC->Draw("genZPt:genZMass>>hMC_genzPtM",hiHFRange);

   tData->Draw("muDeta>>hData_muDeta",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("muDeta>>hMC_muDeta",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("genMuDeta>>hMC_genMuDeta",zMassRange&&zPtRange&&hiHFRange);

   tData->Draw("muDphi>>hData_muDphi",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("muDphi>>hMC_muDphi",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("genMuDphi>>hMC_genMuDphi",zMassRange&&zPtRange&&hiHFRange);

   tData->Draw("muDR>>hData_muDR",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("muDR>>hMC_muDR",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("genMuDR>>hMC_genMuDR",zMassRange&&zPtRange&&hiHFRange);

   tData->Draw("muDphiS>>hData_muDphiS",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("muDphiS>>hMC_muDphiS",zMassRange&&zPtRange&&hiHFRange);
   tMC->Draw("genMuDphiS>>hMC_genMuDphiS",zMassRange&&zPtRange&&hiHFRange);

   //int countD = tData->GetEntries(zMassRange&&zPtRange&&hiHFRange);
   //std::cout<<"Data = "<<countD<<std::endl;
   //int countM = tMC->GetEntries(zMassRange&&zPtRange&&hiHFRange);
   //std::cout<<"MC = "<<countM<<std::endl;

   int countD = hData->GetEntries();
   std::cout<<"Data = "<<countD<<std::endl;
   int countM = hMC->GetEntries();
   std::cout<<"MC = "<<countM<<std::endl;

   
   hMC->SetMarkerStyle(24);
   //hMCSame->SetMarkerStyle(24);
   //hMCSame->SetLineColor(2);
   //hMCSame->SetMarkerColor(2);

   hData_eta->Sumw2();
   hMC_eta->Sumw2();
   hData_phi->Sumw2();
   hMC_phi->Sumw2();
   hData_pt->Sumw2();
   hMC_pt->Sumw2();

   hMC_muPtRatio1->Sumw2();
   hMC_muPtRatio2->Sumw2();

   hData_muPt1->Sumw2();
   hMC_muPt1->Sumw2();
   hData_muPt2->Sumw2();
   hMC_muPt2->Sumw2();

   hMC_genMuPt1->Sumw2();
   hMC_genMuPt2->Sumw2();

   hData_muDeta->Sumw2();
   hData_muDphi->Sumw2();
   hData_muDR->Sumw2();
   hData_muDphiS->Sumw2();
   hMC_muDeta->Sumw2();
   hMC_muDphi->Sumw2(); 
   hMC_muDR->Sumw2(); 
   hMC_muDphiS->Sumw2();

   hMC_genMuDeta->Sumw2();
   hMC_genMuDphi->Sumw2();
   hMC_genMuDR->Sumw2();
   hMC_genMuDphiS->Sumw2();

   hData_muPt12->Sumw2();
   hMC_muPt12->Sumw2();
   hMC_genMuPt12->Sumw2();

   hData_muPt12_dr0->Sumw2();
   hMC_muPt12_dr0->Sumw2();
   hMC_genMuPt12_dr0->Sumw2();

   hData_zMDR->Sumw2();
   hMC_zMDR->Sumw2();
   hMC_genzMDR->Sumw2();

   hData_zPtM->Sumw2();
   hMC_zPtM->Sumw2();
   hMC_genzPtM->Sumw2();
   
   hMC_eta->SetMarkerStyle(24);
   hMC_phi->SetMarkerStyle(24);
   hMC_pt->SetMarkerStyle(24);
   hMC_muPt1->SetMarkerStyle(24);
   hMC_muPt2->SetMarkerStyle(24);

   hMC_genMuPt1->SetMarkerStyle(25);
   hMC_genMuPt2->SetMarkerStyle(25);

   hMC_genMuDeta->SetMarkerStyle(25);
   hMC_genMuDphi->SetMarkerStyle(25);
   hMC_genMuDR->SetMarkerStyle(25);
   hMC_genMuDphiS->SetMarkerStyle(25);

   hMC_muDeta->SetMarkerStyle(24);
   hMC_muDphi->SetMarkerStyle(24);
   hMC_muDR->SetMarkerStyle(24);
   hMC_muDphiS->SetMarkerStyle(24);

   
//   hData->Draw("e");
//   hMC->Draw("e same");
 //  hDataSame->Draw("same");
//   hMCSame->Draw("same");

   //hData->SetMarkerStyle(kFullCircle);
   hData->SetMarkerColor(kBlack);
   //hData->SetMarkerSize(2);
   hData->SetLineColor(kBlack);
   
   //hMC->SetMarkerStyle(kFullSquare);
   hMC->SetMarkerColor(kRed);
   //hMC->SetMarkerSize(2);
   hMC->SetLineColor(kRed);

   hMC -> SetLineWidth(2);
   hData->SetLineWidth(2);

   hData_eta->SetMarkerColor(kBlack);
   hMC_eta->SetMarkerColor(kRed);
   hData_phi->SetMarkerColor(kBlack);
   hMC_phi->SetMarkerColor(kRed);
   hData_pt->SetMarkerColor(kBlack);
   hMC_pt->SetMarkerColor(kRed);

   hMC_muPtRatio1->SetMarkerColor(kBlack);
   hMC_muPtRatio2->SetMarkerColor(kBlack);

   hMC_muPtRatio1->SetLineColor(kBlack);
   hMC_muPtRatio2->SetLineColor(kBlack);

   hData_muPt1->SetMarkerColor(kBlack);
   //hMC_muPt1->SetMarkerColor(kRed);
   hData_muPt2->SetMarkerColor(kBlack);
   //hMC_muPt2->SetMarkerColor(kRed);

   hData_muDeta->SetMarkerColor(kBlack);
   hData_muDphi->SetMarkerColor(kBlack);
   hData_muDR->SetMarkerColor(kBlack);
   hData_muDphiS->SetMarkerColor(kBlack);

   hData_eta->SetLineColor(kBlack);
   hMC_eta->SetLineColor(kRed);
   hData_phi->SetLineColor(kBlack);
   hMC_phi->SetLineColor(kRed);
   hData_pt->SetLineColor(kBlack);
   hMC_pt->SetLineColor(kRed);

   hData_muPt1->SetLineColor(kBlack);
   //hMC_muPt1->SetLineColor(kRed);
   hData_muPt2->SetLineColor(kBlack);
   //hMC_muPt2->SetLineColor(kRed);

   hData_muDeta->SetLineColor(kBlack);
   hData_muDphi->SetLineColor(kBlack);
   hData_muDR->SetLineColor(kBlack);
   hData_muDphiS->SetLineColor(kBlack);

   //hMC_muDeta->SetLineColor(kRed);
   //hMC_muDphi->SetLineColor(kRed);
   //hMC_muDR->SetLineColor(kRed);
   //hMC_muDphiS->SetLineColor(kRed);

   hData_muPt1->SetLineWidth(2);
   hData_muPt2->SetLineWidth(2);

   hData_muDeta->SetLineWidth(2);
   hData_muDphi->SetLineWidth(2);
   hData_muDR->SetLineWidth(2);
   hData_muDphiS->SetLineWidth(2);

   hMC_muPtRatio1->SetLineWidth(2);
   hMC_muPtRatio2->SetLineWidth(2);

   hMC_genMuPt1->SetLineColor(TColor::GetColor("#004488"));
   hMC_muPt1->SetLineColor(TColor::GetColor("#994455"));
   hMC_genMuPt1->SetFillColor(TColor::GetColor("#6699CC"));
   hMC_muPt1->SetFillColor(TColor::GetColor("#EE99AA"));
   hMC_genMuPt1->SetFillStyle(3345); hMC_genMuPt1->SetLineWidth(3);
   hMC_muPt1->SetFillStyle(3354); hMC_muPt1->SetLineWidth(3);

   hMC_genMuPt2->SetLineColor(TColor::GetColor("#004488"));
   hMC_muPt2->SetLineColor(TColor::GetColor("#994455"));
   hMC_genMuPt2->SetFillColor(TColor::GetColor("#6699CC"));
   hMC_muPt2->SetFillColor(TColor::GetColor("#EE99AA"));
   hMC_genMuPt2->SetFillStyle(3345); hMC_genMuPt2->SetLineWidth(3);
   hMC_muPt2->SetFillStyle(3354); hMC_muPt2->SetLineWidth(3);

   hMC_genMuDeta->SetLineColor(TColor::GetColor("#004488"));
   hMC_muDeta->SetLineColor(TColor::GetColor("#994455"));
   hMC_genMuDeta->SetFillColor(TColor::GetColor("#6699CC"));
   hMC_muDeta->SetFillColor(TColor::GetColor("#EE99AA"));
   hMC_genMuDeta->SetFillStyle(3345); hMC_genMuDeta->SetLineWidth(3);
   hMC_muDeta->SetFillStyle(3354); hMC_muDeta->SetLineWidth(3);

   hMC_genMuDphi->SetLineColor(TColor::GetColor("#004488"));
   hMC_muDphi->SetLineColor(TColor::GetColor("#994455"));
   hMC_genMuDphi->SetFillColor(TColor::GetColor("#6699CC"));
   hMC_muDphi->SetFillColor(TColor::GetColor("#EE99AA"));
   hMC_genMuDphi->SetFillStyle(3345); hMC_genMuDphi->SetLineWidth(3);
   hMC_muDphi->SetFillStyle(3354); hMC_muDphi->SetLineWidth(3);

   hMC_genMuDR->SetLineColor(TColor::GetColor("#004488"));
   hMC_muDR->SetLineColor(TColor::GetColor("#994455"));
   hMC_genMuDR->SetFillColor(TColor::GetColor("#6699CC"));
   hMC_muDR->SetFillColor(TColor::GetColor("#EE99AA"));
   hMC_genMuDR->SetFillStyle(3345); hMC_genMuDR->SetLineWidth(3);
   hMC_muDR->SetFillStyle(3354); hMC_muDR->SetLineWidth(3);

   hMC_genMuDphiS->SetLineColor(TColor::GetColor("#004488"));
   hMC_muDphiS->SetLineColor(TColor::GetColor("#994455"));
   hMC_genMuDphiS->SetFillColor(TColor::GetColor("#6699CC"));
   hMC_muDphiS->SetFillColor(TColor::GetColor("#EE99AA"));
   hMC_genMuDphiS->SetFillStyle(3345); hMC_genMuDphiS->SetLineWidth(3);
   hMC_muDphiS->SetFillStyle(3354); hMC_muDphiS->SetLineWidth(3);

   //hMC->GetXaxis()->SetTitleSize(48);
   //hMC->GetXaxis()->SetTitleFont(43);

   //hData->Scale(1./hData->GetEntries());
   //hMC->Scale(1./hMC->GetEntries());

   hData->Scale(1./hData->Integral("width"));
   hMC->Scale(1./hMC->Integral("width"));

   double max1 = hMC->GetMaximum();
   double max2 = hData->GetMaximum();
   
   if(max1<max2) hData->Draw();
   else hMC->Draw();
   hMC->Draw("same");
   hData->Draw("same");

   hMC->SetXTitle("M_{#mu#mu} (GeV)");
   hData->SetXTitle("M_{#mu#mu} (GeV)");
   //hMC->SetLineColor(2);
   //hMC->SetMarkerColor(2);

   hMC->SetMinimum(0);
   hData->SetMinimum(0);

/*
   hData_eta->Scale(1./hData_eta->GetEntries());
   hData_phi->Scale(1./hData_phi->GetEntries());
   hData_pt->Scale(1./hData_pt->GetEntries());
   hMC_eta->Scale(1./hMC_eta->GetEntries());
   hMC_phi->Scale(1./hMC_phi->GetEntries());
   hMC_pt->Scale(1./hMC_pt->GetEntries());
*/

   hData_eta->Scale(1./hData_eta->Integral("width"));
   hData_phi->Scale(1./hData_phi->Integral("width"));
   hData_pt->Scale(1./hData_pt->Integral("width"));
   hMC_eta->Scale(1./hMC_eta->Integral("width"));
   hMC_phi->Scale(1./hMC_phi->Integral("width"));
   hMC_pt->Scale(1./hMC_pt->Integral("width"));

   hData_muPt1->Scale(1./hData_muPt1->Integral("width"));
   hMC_muPt1->Scale(1./hMC_muPt1->Integral("width"));
   hData_muPt2->Scale(1./hData_muPt2->Integral("width"));
   hMC_muPt2->Scale(1./hMC_muPt2->Integral("width"));

   hMC_genMuPt1->Scale(1./hMC_genMuPt1->Integral("width"));
   hMC_genMuPt2->Scale(1./hMC_genMuPt2->Integral("width"));

   hData_muDeta->Scale(1./hData_muDeta->Integral("width"));
   hData_muDphi->Scale(1./hData_muDphi->Integral("width"));
   hData_muDR->Scale(1./hData_muDR->Integral("width"));
   hData_muDphiS->Scale(1./hData_muDphiS->Integral("width"));
   hMC_muDeta->Scale(1./hMC_muDeta->Integral("width"));
   hMC_muDphi->Scale(1./hMC_muDphi->Integral("width"));
   hMC_muDR->Scale(1./hMC_muDR->Integral("width"));
   hMC_muDphiS->Scale(1./hMC_muDphiS->Integral("width"));
   hMC_genMuDeta->Scale(1./hMC_genMuDeta->Integral("width"));
   hMC_genMuDphi->Scale(1./hMC_genMuDphi->Integral("width"));
   hMC_genMuDR->Scale(1./hMC_genMuDR->Integral("width"));
   hMC_genMuDphiS->Scale(1./hMC_genMuDphiS->Integral("width"));

   hData_muPt12->Scale(1./hData_muPt12->Integral("width"));
   hMC_muPt12->Scale(1./hMC_muPt12->Integral("width"));
   hMC_genMuPt12->Scale(1./hMC_genMuPt12->Integral("width"));

   hData_muPt12_dr0->Scale(1./hData_muPt12_dr0->Integral("width"));
   hMC_muPt12_dr0->Scale(1./hMC_muPt12_dr0->Integral("width"));
   hMC_genMuPt12_dr0->Scale(1./hMC_genMuPt12_dr0->Integral("width"));

   hData_zMDR->Scale(1./hData_zMDR->Integral("width"));
   hMC_zMDR->Scale(1./hMC_zMDR->Integral("width"));
   hMC_genzMDR->Scale(1./hMC_genzMDR->Integral("width"));
   hData_zPtM->Scale(1./hData_zPtM->Integral("width"));
   hMC_zPtM->Scale(1./hMC_zPtM->Integral("width"));
   hMC_genzPtM->Scale(1./hMC_genzPtM->Integral("width"));

   hMC_muPtRatio1->Scale(1./hMC_muPtRatio1->Integral("width"));
   hMC_muPtRatio2->Scale(1./hMC_muPtRatio2->Integral("width"));

   ////style();
   
//   TF1 *f = new TF1("f","[0]+[1]*x+[2]*TMath::BreitWigner(x, [3], [4])");
   TF1 *f = new TF1("f",bwfun,81.2,101.2,7);
   f->SetParNames("Mean","Width","Area","GSigma","BkgArea","ExpCont","ExpShift");
   f->SetParameters(91,4.6,3,0.1,0,0,0,0);
   f->SetLineStyle(2);
   //f->SetLineColor(1);
   f->SetLineColor(kBlack);
   
 //  TF1 *f2 = new TF1("f2","[0]+[1]*x+[2]*TMath::BreitWigner(x, [3], [4])");
//   f2->SetParameters(0,0,1,91.2,1);

   TF1 *f2 = new TF1("f2",bwfun,81.2,101.2,7);
   f2->SetParNames("Mean","Width","Area","GSigma","BkgArea","ExpCont","ExpShift");
   f2->SetParameters(91,4.6,3,0.1,0,0,0,0);
   f2->SetLineStyle(2);
   //f2->SetLineColor(2);
   f2->SetLineColor(kRed);
   
   hData->Fit("f","LL");
   hData->Fit("f","");
   hData->Fit("f","LL m");
   hData->Fit("f","");
   hData->Fit("f","LL m");
   hData->Fit("f","");
   hData->Fit("f","");
   hData->Fit("f","");
   hData->Fit("f","");
   hData->Fit("f","");
 
   hMC->Fit("f2","LL m");
   hMC->Fit("f2","");
   hMC->Fit("f2","LL");
   hMC->Fit("f2","");
   hMC->Fit("f2","LL m");
   hMC->Fit("f2","");
   hMC->Fit("f2","");
   hMC->Fit("f2","");
   hMC->Fit("f2","");
   hMC->Fit("f2","");

/*
   std::cout<<"UnbinnedFit 1"<<std::endl;
   tData->UnbinnedFit("f", "zMass", zPtRange&&hiHFRange,"");
   //tMC->UnbinnedFit("f2", "zMass", zPtRange&&hiHFRange,"");

   std::cout<<"UnbinnedFit 2"<<std::endl;
   tData->UnbinnedFit("f", "zMass", zPtRange&&hiHFRange,"E");
   //tMC->UnbinnedFit("f2", "zMass", zPtRange&&hiHFRange,"E");

   std::cout<<"UnbinnedFit 3"<<std::endl;
   tData->UnbinnedFit("f", "zMass", zPtRange&&hiHFRange,"E M");
   //tMC->UnbinnedFit("f2", "zMass", zPtRange&&hiHFRange,"E M");
  */

   h->SetBinContent(1,f->GetParameter(0));
   h->SetBinContent(2,f->GetParError(0));
   h->SetBinContent(3,f2->GetParameter(0));
   h->SetBinContent(4,f2->GetParError(0));

   h->SetBinContent(5,f->GetParameter(1));
   h->SetBinContent(6,f->GetParError(1));
   h->SetBinContent(7,f2->GetParameter(1));
   h->SetBinContent(8,f2->GetParError(1));

   h->SetBinContent(9, hData->GetMean());
   h->SetBinContent(10,hMC->GetMean());
   h->SetBinContent(11,hData->GetStdDev());
   h->SetBinContent(12,hMC->GetStdDev());

   TLegend leg(0.58,0.78,0.98,0.9);
   leg.AddEntry(hMC ,"Monte Carlo: DYLL","lep");
   leg.AddEntry(hData ,Form("Data: %s",typeofdatatext),"lep");
   leg.SetFillColorAlpha(kWhite,0);
   leg.SetLineColor(kBlack);
   leg.SetLineWidth(1);
   leg.Draw();

   TLegend legMuPt(0.58,0.68,0.98,0.9);
   legMuPt.AddEntry(hMC_genMuPt1 ,"Monte Carlo: GEN level","lf");
   legMuPt.AddEntry(hMC_muPt1 ,"Monte Carlo: RECO","lf");
   legMuPt.AddEntry(hData_muPt1 ,Form("Data: %s",typeofdatatext),"lep");
   legMuPt.SetFillColorAlpha(kWhite,0);
   legMuPt.SetLineColor(kBlack);
   legMuPt.SetLineWidth(1);
   //legMuPt.Draw();

   //TLatex *pt = new TLatex(0.18,0.82,Form("%d < Centrality < %d",centL,centH));
   TLatex *pt = new TLatex(0.18,0.82,Form("%.0f %%< Centrality < %.0f %%",cent_diff[centL],cent_diff[centH]));
   pt->SetTextFont(42);
   //std::cout<<"TextSize = "<<pt->GetTextSize()<<std::endl;
   pt->SetTextSize(0.03);
   pt->SetNDC(kTRUE);
   pt->Draw();
   TLatex *pt2 = new TLatex(0.18,0.88,Form("%.1f < Z p_{T} < %.1f",ptL,ptH));
   pt2->SetTextFont(42);
   //std::cout<<"TextSize2 = "<<pt2->GetTextSize()<<std::endl;
   pt2->SetTextSize(0.03);
   pt2->SetNDC(kTRUE);
   pt2->Draw();

   TLatex *ptp = new TLatex(0.18,0.85,Form("%.0f %%< Centrality < %.0f %%",cent_diff[centL],cent_diff[centH]));
   ptp->SetTextFont(42);
   ptp->SetTextSize(0.025);
   ptp->SetNDC(kTRUE);
   //ptp->Draw();
   TLatex *ptp2 = new TLatex(0.18,0.89,Form("%.1f < Z p_{T} < %.1f",ptL,ptH));
   ptp2->SetTextFont(42);
   ptp2->SetTextSize(0.025);
   ptp2->SetNDC(kTRUE);
   //ptp2->Draw();

   TLatex *ptpp1 = new TLatex(0.18,0.88,Form("%.1f < p_{T,#mu_{1}} < %.1f",ptL,ptH));
   ptpp1->SetTextFont(42);
   ptpp1->SetTextSize(0.03);
   ptpp1->SetNDC(kTRUE);

   TLatex *ptpp2 = new TLatex(0.18,0.88,Form("%.1f < p_{T,#mu_{2}} < %.1f",ptL,ptH));
   ptpp2->SetTextFont(42);
   ptpp2->SetTextSize(0.03);
   ptpp2->SetNDC(kTRUE);

   TLatex *ptN = new TLatex(0.6,0.97,Form("N_{MC} = %d, N_{Data} = %d",countM,countD));
   ptN->SetTextFont(42);
   ptN->SetTextSize(0.03);
   ptN->SetNDC(kTRUE);
   ptN->Draw();

   gSystem->Exec(Form("mkdir -p figs/mass/%s",typeofdata));

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   ////style();
   
   max1 = hMC_eta->GetMaximum();
   max2 = hData_eta->GetMaximum();
   
   if(max1<max2) hData_eta->Draw();
   else hMC_eta->Draw();
   hMC_eta->Draw("same");
   hData_eta->Draw("same");

   hMC_eta->SetXTitle("#eta_{Z}");
   hData_eta->SetXTitle("#eta_{Z}");


   leg.Draw();
   pt->Draw();
   pt2->Draw();
   hMC_eta->SetMinimum(0);
   hData_eta->SetMinimum(0);
   //if (ptL==0&&ptH==200) hMC_eta->SetMaximum(0.1); 
   //else if(ptL==0&&ptH==20) hMC_eta->SetMaximum(0.1); 
   //else if(ptL==80&&ptH==100) hMC_eta->SetMaximum(0.2); 
   //else if(ptL==0&&ptH==2000) hMC_eta->SetMaximum(0.1); 
   //else hMC_eta->SetMaximum(0.15);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_eta.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   ////style();

   max1 = hMC_phi->GetMaximum();
   max2 = hData_phi->GetMaximum();
   
   if(max1<max2) hData_phi->Draw();
   else hMC_phi->Draw();
   hMC_phi->Draw("same");
   hData_phi->Draw("same");

   hData_phi->SetXTitle("#phi_{Z}");
   hMC_phi->SetXTitle("#phi_{Z}");

   leg.Draw();
   pt->Draw();
   pt2->Draw();
   hMC_phi->SetMinimum(0);
   hData_phi->SetMinimum(0);
   //hMC_phi->SetMaximum(0.08);
   //if(ptL==80&&ptH==100) hMC_phi->SetMaximum(0.1);
   //else if(ptL==0&&ptH==2000) hMC_phi->SetMaximum(0.05);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_phi.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   max1 = hMC_muPt1->GetMaximum();
   max2 = hData_muPt1->GetMaximum();
   double max3 = hMC_genMuPt1->GetMaximum();
   
   if(max1<max2&&max3<max2) hData_muPt1->Draw();
   else if(max1<max3&&max2<max3) hMC_genMuPt1->Draw("hist");
   else hMC_muPt1->Draw("hist");

   hMC_genMuPt1->Draw("hist same");
   hMC_muPt1->Draw("hist same");
   hData_muPt1->Draw("same");

   hData_muPt1->SetXTitle("p_{T,#mu_{1}} (GeV)");
   hMC_muPt1->SetXTitle("p_{T,#mu_{1}} (GeV)");
   hMC_genMuPt1->SetXTitle("p_{T,#mu_{1}} (GeV)");

   legMuPt.Draw();
   ptp->Draw();
   ptp2->Draw();
   hMC_muPt1->SetMinimum(0);
   hData_muPt1->SetMinimum(0);
   hMC_genMuPt1->SetMinimum(0);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muPt1.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   max1 = hMC_muPt2->GetMaximum();
   max2 = hData_muPt2->GetMaximum();
   max3 = hMC_genMuPt2->GetMaximum();

   if(max1<max2&&max3<max2) hData_muPt2->Draw();
   else if(max1<max3&&max2<max3) hMC_genMuPt2->Draw("hist");
   else hMC_muPt2->Draw("hist");
   
   hMC_genMuPt2->Draw("hist same");
   hMC_muPt2->Draw("hist same");
   hData_muPt2->Draw("same");
   
   hData_muPt2->SetXTitle("p_{T,#mu_{2}} (GeV)");
   hMC_muPt2->SetXTitle("p_{T,#mu_{2}} (GeV)");
   hMC_genMuPt2->SetXTitle("p_{T,#mu_{2}} (GeV)");

   legMuPt.Draw();
   ptp->Draw();
   ptp2->Draw();
   hMC_muPt2->SetMinimum(0);
   hData_muPt2->SetMinimum(0);
   hMC_genMuPt2->SetMinimum(0);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muPt2.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   max1 = hMC_muDeta->GetMaximum();
   max2 = hData_muDeta->GetMaximum();
   max3 = hMC_genMuDeta->GetMaximum();

   if(max1<max2&&max3<max2) hData_muDeta->Draw();
   else if(max1<max3&&max2<max3) hMC_genMuDeta->Draw("hist");
   else hMC_muDeta->Draw("hist");
   
   hMC_genMuDeta->Draw("hist same");
   hMC_muDeta->Draw("hist same");
   hData_muDeta->Draw("same");
   
   hData_muDeta->SetXTitle("#Delta#eta_{#mu#mu}");
   hMC_muDeta->SetXTitle("#Delta#eta_{#mu#mu}");
   hMC_genMuDeta->SetXTitle("#Delta#eta_{#mu#mu}");

   legMuPt.Draw();
   ptp->Draw();
   ptp2->Draw();
   hMC_muDeta->SetMinimum(0);
   hData_muDeta->SetMinimum(0);
   hMC_genMuDeta->SetMinimum(0);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muDeta.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   max1 = hMC_muDphi->GetMaximum();
   max2 = hData_muDphi->GetMaximum();
   max3 = hMC_genMuDphi->GetMaximum();

   if(max1<max2&&max3<max2) hData_muDphi->Draw();
   else if(max1<max3&&max2<max3) hMC_genMuDphi->Draw("hist");
   else hMC_muDphi->Draw("hist");
   
   hMC_genMuDphi->Draw("hist same");
   hMC_muDphi->Draw("hist same");
   hData_muDphi->Draw("same");
   
   hData_muDphi->SetXTitle("#Delta#phi_{#mu#mu}");
   hMC_muDphi->SetXTitle("#Delta#phi_{#mu#mu}");
   hMC_genMuDphi->SetXTitle("#Delta#phi_{#mu#mu}");

   legMuPt.Draw();
   ptp->Draw();
   ptp2->Draw();
   hMC_muDphi->SetMinimum(0);
   hData_muDphi->SetMinimum(0);
   hMC_genMuDphi->SetMinimum(0);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muDphi.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   max1 = hMC_muDR->GetMaximum();
   max2 = hData_muDR->GetMaximum();
   max3 = hMC_genMuDR->GetMaximum();

   if(max1<max2&&max3<max2) hData_muDR->Draw();
   else if(max1<max3&&max2<max3) hMC_genMuDR->Draw("hist");
   else hMC_muDR->Draw("hist");
   
   hMC_genMuDR->Draw("hist same");
   hMC_muDR->Draw("hist same");
   hData_muDR->Draw("same");
   
   hData_muDR->SetXTitle("#DeltaR_{#mu#mu}");
   hMC_muDR->SetXTitle("#DeltaR_{#mu#mu}");
   hMC_genMuDR->SetXTitle("#DeltaR_{#mu#mu}");

   legMuPt.Draw();
   ptp->Draw();
   ptp2->Draw();
   hMC_muDR->SetMinimum(0);
   hData_muDR->SetMinimum(0);
   hMC_genMuDR->SetMinimum(0);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muDR.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   max1 = hMC_muDphiS->GetMaximum();
   max2 = hData_muDphiS->GetMaximum();
   max3 = hMC_genMuDphiS->GetMaximum();

   if(max1<max2&&max3<max2) hData_muDphiS->Draw();
   else if(max1<max3&&max2<max3) hMC_genMuDphiS->Draw("hist");
   else hMC_muDphiS->Draw("hist");
   
   hMC_genMuDphiS->Draw("hist same");
   hMC_muDphiS->Draw("hist same");
   hData_muDphiS->Draw("same");
   
   hData_muDphiS->SetXTitle("#Delta#phi*_{#mu#mu}");
   hMC_muDphiS->SetXTitle("#Delta#phi*_{#mu#mu}");
   hMC_genMuDphiS->SetXTitle("#Delta#phi*_{#mu#mu}");

   legMuPt.Draw();
   ptp->Draw();
   ptp2->Draw();
   hMC_muDphiS->SetMinimum(0);
   hData_muDphiS->SetMinimum(0);
   hMC_genMuDphiS->SetMinimum(0);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muDphiS.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();



   if(ptL==0&&(ptH==200||ptH==2000)){
      ////style();

      max1 = hMC_pt->GetMaximum();
      max2 = hData_pt->GetMaximum();
      
      if(max1<max2) hData_pt->Draw();
      else hMC_pt->Draw();
      hMC_pt->Draw("same");
      hData_pt->Draw("same");
   
      hData_pt->SetXTitle("Z p_{T} (GeV)");
      hMC_pt->SetXTitle("Z p_{T} (GeV)");

      hMC_pt->GetXaxis()->SetLimits(0,200);
      hData_pt->GetXaxis()->SetLimits(0,200);

      hMC_pt->SetMinimum(0);
      hData_pt->SetMinimum(0);
      //hMC_pt->SetMaximum(0.3);
      leg.Draw();
      pt->Draw();
      pt->SetY(0.88);
      //pt2->Draw();
   
      ptN->Draw();

      c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_pt.png",typeofdata,typeofdata,cent_diff[centL],cent_diff[centH])); 
      hMC_pt->SetMinimum(0.00001);
      hData_pt->SetMinimum(0.00001);
      c->SetLogy(1);
      //hMC_pt->Draw("same");
      //hMC_pt->Draw("axis same");

      c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_pt_log.png",typeofdata,typeofdata,cent_diff[centL],cent_diff[centH])); 
      c->SetLogy(0);
      c->Clear();
      pt->SetY(0.82);
   }

   c->SetCanvasSize(1400,800);

   TLatex *ptRm1 = new TLatex(0.18,0.76,Form("Mean = %.3f",hMC_muPtRatio1->GetMean()));
   TLatex *ptRs1 = new TLatex(0.18,0.70,Form("Std Dev = %.3f",hMC_muPtRatio1->GetStdDev()));
   TLatex *ptRm2 = new TLatex(0.18,0.76,Form("Mean = %.3f",hMC_muPtRatio2->GetMean()));
   TLatex *ptRs2 = new TLatex(0.18,0.70,Form("Std Dev = %.3f",hMC_muPtRatio2->GetStdDev()));

   ptRm1->SetTextFont(42);
   ptRm1->SetTextSize(0.03);
   ptRm1->SetNDC(kTRUE);

   ptRs1->SetTextFont(42);
   ptRs1->SetTextSize(0.03);
   ptRs1->SetNDC(kTRUE);

   ptRm2->SetTextFont(42);
   ptRm2->SetTextSize(0.03);
   ptRm2->SetNDC(kTRUE);

   ptRs2->SetTextFont(42);
   ptRs2->SetTextSize(0.03);
   ptRs2->SetNDC(kTRUE);

   c->Divide(2);
   hMC_muPtRatio1->SetMinimum(0);
   hMC_muPtRatio2->SetMinimum(0);

   c->cd(1);
   hMC_muPtRatio1->Draw();
   hMC_muPtRatio1->SetXTitle("Reco p_{T,#mu_{1}}/Gen p_{T,#mu_{1}} (GeV)");

   pt->Draw();
   ptpp1->Draw();
   ptRm1->Draw();
   ptRs1->Draw();

   c->cd(2);
   hMC_muPtRatio2->Draw();
   hMC_muPtRatio2->SetXTitle("Reco p_{T,#mu_{2}}/Gen p_{T,#mu_{2}} (GeV)");

   pt->Draw();
   ptpp2->Draw();
   ptRm2->Draw();
   ptRs2->Draw();

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muPtRatio.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 

   c->Clear();

   c->SetCanvasSize(2100,800);
   c->Divide(3);

   hData_muPt12->SetMinimum(0);
   hMC_muPt12->SetMinimum(0);
   hMC_genMuPt12->SetMinimum(0);

   c->cd(1);

   hData_muPt12->Draw("COLZ");
   hData_muPt12->GetYaxis()->SetTitle("Data p_{T,#mu_{1}} (GeV)");
   hData_muPt12->GetXaxis()->SetTitle("Data p_{T,#mu_{2}} (GeV)");
   hData_muPt12->GetXaxis()->SetNdivisions(50205,kFALSE);
   
   ptp->Draw();
   ptp2->Draw();

   c->cd(2);

   hMC_muPt12->Draw("COLZ");
   hMC_muPt12->GetYaxis()->SetTitle("MC RECO p_{T,#mu_{1}} (GeV)");
   hMC_muPt12->GetXaxis()->SetTitle("MC RECO p_{T,#mu_{2}} (GeV)");
   hMC_muPt12->GetXaxis()->SetNdivisions(50205,kFALSE);
   
   c->cd(3);

   hMC_genMuPt12->Draw("COLZ");
   hMC_genMuPt12->GetYaxis()->SetTitle("MC GEN p_{T,#mu_{1}} (GeV)");
   hMC_genMuPt12->GetXaxis()->SetTitle("MC GEN p_{T,#mu_{2}} (GeV)");
   hMC_genMuPt12->GetXaxis()->SetNdivisions(50205,kFALSE);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muPt12.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   TLatex *ptdr0 = new TLatex(0.18,0.82,"#DeltaR_{#mu#mu}<1");
   ptdr0->SetTextFont(42);
   ptdr0->SetTextSize(0.025);
   ptdr0->SetNDC(kTRUE);

   c->Divide(3);

   hData_muPt12_dr0->SetMinimum(0);
   hMC_muPt12_dr0->SetMinimum(0);
   hMC_genMuPt12_dr0->SetMinimum(0);

   c->cd(1);

   hData_muPt12_dr0->Draw("COLZ");
   hData_muPt12_dr0->GetYaxis()->SetTitle("Data p_{T,#mu_{1}} (GeV)");
   hData_muPt12_dr0->GetXaxis()->SetTitle("Data p_{T,#mu_{2}} (GeV)");
   hData_muPt12_dr0->GetXaxis()->SetNdivisions(50205,kFALSE);
   
   ptp->Draw();
   ptp2->Draw();
   ptdr0->Draw();

   c->cd(2);

   hMC_muPt12_dr0->Draw("COLZ");
   hMC_muPt12_dr0->GetYaxis()->SetTitle("MC RECO p_{T,#mu_{1}} (GeV)");
   hMC_muPt12_dr0->GetXaxis()->SetTitle("MC RECO p_{T,#mu_{2}} (GeV)");
   hMC_muPt12_dr0->GetXaxis()->SetNdivisions(50205,kFALSE);
   
   c->cd(3);

   hMC_genMuPt12_dr0->Draw("COLZ");
   hMC_genMuPt12_dr0->GetYaxis()->SetTitle("MC GEN p_{T,#mu_{1}} (GeV)");
   hMC_genMuPt12_dr0->GetXaxis()->SetTitle("MC GEN p_{T,#mu_{2}} (GeV)");
   hMC_genMuPt12_dr0->GetXaxis()->SetNdivisions(50205,kFALSE);

   //ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muPt12_dr0.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();


   c->Divide(3);

   hData_zMDR->SetMinimum(0);
   hMC_zMDR->SetMinimum(0);
   hMC_genzMDR->SetMinimum(0);

   c->cd(1);

   hData_zMDR->Draw("COLZ");
   hData_zMDR->GetYaxis()->SetTitle("Data M_{#mu#mu} (GeV)");
   hData_zMDR->GetXaxis()->SetTitle("Data #DeltaR_{#mu#mu}");
   hData_zMDR->GetXaxis()->SetNdivisions(50205,kFALSE);
   
   ptp->Draw();
   ptp2->Draw();

   c->cd(2);

   hMC_zMDR->Draw("COLZ");
   hMC_zMDR->GetYaxis()->SetTitle("MC RECO M_{#mu#mu} (GeV)");
   hMC_zMDR->GetXaxis()->SetTitle("MC RECO #DeltaR_{#mu#mu}");
   hMC_zMDR->GetXaxis()->SetNdivisions(50205,kFALSE);
   
   c->cd(3);

   hMC_genzMDR->Draw("COLZ");
   hMC_genzMDR->GetYaxis()->SetTitle("MC GEN M_{#mu#mu} (GeV)");
   hMC_genzMDR->GetXaxis()->SetTitle("MC GEN #DeltaR_{#mu#mu}");
   hMC_genzMDR->GetXaxis()->SetNdivisions(50205,kFALSE);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_zMDR.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   if(ptL==0&&(ptH==200||ptH==2000)){
      c->Divide(3);
   
      //hData_zPtM->SetMinimum(0);
      //hMC_zPtM->SetMinimum(0);
      //hMC_genzPtM->SetMinimum(0);

      c->cd(1);
   
      hData_zPtM->Draw("COLZ");
      hData_zPtM->GetYaxis()->SetTitle("Data Z p_{T} (GeV)");
      hData_zPtM->GetXaxis()->SetTitle("Data M_{#mu#mu} (GeV)");
      hData_zPtM->GetXaxis()->SetNdivisions(50205,kFALSE);
      
      ptp->Draw();
      //ptp2->Draw();
   
      c->cd(2);
   
      hMC_zPtM->Draw("COLZ");
      hMC_zPtM->GetYaxis()->SetTitle("MC RECO Z p_{T} (GeV)");
      hMC_zPtM->GetXaxis()->SetTitle("MC RECO M_{#mu#mu} (GeV)");
      hMC_zPtM->GetXaxis()->SetNdivisions(50205,kFALSE);
      
      c->cd(3);
   
      hMC_genzPtM->Draw("COLZ");
      hMC_genzPtM->GetYaxis()->SetTitle("MC GEN Z p_{T} (GeV)");
      hMC_genzPtM->GetXaxis()->SetTitle("MC GEN M_{#mu#mu} (GeV)");
      hMC_genzPtM->GetXaxis()->SetNdivisions(50205,kFALSE);
   
      ptN->Draw();
   
      c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_zPtM.png",typeofdata,typeofdata,cent_diff[centL],cent_diff[centH])); 
      c->Clear();
   }

   //hData->Reset();
   //hMC->Reset();

   //hData_eta->Reset();
   //hMC_eta->Reset();
   //hData_phi->Reset();
   //hMC_phi->Reset();

   //hData_pt->Reset();
   //hMC_pt->Reset();

   return h;

/* Seems useless
   delete hData; hData=NULL;
   delete hMC; hMC=NULL;
   delete hData_eta; hData_eta=NULL;
   delete hMC_eta; hMC_eta=NULL;
   delete hData_phi; hData_phi=NULL;
   delete hMC_phi; hMC_phi=NULL;
   delete hData_pt; hData_pt=NULL;
   delete hMC_pt; hMC_pt=NULL;
   delete hMC_muPt1; hMC_muPt1=NULL;
   delete hData_muPt1; hData_muPt1=NULL;
   delete hMC_muPt2; hMC_muPt2=NULL;
   delete hData_muPt2; hData_muPt2=NULL;
   delete hMC_genMuPt1; hMC_genMuPt1=NULL;
   delete hMC_genMuPt2; hMC_genMuPt2=NULL;

   delete gROOT->FindObject("hData");
   delete gROOT->FindObject("hMC");
   delete gROOT->FindObject("hData_eta");
   delete gROOT->FindObject("hMC_eta");
   delete gROOT->FindObject("hData_phi");
   delete gROOT->FindObject("hMC_phi");
   delete gROOT->FindObject("hData_pt");
   delete gROOT->FindObject("hMC_pt");
   delete gROOT->FindObject("hMC_muPt1");
   delete gROOT->FindObject("hData_muPt1");
   delete gROOT->FindObject("hMC_muPt2");
   delete gROOT->FindObject("hData_muPt2");
   delete gROOT->FindObject("hMC_genMuPt1");
   delete gROOT->FindObject("hMC_genMuPt2");
*/
}

void loop()
{
   TH1D *hDataMass = new TH1D("hDataMass","",5,0,100);
   TH1D *hMCMass = new TH1D("hMCMass","",5,0,100);

   TH1D *hDataWidth = new TH1D("hDataWidth","",5,0,100);
   TH1D *hMCWidth = new TH1D("hMCWidth","",5,0,100);

   TH1D *hDataMass_1 = new TH1D("hDataMass_1","",5,0,100);
   TH1D *hMCMass_1 = new TH1D("hMCMass_1","",5,0,100);

   TH1D *hDataWidth_1 = new TH1D("hDataWidth_1","",5,0,100);
   TH1D *hMCWidth_1 = new TH1D("hMCWidth_1","",5,0,100);
   
   for (int i=1;i<=hDataMass->GetNbinsX();i++)
   {
      TH1D *h = ZmassAnalysis_single(hDataMass->GetBinLowEdge(i),hDataMass->GetBinLowEdge(i+1));
      hDataMass->SetBinContent(i,h->GetBinContent(1));
      hDataMass->SetBinError(i,h->GetBinContent(2));
      hMCMass->SetBinContent(i,h->GetBinContent(3));
      hMCMass->SetBinError(i,h->GetBinContent(4));

      hDataWidth->SetBinContent(i,h->GetBinContent(5));
      hDataWidth->SetBinError(i,h->GetBinContent(6));
      hMCWidth->SetBinContent(i,h->GetBinContent(7));
      hMCWidth->SetBinError(i,h->GetBinContent(8));

      hDataMass_1->SetBinContent(i,h->GetBinContent(9));
      hMCMass_1->SetBinContent(i,h->GetBinContent(10));
      hDataWidth_1->SetBinContent(i,h->GetBinContent(11));
      hMCWidth_1->SetBinContent(i,h->GetBinContent(12));
   }

   

   TCanvas *c1 = new TCanvas("c1","",800,800);
   //style();
   
   hDataMass->SetLineWidth(2);
   hMCMass->SetLineWidth(2);

   hMCMass_1->SetMarkerStyle(24);
   hDataMass_1->SetMarkerStyle(25);

   hMCMass_1->SetMarkerColor(kRed);
   hDataMass_1->SetMarkerColor(kBlack);

   hDataMass->SetXTitle("Z p_{T} (GeV)");
   hDataMass->SetYTitle("M_{#mu#mu} (GeV)");
   hDataMass->Draw("e");
   hDataMass->SetMinimum(88);
   hDataMass->SetMaximum(94);
   hMCMass->SetLineColor(kRed);
   hMCMass->SetMarkerColor(kRed);
   hDataMass->SetLineColor(kBlack);
   hDataMass->SetMarkerColor(kBlack);
   hMCMass->Draw("e same");

   hMCMass_1->Draw("p same");
   hDataMass_1->Draw("p same");

   TLegend leg(0.18,0.68,0.58,0.9);
   leg.AddEntry(hMCMass ,"Monte Carlo fitting","lep");
   leg.AddEntry(hDataMass ,Form("%s data fitting",typeofdatatext),"lep");
   leg.AddEntry(hMCMass_1 ,"Monte Carlo mean","p");
   leg.AddEntry(hDataMass_1 ,Form("%s data mean",typeofdatatext),"p");
   leg.SetFillColorAlpha(kWhite,0);
   leg.SetLineColor(kBlack);
   leg.SetLineWidth(1);
   leg.Draw();

   c1->SaveAs(Form("figs/mass/%s/Zmass_%s_loop.png",typeofdata,typeofdata)); 
   c1->Clear();

   ////style();

   hDataWidth->SetLineWidth(2);
   hMCWidth->SetLineWidth(2);

   hMCWidth_1->SetMarkerStyle(24);
   hDataWidth_1->SetMarkerStyle(25);

   hMCWidth_1->SetMarkerColor(kRed);
   hDataWidth_1->SetMarkerColor(kBlack);

   hDataWidth->SetXTitle("Z p_{T} (GeV)");
   hDataWidth->SetYTitle("M_{#mu#mu} width (GeV)");
   
   hDataWidth->SetMinimum(0);
   hDataWidth->SetMaximum(6);
   hMCWidth->SetMinimum(0);
   hMCWidth->SetMaximum(6);

   hDataWidth_1->SetMinimum(0);
   hDataWidth_1->SetMaximum(6);
   hMCWidth_1->SetMinimum(0);
   hMCWidth_1->SetMaximum(6);

   hDataWidth->Draw("e");

   hMCWidth->SetLineColor(kRed);
   hMCWidth->SetMarkerColor(kRed);
   hDataWidth->SetLineColor(kBlack);
   hDataWidth->SetMarkerColor(kBlack);
   hMCWidth->Draw("e same");

   hMCWidth_1->Draw("p same");
   hDataWidth_1->Draw("p same");

   TLegend leg1(0.18,0.68,0.58,0.9);
   leg1.AddEntry(hMCWidth ,"Monte Carlo fitting","lep");
   leg1.AddEntry(hDataWidth ,Form("%s data fitting",typeofdatatext),"lep");
   leg1.AddEntry(hMCWidth_1 ,"Monte Carlo StdDevs","p");
   leg1.AddEntry(hDataWidth_1 ,Form("%s data StdDevs",typeofdatatext),"p");

   leg1.SetFillColorAlpha(kWhite,0);
   leg1.SetLineColor(kBlack);
   leg1.SetLineWidth(1);
   leg1.Draw();

   c1->SaveAs(Form("figs/mass/%s/ZmassWidth_%s_loop.png",typeofdata,typeofdata)); 
   c1->Clear();



/* Seems useless
   delete hDataMass; hDataMass=NULL;
   delete hMCMass; hMCMass=NULL;
   delete hDataWidth; hDataWidth=NULL;
   delete hMCWidth; hMCWidth=NULL;
   delete hDataMass_1; hDataMass_1=NULL;
   delete hMCMass_1; hMCMass_1=NULL;
   delete hDataWidth_1; hDataWidth_1=NULL;
   delete hMCWidth_1; hMCWidth_1=NULL;

   delete gROOT->FindObject("hDataMass");
   delete gROOT->FindObject("hMCMass");
   delete gROOT->FindObject("hDataWidth");
   delete gROOT->FindObject("hMCWidth");
   delete gROOT->FindObject("hDataMass_1");
   delete gROOT->FindObject("hMCMass_1");
   delete gROOT->FindObject("hDataWidth_1");
   delete gROOT->FindObject("hMCWidth_1");
   */

}  


void loopHiBin()
{
   //TH1D *hDataMass = new TH1D("hDataMass","",10,0,200);
   //TH1D *hMCMass = new TH1D("hMCMass","",10,0,200);

   TH1D *hDataMass = new TH1D("hDataMass","",4,cent_diff);
   TH1D *hMCMass = new TH1D("hMCMass","",4,cent_diff);

   TH1D *hDataWidth = new TH1D("hDataWidth","",4,cent_diff);
   TH1D *hMCWidth = new TH1D("hMCWidth","",4,cent_diff);

   TH1D *hDataMass_1 = new TH1D("hDataMass_1","",4,cent_diff);
   TH1D *hMCMass_1 = new TH1D("hMCMass_1","",4,cent_diff);

   TH1D *hDataWidth_1 = new TH1D("hDataWidth_1","",4,cent_diff);
   TH1D *hMCWidth_1 = new TH1D("hMCWidth_1","",4,cent_diff);


   //TH1F *hnew = new TH1F("hnew","rebinned",k,cent_diff);
   
   //for (int i=1;i<=hDataMass->GetNbinsX();i++)
   for (int i=1;i<=4;i++)
   {
      //TH1D *h = ZmassAnalysis_single(0,200,hDataMass->GetBinLowEdge(i),hDataMass->GetBinLowEdge(i+1));
      TH1D *h = ZmassAnalysis_single(0,2000,i-1,i); 
      hDataMass->SetBinContent(i,h->GetBinContent(1));
      hDataMass->SetBinError(i,h->GetBinContent(2));
      hMCMass->SetBinContent(i,h->GetBinContent(3));
      hMCMass->SetBinError(i,h->GetBinContent(4));

      hDataWidth->SetBinContent(i,h->GetBinContent(5));
      hDataWidth->SetBinError(i,h->GetBinContent(6));
      hMCWidth->SetBinContent(i,h->GetBinContent(7));
      hMCWidth->SetBinError(i,h->GetBinContent(8));

      hDataMass_1->SetBinContent(i,h->GetBinContent(9));
      hMCMass_1->SetBinContent(i,h->GetBinContent(10));
      hDataWidth_1->SetBinContent(i,h->GetBinContent(11));
      hMCWidth_1->SetBinContent(i,h->GetBinContent(12));

      //std::cout<<"i = "<<i<<std::endl;
      //std::cout<<"hDataMass = "<<h->GetBinContent(1)<<"+-"<<h->GetBinContent(2)<<std::endl;
      //std::cout<<"hMCMass = "<<h->GetBinContent(3)<<"+-"<<h->GetBinContent(4)<<std::endl;
      //std::cout<<"hDataWidth = "<<h->GetBinContent(5)<<"+-"<<h->GetBinContent(6)<<std::endl;
      //std::cout<<"hMCWidth = "<<h->GetBinContent(7)<<"+-"<<h->GetBinContent(8)<<std::endl;
//
      //std::cout<<"Get: hDataMass = "<<hDataMass->GetBinContent(i)<<"+-"<<hDataMass->GetBinError(i)<<std::endl;
      //std::cout<<"Get: hMCMass = "<<hMCMass->GetBinContent(i)<<"+-"<<hMCMass->GetBinError(i)<<std::endl;
      //std::cout<<"Get: hDataWidth = "<<hDataWidth->GetBinContent(i)<<"+-"<<hDataWidth->GetBinError(i)<<std::endl;
      //std::cout<<"Get: hMCWidth = "<<hMCWidth->GetBinContent(i)<<"+-"<<hMCWidth->GetBinError(i)<<std::endl;
   }
   

   TCanvas *c2 = new TCanvas("c2","",800,800);
   //style();

   hDataMass->SetLineWidth(2);
   hMCMass->SetLineWidth(2);

   hMCMass_1->SetMarkerStyle(24);
   hDataMass_1->SetMarkerStyle(25);

   hMCMass_1->SetMarkerColor(kRed);
   hDataMass_1->SetMarkerColor(kBlack);

   //hDataMass->SetXTitle("Centrality Bin");
   hDataMass->SetXTitle("Centrality (%)");
   hDataMass->SetYTitle("M_{#mu#mu} (GeV)");
   hDataMass->Draw("e");
   hDataMass->SetMinimum(90);
   hDataMass->SetMaximum(92);
   hMCMass->SetLineColor(kRed);
   hMCMass->SetMarkerColor(kRed);
   hDataMass->SetLineColor(kBlack);
   hDataMass->SetMarkerColor(kBlack);
   hMCMass->Draw("e same");

   hMCMass_1->Draw("p same");
   hDataMass_1->Draw("p same");

   TLegend leg(0.18,0.68,0.58,0.9);
   leg.AddEntry(hMCMass ,"Monte Carlo fitting","lep");
   leg.AddEntry(hDataMass ,Form("%s data fitting",typeofdatatext),"lep");
   leg.AddEntry(hMCMass_1 ,"Monte Carlo mean","p");
   leg.AddEntry(hDataMass_1 ,Form("%s data mean",typeofdatatext),"p");
   leg.SetFillColorAlpha(kWhite,0);
   leg.SetLineColor(kBlack);
   leg.SetLineWidth(1);
   leg.Draw();

   c2->SaveAs(Form("figs/mass/%s/Zmass_%s_loopHiBin.png",typeofdata,typeofdata)); 
   c2->Clear();

   ////style();

   hDataWidth->SetLineWidth(2);
   hMCWidth->SetLineWidth(2);

   hMCWidth_1->SetMarkerStyle(24);
   hDataWidth_1->SetMarkerStyle(25);

   hMCWidth_1->SetMarkerColor(kRed);
   hDataWidth_1->SetMarkerColor(kBlack);

   //hDataMass->SetXTitle("Centrality Bin");
   hDataWidth->SetXTitle("Centrality (%)");
   hDataWidth->SetYTitle("M_{#mu#mu} width (GeV)");
   hDataWidth->SetMinimum(0);
   hDataWidth->SetMaximum(8);
   hMCWidth->SetMinimum(0);
   hMCWidth->SetMaximum(8);
   hDataWidth_1->SetMinimum(0);
   hDataWidth_1->SetMaximum(8);
   hMCWidth_1->SetMinimum(0);
   hMCWidth_1->SetMaximum(8);

   hDataWidth->Draw("e");
   hMCWidth->SetLineColor(kRed);
   hMCWidth->SetMarkerColor(kRed);
   hDataWidth->SetLineColor(kBlack);
   hDataWidth->SetMarkerColor(kBlack);
   hMCWidth->Draw("e same");

   hMCWidth_1->Draw("p same");
   hDataWidth_1->Draw("p same");

   TLegend leg1(0.18,0.68,0.58,0.9);
   leg1.AddEntry(hMCWidth ,"Monte Carlo fitting","lep");
   leg1.AddEntry(hDataWidth ,Form("%s data fitting",typeofdatatext),"lep");
   leg1.AddEntry(hMCWidth_1 ,"Monte Carlo StdDevs","p");
   leg1.AddEntry(hDataWidth_1 ,Form("%s data StdDevs",typeofdatatext),"p");
   leg1.SetFillColorAlpha(kWhite,0);
   leg1.SetLineColor(kBlack);
   leg1.SetLineWidth(1);
   leg1.Draw();

   c2->SaveAs(Form("figs/mass/%s/ZmassWidth_%s_loopHiBin.png",typeofdata,typeofdata)); 
   c2->Clear();
}  


void ZmassAnalysis(){

   style();

   ZmassAnalysis_single();
   loop();
   loopHiBin();
   //ZmassAnalysis_single(0,20);

}
