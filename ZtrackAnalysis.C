#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <cmath>
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
const char *typeofdatatext = "double muon";

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

void ZtrackAnalysis_single(double ptL=0,double ptH=2000,int centL=0,int centH=4,double TptL=0,double TptH=10000)
{

   TCanvas *c = new TCanvas("c","",800,800);

   int binnum = 40;

   TH1D *hData_eta0 = new TH1D("hData_eta0","",binnum,0,6.3);//without leptons
   TH1D *hMC_eta0 = new TH1D("hMC_eta0","",binnum,0,6.3);

   TH1D *hData_phi0 = new TH1D("hData_phi0","",binnum,-3.1415926,3.1415926);
   TH1D *hMC_phi0 = new TH1D("hMC_phi0","",binnum,-3.1415926,3.1415926);


   TH1D *hData_phi0s_0 = new TH1D("hData_phi0s_0","",binnum,-3.1415926,3.1415926);
   TH1D *hMC_phi0s_0 = new TH1D("hMC_phi0s_0","",binnum,-3.1415926,3.1415926);
   TH1D *hData_phi0s_1 = new TH1D("hData_phi0s_1","",binnum,-3.1415926,3.1415926);
   TH1D *hMC_phi0s_1 = new TH1D("hMC_phi0s_1","",binnum,-3.1415926,3.1415926);
   TH1D *hData_phi0s_2 = new TH1D("hData_phi0s_2","",binnum,-3.1415926,3.1415926);
   TH1D *hMC_phi0s_2 = new TH1D("hMC_phi0s_2","",binnum,-3.1415926,3.1415926);
   TH1D *hData_phi0s_3 = new TH1D("hData_phi0s_3","",binnum,-3.1415926,3.1415926);
   TH1D *hMC_phi0s_3 = new TH1D("hMC_phi0s_3","",binnum,-3.1415926,3.1415926);
   TH1D *hData_phi0s_4 = new TH1D("hData_phi0s_4","",binnum,-3.1415926,3.1415926);
   TH1D *hMC_phi0s_4 = new TH1D("hMC_phi0s_4","",binnum,-3.1415926,3.1415926);
   TH1D *hData_phi0s_5 = new TH1D("hData_phi0s_5","",binnum,-3.1415926,3.1415926);
   TH1D *hMC_phi0s_5 = new TH1D("hMC_phi0s_5","",binnum,-3.1415926,3.1415926);

   TH1D *hData_eta0s_0 = new TH1D("hData_eta0s_0","",binnum,0,6.3);
   TH1D *hMC_eta0s_0 = new TH1D("hMC_eta0s_0","",binnum,0,6.3);
   TH1D *hData_eta0s_1 = new TH1D("hData_eta0s_1","",binnum,0,6.3);
   TH1D *hMC_eta0s_1 = new TH1D("hMC_eta0s_1","",binnum,0,6.3);
   TH1D *hData_eta0s_2 = new TH1D("hData_eta0s_2","",binnum,0,6.3);
   TH1D *hMC_eta0s_2 = new TH1D("hMC_eta0s_2","",binnum,0,6.3);
   TH1D *hData_eta0s_3 = new TH1D("hData_eta0s_3","",binnum,0,6.3);
   TH1D *hMC_eta0s_3 = new TH1D("hMC_eta0s_3","",binnum,0,6.3);

   TH2D *hData_etaphi0 = new TH2D("hData_etaphi0","",250,0,6.3,250,-3.1415926,3.1415926);
   TH2D *hMC_etaphi0 = new TH2D("hMC_etaphi0","",250,0,6.3,250,-3.1415926,3.1415926);


   TChain *tMC = new TChain("t");
   tMC->Add("/eos/cms/store/group/phys_heavyions_ops/pchou/MC/*.root?#t");
 
   TChain *tData = new TChain("t");
   tData->Add(Form("/eos/cms/store/group/phys_heavyions_ops/pchou/%s/*.root?#t",typeofdata));

   TCut sel_reg = Form("zPt[0]>%f&&zPt[0]<%f&&hiHF<=%.4f&&hiHF>%.4f&&trackPt>%f&&trackPt<%f",ptL,ptH,hf_diff[centL],hf_diff[centH],TptL,TptH);


   TCut nolepton = "abs(trackPDFId)!=13&&abs(trackPDFId)!=11";

   tData->Draw("trackDeta>>hData_eta0",sel_reg&&nolepton);
   tMC->Draw("trackDeta>>hMC_eta0",sel_reg&&nolepton);
   tData->Draw("trackDphi>>hData_phi0",sel_reg&&nolepton);
   tMC->Draw("trackDphi>>hMC_phi0",sel_reg&&nolepton);
   tMC->Draw("trackDphi:trackDeta>>hMC_etaphi0",sel_reg&&nolepton);
   tData->Draw("trackDphi:trackDeta>>hData_etaphi0",sel_reg&&nolepton);


   tData->Draw("trackDphi>>hData_phi0s_0",sel_reg&&nolepton&&"trackDeta>=0&&trackDeta<3.1415926/2");
   tMC->Draw("trackDphi>>hMC_phi0s_0",sel_reg&&nolepton&&"trackDeta>=0&&trackDeta<3.1415926/2");

   tData->Draw("trackDphi>>hData_phi0s_1",sel_reg&&nolepton&&"trackDeta>=0&&trackDeta<3.1415926/8");
   tMC->Draw("trackDphi>>hMC_phi0s_1",sel_reg&&nolepton&&"trackDeta>=0&&trackDeta<3.1415926/8");
   tData->Draw("trackDphi>>hData_phi0s_2",sel_reg&&nolepton&&"trackDeta>=3.1415926/8&&trackDeta<3.1415926/4");
   tMC->Draw("trackDphi>>hMC_phi0s_2",sel_reg&&nolepton&&"trackDeta>=3.1415926/8&&trackDeta<3.1415926/4");
   tData->Draw("trackDphi>>hData_phi0s_3",sel_reg&&nolepton&&"trackDeta>=3.1415926/4&&trackDeta<3.1415926/2");
   tMC->Draw("trackDphi>>hMC_phi0s_3",sel_reg&&nolepton&&"trackDeta>=3.1415926/4&&trackDeta<3.1415926/2");
   tData->Draw("trackDphi>>hData_phi0s_4",sel_reg&&nolepton&&"trackDeta>=3.1415926/2&&trackDeta<3.1415926");
   tMC->Draw("trackDphi>>hMC_phi0s_4",sel_reg&&nolepton&&"trackDeta>=3.1415926/2&&trackDeta<3.1415926");
   tData->Draw("trackDphi>>hData_phi0s_5",sel_reg&&nolepton&&"trackDeta>=3.1415926&&trackDeta<3.1415926*2");
   tMC->Draw("trackDphi>>hMC_phi0s_5",sel_reg&&nolepton&&"trackDeta>=3.1415926&&trackDeta<3.1415926*2");

   tData->Draw("trackDeta>>hData_eta0s_0",sel_reg&&nolepton&&"abs(trackDphi)>=0&&abs(trackDphi)<3.1415926/4");
   tMC->Draw("trackDeta>>hMC_eta0s_0",sel_reg&&nolepton&&"abs(trackDphi)>=0&&abs(trackDphi)<3.1415926/4");

   tData->Draw("trackDeta>>hData_eta0s_1",sel_reg&&nolepton&&"abs(trackDphi)>=0&&abs(trackDphi)<3.1415926/4");
   tMC->Draw("trackDeta>>hMC_eta0s_1",sel_reg&&nolepton&&"abs(trackDphi)>=0&&abs(trackDphi)<3.1415926/4");
   tData->Draw("trackDeta>>hData_eta0s_2",sel_reg&&nolepton&&"abs(trackDphi)>=3.1415926/4&&abs(trackDphi)<3.1415926/2");
   tMC->Draw("trackDeta>>hMC_eta0s_2",sel_reg&&nolepton&&"abs(trackDphi)>=3.1415926/4&&abs(trackDphi)<3.1415926/2");
   tData->Draw("trackDeta>>hData_eta0s_3",sel_reg&&nolepton&&"abs(trackDphi)>=3.1415926/2&&abs(trackDphi)<3.1415926");
   tMC->Draw("trackDeta>>hMC_eta0s_3",sel_reg&&nolepton&&"abs(trackDphi)>=3.1415926/2&&abs(trackDphi)<3.1415926");

   int countD0 = hData_eta0->GetEntries();
   std::cout<<"Data 0 = "<<countD0<<std::endl;
   int countM0 = hMC_eta0->GetEntries();
   std::cout<<"MC 0 = "<<countM0<<std::endl;

   int countDs = hData_phi0s_0->GetEntries();
   std::cout<<"Data 1 = "<<countDs<<std::endl;
   int countMs = hMC_phi0s_0->GetEntries();
   std::cout<<"MC 1 = "<<countMs<<std::endl;

   int countDe = hData_eta0s_0->GetEntries();
   std::cout<<"Data 1 = "<<countDe<<std::endl;
   int countMe = hMC_eta0s_0->GetEntries();
   std::cout<<"MC 1 = "<<countMe<<std::endl;

   hData_eta0->Sumw2();
   hMC_eta0->Sumw2();
   hData_phi0->Sumw2();
   hMC_phi0->Sumw2();

   hMC_phi0s_0->Sumw2();
   hMC_phi0s_1->Sumw2();
   hMC_phi0s_2->Sumw2();
   hMC_phi0s_3->Sumw2();
   hMC_phi0s_4->Sumw2();
   hMC_phi0s_5->Sumw2();
   hData_phi0s_0->Sumw2();
   hData_phi0s_1->Sumw2();
   hData_phi0s_2->Sumw2();
   hData_phi0s_3->Sumw2();
   hData_phi0s_4->Sumw2();
   hData_phi0s_5->Sumw2();

   hMC_eta0s_0->Sumw2();
   hMC_eta0s_1->Sumw2();
   hMC_eta0s_2->Sumw2();
   hMC_eta0s_3->Sumw2();
   hData_eta0s_0->Sumw2();
   hData_eta0s_1->Sumw2();
   hData_eta0s_2->Sumw2();
   hData_eta0s_3->Sumw2();

   hMC_eta0->SetMarkerStyle(24);
   hMC_phi0->SetMarkerStyle(24);

   hMC_phi0s_0->SetMarkerStyle(24);
   hMC_phi0s_1->SetMarkerStyle(24);
   hMC_phi0s_2->SetMarkerStyle(25);
   hMC_phi0s_3->SetMarkerStyle(27);
   hMC_phi0s_4->SetMarkerStyle(28);
   hMC_phi0s_5->SetMarkerStyle(42);
   hData_phi0s_0->SetMarkerStyle(24);
   hData_phi0s_1->SetMarkerStyle(24);
   hData_phi0s_2->SetMarkerStyle(25);
   hData_phi0s_3->SetMarkerStyle(27);
   hData_phi0s_4->SetMarkerStyle(28);
   hData_phi0s_5->SetMarkerStyle(42);

   hMC_eta0s_0->SetMarkerStyle(24);
   hMC_eta0s_1->SetMarkerStyle(24);
   hMC_eta0s_2->SetMarkerStyle(25);
   hMC_eta0s_3->SetMarkerStyle(27);
   hData_eta0s_0->SetMarkerStyle(24);
   hData_eta0s_1->SetMarkerStyle(24);
   hData_eta0s_2->SetMarkerStyle(25);
   hData_eta0s_3->SetMarkerStyle(27);

   hData_eta0->SetMarkerColor(kBlack);
   hMC_eta0->SetMarkerColor(kRed);
   hData_phi0->SetMarkerColor(kBlack);
   hMC_phi0->SetMarkerColor(kRed);

   hData_eta0->SetLineColor(kBlack);
   hMC_eta0->SetLineColor(kRed);
   hData_phi0->SetLineColor(kBlack);
   hMC_phi0->SetLineColor(kRed);

   hMC_phi0s_0->SetMarkerColor(kRed);
   hMC_phi0s_1->SetMarkerColor(TColor::GetColor("#FF1F5B"));
   hMC_phi0s_2->SetMarkerColor(TColor::GetColor("#00CD6C"));
   hMC_phi0s_3->SetMarkerColor(TColor::GetColor("#009ADE"));
   hMC_phi0s_4->SetMarkerColor(TColor::GetColor("#AF58BA"));
   hMC_phi0s_5->SetMarkerColor(TColor::GetColor("#F28522"));
   hData_phi0s_0->SetMarkerColor(kBlack);
   hData_phi0s_1->SetMarkerColor(TColor::GetColor("#FF1F5B"));
   hData_phi0s_2->SetMarkerColor(TColor::GetColor("#00CD6C"));
   hData_phi0s_3->SetMarkerColor(TColor::GetColor("#009ADE"));
   hData_phi0s_4->SetMarkerColor(TColor::GetColor("#AF58BA"));
   hData_phi0s_5->SetMarkerColor(TColor::GetColor("#F28522"));

   hMC_eta0s_0->SetMarkerColor(kRed);
   hMC_eta0s_1->SetMarkerColor(TColor::GetColor("#FF1F5B"));
   hMC_eta0s_2->SetMarkerColor(TColor::GetColor("#00CD6C"));
   hMC_eta0s_3->SetMarkerColor(TColor::GetColor("#009ADE"));
   hData_eta0s_0->SetMarkerColor(kBlack);
   hData_eta0s_1->SetMarkerColor(TColor::GetColor("#FF1F5B"));
   hData_eta0s_2->SetMarkerColor(TColor::GetColor("#00CD6C"));
   hData_eta0s_3->SetMarkerColor(TColor::GetColor("#009ADE"));

   hMC_phi0s_0->SetLineColor(kRed);
   hData_phi0s_1->SetLineColor(kBlack);

   hMC_phi0s_1->SetLineColor(TColor::GetColor("#FF1F5B"));
   hMC_phi0s_2->SetLineColor(TColor::GetColor("#00CD6C"));
   hMC_phi0s_3->SetLineColor(TColor::GetColor("#009ADE"));
   hMC_phi0s_4->SetLineColor(TColor::GetColor("#AF58BA"));
   hMC_phi0s_5->SetLineColor(TColor::GetColor("#F28522"));
   hData_phi0s_1->SetLineColor(TColor::GetColor("#FF1F5B"));
   hData_phi0s_2->SetLineColor(TColor::GetColor("#00CD6C"));
   hData_phi0s_3->SetLineColor(TColor::GetColor("#009ADE"));
   hData_phi0s_4->SetLineColor(TColor::GetColor("#AF58BA"));
   hData_phi0s_5->SetLineColor(TColor::GetColor("#F28522"));


   hMC_eta0s_0->SetLineColor(kRed);
   hData_eta0s_1->SetLineColor(kBlack);

   hMC_eta0s_1->SetLineColor(TColor::GetColor("#FF1F5B"));
   hMC_eta0s_2->SetLineColor(TColor::GetColor("#00CD6C"));
   hMC_eta0s_3->SetLineColor(TColor::GetColor("#009ADE"));
   hData_eta0s_1->SetLineColor(TColor::GetColor("#FF1F5B"));
   hData_eta0s_2->SetLineColor(TColor::GetColor("#00CD6C"));
   hData_eta0s_3->SetLineColor(TColor::GetColor("#009ADE"));

   hData_eta0->Scale(1./hData_eta0->GetEntries());
   hData_phi0->Scale(1./hData_phi0->GetEntries());
   hMC_eta0->Scale(1./hMC_eta0->GetEntries());
   hMC_phi0->Scale(1./hMC_phi0->GetEntries());

   hData_etaphi0->Scale(1./hData_etaphi0->GetEntries());
   hMC_etaphi0->Scale(1./hMC_etaphi0->GetEntries());

   hMC_phi0s_0->Scale(1./hMC_phi0s_0->GetEntries());
   hData_phi0s_0->Scale(1./hData_phi0s_0->GetEntries());

   hMC_phi0s_1->Scale(1./hMC_phi0s_1->GetEntries());
   hMC_phi0s_2->Scale(1./hMC_phi0s_2->GetEntries());
   hMC_phi0s_3->Scale(1./hMC_phi0s_3->GetEntries());
   hMC_phi0s_4->Scale(1./hMC_phi0s_4->GetEntries());
   hMC_phi0s_5->Scale(1./hMC_phi0s_5->GetEntries());
   hData_phi0s_1->Scale(1./hData_phi0s_1->GetEntries());
   hData_phi0s_2->Scale(1./hData_phi0s_2->GetEntries());
   hData_phi0s_3->Scale(1./hData_phi0s_3->GetEntries());
   hData_phi0s_4->Scale(1./hData_phi0s_4->GetEntries());
   hData_phi0s_5->Scale(1./hData_phi0s_5->GetEntries());

   hMC_eta0s_0->Scale(1./hMC_eta0s_0->GetEntries());
   hData_eta0s_0->Scale(1./hData_eta0s_0->GetEntries());

   hMC_eta0s_1->Scale(1./hMC_eta0s_1->GetEntries());
   hMC_eta0s_2->Scale(1./hMC_eta0s_2->GetEntries());
   hMC_eta0s_3->Scale(1./hMC_eta0s_3->GetEntries());
   hData_eta0s_1->Scale(1./hData_eta0s_1->GetEntries());
   hData_eta0s_2->Scale(1./hData_eta0s_2->GetEntries());
   hData_eta0s_3->Scale(1./hData_eta0s_3->GetEntries());

   hMC_phi0s_0->SetMarkerSize(0.5);
   hData_phi0s_0->SetMarkerSize(0.5);

   hMC_phi0s_1->SetMarkerSize(0.5);
   hMC_phi0s_2->SetMarkerSize(0.5);
   hMC_phi0s_3->SetMarkerSize(0.5);
   hMC_phi0s_4->SetMarkerSize(0.5);
   hMC_phi0s_5->SetMarkerSize(0.5);
   hData_phi0s_1->SetMarkerSize(0.5);
   hData_phi0s_2->SetMarkerSize(0.5);
   hData_phi0s_3->SetMarkerSize(0.5);
   hData_phi0s_4->SetMarkerSize(0.5);
   hData_phi0s_5->SetMarkerSize(0.5);

   hMC_eta0s_0->SetMarkerSize(0.5);
   hData_eta0s_0->SetMarkerSize(0.5);

   hMC_eta0s_1->SetMarkerSize(0.5);
   hMC_eta0s_2->SetMarkerSize(0.5);
   hMC_eta0s_3->SetMarkerSize(0.5);
   hData_eta0s_1->SetMarkerSize(0.5);
   hData_eta0s_2->SetMarkerSize(0.5);
   hData_eta0s_3->SetMarkerSize(0.5);

   TLegend leg(0.58,0.78,0.98,0.9);
   leg.AddEntry(hMC_eta0 ,"Monte Carlo: DYLL","lep");
   leg.AddEntry(hData_eta0 ,Form("Data: %s",typeofdatatext),"lep");
   leg.SetFillColorAlpha(kWhite,0);
   leg.SetLineColor(kBlack);
   leg.SetLineWidth(1);

   TLatex *pt = new TLatex(0.18,0.88,Form("%.0f %%< Centrality < %.0f %%",cent_diff[centL],cent_diff[centH]));
   pt->SetTextFont(42);
   pt->SetTextSize(0.03);
   pt->SetNDC(kTRUE);

   TLatex *pt2 = new TLatex(0.18,0.82,Form("%.1f < Z p_{T} < %.1f",ptL,ptH));
   pt2->SetTextFont(42);
   pt2->SetTextSize(0.03);
   pt2->SetNDC(kTRUE);

   TLatex *pt3 = new TLatex(0.18,0.76,Form("%.1f < Track p_{T} < %.1f",TptL,TptH));
   pt3->SetTextFont(42);
   pt3->SetTextSize(0.03);
   pt3->SetNDC(kTRUE);

   TLatex *ptN0 = new TLatex(0.5,0.97,Form("N_{MC} = %d, N_{Data} = %d",countM0,countD0));
   ptN0->SetTextFont(42);
   ptN0->SetTextSize(0.03);
   ptN0->SetNDC(kTRUE);

   TLatex *ptNs = new TLatex(0.5,0.97,Form("N_{MC} = %d, N_{Data} = %d",countMs,countDs));
   ptNs->SetTextFont(42);
   ptNs->SetTextSize(0.03);
   ptNs->SetNDC(kTRUE);

   TLatex *ptNe = new TLatex(0.5,0.97,Form("N_{MC} = %d, N_{Data} = %d",countMe,countDe));
   ptNe->SetTextFont(42);
   ptNe->SetTextSize(0.03);
   ptNe->SetNDC(kTRUE);

   TLatex *ptNs1 = new TLatex(0.15,0.97,"|#Delta#eta| < #pi/2");
   ptNs1->SetTextFont(42);
   ptNs1->SetTextSize(0.03);
   ptNs1->SetNDC(kTRUE);

   TLatex *ptNe1 = new TLatex(0.15,0.97,"|#Delta#phi| < #pi/4");
   ptNe1->SetTextFont(42);
   ptNe1->SetTextSize(0.03);
   ptNe1->SetNDC(kTRUE);

   // == Start drawing == //

   gSystem->Exec(Form("mkdir -p figs/track/%s",typeofdata));

   // Draw eta 

   double max1 = hMC_eta0->GetMaximum();
   double max2 = hData_eta0->GetMaximum();
   
   if(max1<max2) hData_eta0->Draw();
   else hMC_eta0->Draw();
   hMC_eta0->Draw("same");
   hData_eta0->Draw("same");

   hMC_eta0->SetXTitle("|#Delta#eta_{Z,track}|");
   hData_eta0->SetXTitle("|#Delta#eta_{Z,track}|");

   leg.Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();
   hMC_eta0->SetMinimum(0);
   hMC_eta0->SetMaximum(3.0/binnum); 
   hData_eta0->SetMinimum(0);
   hData_eta0->SetMaximum(3.0/binnum); 

   ptN0->Draw();

   c->SaveAs(Form("figs/track/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_Deta0.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH],TptL,TptH)); 
   c->Clear();

   max1 = hMC_phi0->GetMaximum();
   max2 = hData_phi0->GetMaximum();
   
   if(max1<max2) hData_phi0->Draw();
   else hMC_phi0->Draw();
   hMC_phi0->Draw("same");
   hData_phi0->Draw("same");

   if(max1<max2) max1=max2;

   hMC_phi0->SetXTitle("#Delta#phi_{Z,track}");
   hData_phi0->SetXTitle("#Delta#phi_{Z,track}");

   leg.Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();
   hMC_phi0->SetMinimum(0);
   hMC_phi0->SetMaximum(1.6*max1);
   hData_phi0->SetMinimum(0);
   hData_phi0->SetMaximum(1.6*max1);

   ptN0->Draw();

   c->SaveAs(Form("figs/track/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_Dphi0.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH],TptL,TptH)); 
   c->Clear();

   c->SetCanvasSize(1400,800);
   c->Divide(2);
   c->cd(1);

   hMC_etaphi0->Draw("COLZ");
   hMC_etaphi0->GetYaxis()->SetTitle("MC #Delta#phi_{Z,track}");
   hMC_etaphi0->GetXaxis()->SetTitle("MC |#Delta#eta_{Z,track}|");
   hMC_etaphi0->GetXaxis()->SetNdivisions(50205,kFALSE);

   pt->Draw();
   pt2->Draw();
   pt3->Draw();

   c->cd(2);
   hData_etaphi0->Draw("COLZ");
   hData_etaphi0->GetYaxis()->SetTitle("Data #Delta#phi_{Z,track}");
   hData_etaphi0->GetXaxis()->SetTitle("Data |#Delta#eta_{Z,track}|");
   hData_etaphi0->GetXaxis()->SetNdivisions(50205,kFALSE);

   ptN0->Draw();

   c->SaveAs(Form("figs/track/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_Detaphi0.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH],TptL,TptH)); 
   //c->SetCanvasSize(800,800);
   //c->Divide(1);
   c->Clear();

   TLegend leg_s(0.58,0.7,0.98,0.9);
   leg_s.AddEntry(hMC_phi0s_1 ,"0 < |#Delta#eta| < #pi/8","lep");
   leg_s.AddEntry(hMC_phi0s_2 ,"#pi/8 < |#Delta#eta| < #pi/4","lep");
   leg_s.AddEntry(hMC_phi0s_3 ,"#pi/4 < |#Delta#eta| < #pi/2","lep");
   leg_s.AddEntry(hMC_phi0s_4 ,"#pi/2 < |#Delta#eta| < #pi","lep");
   leg_s.AddEntry(hMC_phi0s_5 ,"#pi < |#Delta#eta| < 2#pi","lep");
   leg_s.SetFillColorAlpha(kWhite,0);
   leg_s.SetLineColor(kBlack);
   leg_s.SetLineWidth(1);

   //c->SetCanvasSize(1400,800);
   c->Divide(2);
   c->cd(1);

   hMC_phi0s_1->Draw();
   hMC_phi0s_2->Draw("same");
   hMC_phi0s_3->Draw("same");
   hMC_phi0s_4->Draw("same");
   hMC_phi0s_5->Draw("same");

   hMC_phi0s_1->SetXTitle("MC #Delta#phi_{Z,track}");

   hMC_phi0s_1->SetMinimum(0);
   max1 = hMC_phi0s_1->GetMaximum();
   hMC_phi0s_1->SetMaximum(1.6*max1);

   leg_s.Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();

   c->cd(2);
   hData_phi0s_1->Draw();
   hData_phi0s_2->Draw("same");
   hData_phi0s_3->Draw("same");
   hData_phi0s_4->Draw("same");
   hData_phi0s_5->Draw("same");

   hData_phi0s_1->SetXTitle("Data #Delta#phi_{Z,track}");
   
   hData_phi0s_1->SetMinimum(0);
   //max2 = hData_phi0s_1->GetMaximum();
   hData_phi0s_1->SetMaximum(1.6*max1);

   leg_s.Draw();
   ptN0->Draw();

   c->SaveAs(Form("figs/track/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_Dphi0s.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH],TptL,TptH)); 

   c->Clear();

   TLegend leg_e(0.58,0.7,0.98,0.9);
   leg_e.AddEntry(hMC_eta0s_1 ,"0 < |#Delta#phi| < #pi/4","lep");
   leg_e.AddEntry(hMC_eta0s_2 ,"#pi/4 < |#Delta#phi| < #pi/2","lep");
   leg_e.AddEntry(hMC_eta0s_3 ,"#pi/2 < |#Delta#phi| < #pi","lep");
   leg_e.SetFillColorAlpha(kWhite,0);
   leg_e.SetLineColor(kBlack);
   leg_e.SetLineWidth(1);

   c->Divide(2);
   c->cd(1);

   hMC_eta0s_1->Draw();
   hMC_eta0s_2->Draw("same");
   hMC_eta0s_3->Draw("same");

   hMC_eta0s_1->SetXTitle("MC #Delta#eta_{Z,track}");

   hMC_eta0s_1->SetMinimum(0);
   max2 = hMC_eta0s_1->GetMaximum();
   hMC_eta0s_1->SetMaximum(1.6*max2);

   leg_e.Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();

   c->cd(2);
   hData_eta0s_1->Draw();
   hData_eta0s_2->Draw("same");
   hData_eta0s_3->Draw("same");

   hData_eta0s_1->SetXTitle("Data #Delta#eta_{Z,track}");

   hData_eta0s_1->SetMinimum(0);
   hData_eta0s_1->SetMaximum(1.6*max2);

   leg_e.Draw();

   ptN0->Draw();

   c->SaveAs(Form("figs/track/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_Deta0s.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH],TptL,TptH)); 
  
   c->SetCanvasSize(800,800);
   c->Divide(1);
   c->Clear();

   
   hMC_phi0s_0->SetMarkerColor(kRed);
   hData_phi0s_0->SetMarkerColor(kBlack);
   hMC_phi0s_0->SetLineColor(kRed);
   hData_phi0s_0->SetLineColor(kBlack);

   max1 = hMC_phi0s_0->GetMaximum();
   max2 = hData_phi0s_0->GetMaximum();
   
   if(max1<max2) hData_phi0s_0->Draw();
   else hMC_phi0s_0->Draw();
   hMC_phi0s_0->Draw("same");
   hData_phi0s_0->Draw("same");

   hMC_phi0s_0->SetXTitle("#Delta#phi_{Z,track}");
   hData_phi0s_0->SetXTitle("#Delta#phi_{Z,track}");

   leg.Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();
   hMC_phi0s_0->SetMinimum(0);
   hMC_phi0s_0->SetMaximum(1.6/binnum);
   hData_phi0s_0->SetMinimum(0);
   hData_phi0s_0->SetMaximum(1.6/binnum);


   ptNs->Draw();
   ptNs1->Draw();

   c->SaveAs(Form("figs/track/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_Dphi0s_1d_2.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH],TptL,TptH)); 
   c->Clear();

   hMC_eta0s_0->SetMarkerColor(kRed);
   hData_eta0s_0->SetMarkerColor(kBlack);
   hMC_eta0s_0->SetLineColor(kRed);
   hData_eta0s_0->SetLineColor(kBlack);

   max1 = hMC_phi0s_0->GetMaximum();
   max2 = hData_phi0s_0->GetMaximum();
   
   if(max1<max2) hData_eta0s_0->Draw();
   else hMC_eta0s_0->Draw();
   hMC_eta0s_0->Draw("same");
   hData_eta0s_0->Draw("same");

   hMC_eta0s_0->SetXTitle("#Delta#eta_{Z,track}");
   hData_eta0s_0->SetXTitle("#Delta#eta_{Z,track}");

   leg.Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();
   hMC_eta0s_0->SetMinimum(0);
   hMC_eta0s_0->SetMaximum(3.2/binnum);
   hData_eta0s_0->SetMinimum(0);
   hData_eta0s_0->SetMaximum(3.2/binnum);

   ptNe->Draw();
   ptNe1->Draw();

   c->SaveAs(Form("figs/track/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_Deta0s_1d.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH],TptL,TptH)); 

}


void ZtrackAnalysis(int casenum){

   style();

   switch (casenum){
      case 0: ZtrackAnalysis_single(); break;
      case 1: ZtrackAnalysis_single(30,40);break;
      case 2: ZtrackAnalysis_single(40,60);break;
      case 3: ZtrackAnalysis_single(60,80);break;
      case 4: ZtrackAnalysis_single(80,100); break;
      case 5: ZtrackAnalysis_single(0,200,0,1);break;
      case 6: ZtrackAnalysis_single(0,200,1,2);break;
      case 7: ZtrackAnalysis_single(0,200,2,3);break;
      case 8: ZtrackAnalysis_single(0,200,3,4);break;
      case 9: ZtrackAnalysis_single(0,2000,0,4,0,10);break;
      case 10: ZtrackAnalysis_single(0,2000,0,4,10,20);break;
      case 11: ZtrackAnalysis_single(0,2000,0,4,20,40);break;
      case 12: ZtrackAnalysis_single(0,2000,0,4,40,70);break;
      case 13: ZtrackAnalysis_single(0,2000,0,4,70,100);break;
      case 14: ZtrackAnalysis_single(0,2000,0,4,100,10000);break;
      case 15: ZtrackAnalysis_single(0,2000,3,4,0,10);break;
      case 16: ZtrackAnalysis_single(0,2000,3,4,10,20);break;
      case 17: ZtrackAnalysis_single(0,2000,3,4,20,40);break;
      case 18: ZtrackAnalysis_single(0,2000,3,4,40,70);break;
      case 19: ZtrackAnalysis_single(0,2000,3,4,70,100);break;
      case 20: ZtrackAnalysis_single(0,2000,3,4,100,10000);break;
      case 21: ZtrackAnalysis_single(0,2000,2,3,0,10);break;
      case 22: ZtrackAnalysis_single(0,2000,2,3,10,20);break;
      case 23: ZtrackAnalysis_single(0,2000,2,3,20,40);break;
      case 24: ZtrackAnalysis_single(0,2000,2,3,40,70);break;
      case 25: ZtrackAnalysis_single(0,2000,2,3,70,100);break;
      case 26: ZtrackAnalysis_single(0,2000,2,3,100,10000);break;
      case 27: ZtrackAnalysis_single(0,2000,2,4,0,10);break;
      case 28: ZtrackAnalysis_single(0,2000,2,4,10,20);break;
      case 29: ZtrackAnalysis_single(0,2000,2,4,20,40);break;
      case 30: ZtrackAnalysis_single(0,2000,2,4,40,70);break;
      case 31: ZtrackAnalysis_single(0,2000,2,4,70,100);break;
      case 32: ZtrackAnalysis_single(0,2000,2,4,100,10000);break;
      case 33: ZtrackAnalysis_single(60,100,0,4,0,10);break;
      case 34: ZtrackAnalysis_single(60,100,0,4,10,20);break;
      case 35: ZtrackAnalysis_single(60,100,0,4,20,40);break;
      case 36: ZtrackAnalysis_single(60,100,0,4,40,70);break;
      case 37: ZtrackAnalysis_single(60,100,0,4,70,100);break;
      case 38: ZtrackAnalysis_single(60,100,0,4,100,10000);break;
   }
}
