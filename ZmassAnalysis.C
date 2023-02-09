#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <cmath>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <RooPlot.h>
#include <TMath.h>
#include <TF1.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TChain.h>
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
   //TH1D *hDataSame = new TH1D("hDataSame","",40,81.2,101.2);
   TH1D *hMC = new TH1D("hMC","",binsize,81.2,101.2);
   //TH1D *hMCSame = new TH1D("hMCSame","",40,81.2,101.2);

   TCanvas *c = new TCanvas("c","",800,800);

   //style();

   // TTree *tData = (TTree*)infData->Get("t");
   //TTree *tMC = (TTree*)infMC->Get("t");

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

   
   TH1D *hData_pt = new TH1D("hData_pt","",binsize,0,200);
   TH1D *hMC_pt = new TH1D("hMC_pt","",binsize,0,200);


   TH1D *hData_muPt1 = new TH1D("hData_muPt1","",40,0,120);
   TH1D *hMC_muPt1 = new TH1D("hMC_muPt1","",40,0,120);

   TH1D *hData_muPt2 = new TH1D("hData_muPt2","",40,0,120);
   TH1D *hMC_muPt2 = new TH1D("hMC_muPt2","",40,0,120);

   TH1D *hMC_genMuPt1 = new TH1D("hMC_genMuPt1","",40,0,120);
   TH1D *hMC_genMuPt2 = new TH1D("hMC_genMuPt2","",40,0,120);
   
   //tData->Draw("zMass>>hData",Form("zPt>%f&&zPt<%f&&hiBin>=%d&&hiBin<%d",ptL,ptH,centL,centH));
   //tData->Draw("zMass>>hDataSame",Form("zPt>%f&&zPt<%f&&hiBin>=%d&&hiBin<%d",ptL,ptH,centL,centH));

   tData->Draw("zMass>>hData",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   //tData->Draw("zMass>>hDataSame",Form("zPt>%f&&zPt<%f&&hiHF<=%f&&hiHF>%f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   hData->Sumw2();
   //hDataSame->Sumw2();

   //hDataSame->SetLineColor(2);
   //hDataSame->SetMarkerColor(2);


   //tMC->Draw("zMass>>hMC",Form("zPt>%f&&zPt<%f&&hiBin>=%d&&hiBin<%d",ptL,ptH,centL,centH));
   //tMC->Draw("zMass>>hMCSame",Form("zPt>%f&&zPt<%f&&hiBin>=%d&&hiBin<%d",ptL,ptH,centL,centH));

   tMC->Draw("zMass>>hMC",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   //tMC->Draw("zMass>>hMCSame",Form("zPt>%f&&zPt<%f&&hiHF<=%f&&hiHF>%f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   hMC->Sumw2();
   //hMCSame->Sumw2();

   tData->Draw("zEta>>hData_eta",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   tMC->Draw("zEta>>hMC_eta",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   tData->Draw("zPhi>>hData_phi",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   tMC->Draw("zPhi>>hMC_phi",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   tData->Draw("zPt>>hData_pt",Form("hiHF<=%.4f&&hiHF>%.4f",hf_diff[centL],hf_diff[centH]));
   tMC->Draw("zPt>>hMC_pt",Form("hiHF<=%.4f&&hiHF>%.4f",hf_diff[centL],hf_diff[centH]));

   tData->Draw("muPt1>>hData_muPt1",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   tMC->Draw("muPt1>>hMC_muPt1",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   tData->Draw("muPt2>>hData_muPt2",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   tMC->Draw("muPt2>>hMC_muPt2",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   tMC->Draw("genMuPt1>>hMC_genMuPt1",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   tMC->Draw("genMuPt2>>hMC_genMuPt2",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   //int countD = tData->GetEntries(Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
   //std::cout<<"Data = "<<countD<<std::endl;
   //int countM = tMC->GetEntries(Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));
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

   hData_muPt1->Sumw2();
   hMC_muPt1->Sumw2();
   hData_muPt2->Sumw2();
   hMC_muPt2->Sumw2();

   hMC_genMuPt1->Sumw2();
   hMC_genMuPt2->Sumw2();

   hMC_eta->SetMarkerStyle(24);
   hMC_phi->SetMarkerStyle(24);
   hMC_pt->SetMarkerStyle(24);
   hMC_muPt1->SetMarkerStyle(24);
   hMC_muPt2->SetMarkerStyle(24);

   hMC_genMuPt1->SetMarkerStyle(25);
   hMC_genMuPt2->SetMarkerStyle(25);
   
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

   hData_muPt1->SetMarkerColor(kBlack);
   //hMC_muPt1->SetMarkerColor(kRed);
   hData_muPt2->SetMarkerColor(kBlack);
   //hMC_muPt2->SetMarkerColor(kRed);

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

   hData_muPt1->SetLineWidth(2);
   hData_muPt2->SetLineWidth(2);

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
 
   hMC->Fit("f2","LL m");
   hMC->Fit("f2","");
   hMC->Fit("f2","LL");
   hMC->Fit("f2","");
   hMC->Fit("f2","LL m");
   hMC->Fit("f2","");
   hMC->Fit("f2","");
   hMC->Fit("f2","");
   hMC->Fit("f2","");

/*
   std::cout<<"UnbinnedFit 1"<<std::endl;
   tData->UnbinnedFit("f", "zMass", Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]),"");
   //tMC->UnbinnedFit("f2", "zMass", Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]),"");

   std::cout<<"UnbinnedFit 2"<<std::endl;
   tData->UnbinnedFit("f", "zMass", Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]),"E");
   //tMC->UnbinnedFit("f2", "zMass", Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]),"E");

   std::cout<<"UnbinnedFit 3"<<std::endl;
   tData->UnbinnedFit("f", "zMass", Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]),"E M");
   //tMC->UnbinnedFit("f2", "zMass", Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]),"E M");
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
   legMuPt.Draw();

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

   TLatex *ptp = new TLatex(0.18,0.84,Form("%.0f %%< Centrality < %.0f %%",cent_diff[centL],cent_diff[centH]));
   ptp->SetTextFont(42);
   ptp->SetTextSize(0.025);
   ptp->SetNDC(kTRUE);
   ptp->Draw();
   TLatex *ptp2 = new TLatex(0.18,0.89,Form("%.1f < Z p_{T} < %.1f",ptL,ptH));
   ptp2->SetTextFont(42);
   ptp2->SetTextSize(0.025);
   ptp2->SetNDC(kTRUE);
   ptp2->Draw();

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

   hData_muPt1->SetXTitle("#mu p_{T} (GeV)");
   hMC_muPt1->SetXTitle("#mu p_{T} (GeV)");
   hMC_genMuPt1->SetXTitle("#mu p_{T} (GeV)");

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
   
   hData_muPt2->SetXTitle("#mu p_{T} (GeV)");
   hMC_muPt2->SetXTitle("#mu p_{T} (GeV)");
   hMC_genMuPt2->SetXTitle("#mu p_{T} (GeV)");

   legMuPt.Draw();
   ptp->Draw();
   ptp2->Draw();
   hMC_muPt2->SetMinimum(0);
   hData_muPt2->SetMinimum(0);
   hMC_genMuPt2->SetMinimum(0);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_muPt2.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
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
