#ifndef __CINT__
#include "RooGlobalFunc.h"
#else
class Roo2DKeysPdf;
#endif
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooBifurGauss.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/DistFunc.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TChain.h"
#include <cmath>
#ifndef __CINT__
#include "RooCFunction1Binding.h" 
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h" 
using namespace RooFit ;


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

//TH1D* ZmassRoofit_single(double ptL=0,double ptH=2000,int centL=0,int centH=200)
TH1D* ZmassRoofit_single(double ptL=0,double ptH=2000,int centL=0,int centH=4)
{

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
   gErrorIgnoreLevel = kError;

   RooRealVar Zmass("Zmass","M_{#mu#mu}",81.2,101.2,"GeV/c^{2}");
   //RooRealVar zPt("zPt","#zPt",0,200,"GeV/c");
   //RooRealVar hiHF("hiHF","hiHF",0,1000000,"");
   RooDataSet dataMC("dataMC","data MC",RooArgSet(Zmass)) ;
   RooDataSet dataData("dataData","data Data",RooArgSet(Zmass)) ;

   TCanvas *c = new TCanvas("c","",800,800);
   
 
   TChain *tMC = new TChain("t");
   tMC->Add("/eos/cms/store/group/phys_heavyions_ops/pchou/MC/*.root?#t");
   
   TChain *tData = new TChain("t");
   tData->Add(Form("/eos/cms/store/group/phys_heavyions_ops/pchou/%s/*.root?#t",typeofdata));

   TH1D *h = new TH1D("h","",100,0,1);

   double TotalEventsMC = 0;
   double TotalEventsData = 0;

   Float_t tagzMassMC, tagzPtMC, taghiHFMC;
   tMC->SetBranchAddress("zMass",&tagzMassMC);
   tMC->SetBranchAddress("zPt", &tagzPtMC);
   tMC->SetBranchAddress("hiHF",&taghiHFMC);

   Float_t tagzMassData, tagzPtData, taghiHFData;
   tData->SetBranchAddress("zMass",&tagzMassData);
   tData->SetBranchAddress("zPt", &tagzPtData);
   tData->SetBranchAddress("hiHF",&taghiHFData);

   int no_of_eventsMC   = (int)tMC->GetEntries();
   int no_of_eventsData = (int)tData->GetEntries();

   for(int idx=0;idx<no_of_eventsMC;idx++)
   { 
     tMC->GetEntry(idx);  
     if(tagzPtMC>ptL&&tagzPtMC<ptH&&taghiHFMC>hf_diff[centH]&&taghiHFMC<=hf_diff[centL]){
         Zmass = tagzMassMC;
         dataMC.add(RooArgSet(Zmass));
         TotalEventsMC++;
     }
   }

   for(int idx=0;idx<no_of_eventsData;idx++)
   { 
     tData->GetEntry(idx);  
     if(tagzPtData>ptL&&tagzPtData<ptH&&taghiHFData>hf_diff[centH]&&taghiHFData<=hf_diff[centL]){
         Zmass = tagzMassData;
         dataData.add(RooArgSet(Zmass));
         TotalEventsData++;
     }
   }

   Zmass.setRange("full_range",81.2,101.2);





   
   TH1D *hData_pt = new TH1D("hData_pt","",40,0,200);
   TH1D *hMC_pt = new TH1D("hMC_pt","",40,0,200);
   
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

   hMC_eta->SetMarkerStyle(24);
   hMC_phi->SetMarkerStyle(24);
   hMC_pt->SetMarkerStyle(24);
   
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

   hData_eta->SetLineColor(kBlack);
   hMC_eta->SetLineColor(kRed);
   hData_phi->SetLineColor(kBlack);
   hMC_phi->SetLineColor(kRed);
   hData_pt->SetLineColor(kBlack);
   hMC_pt->SetLineColor(kRed);

   //hMC->GetXaxis()->SetTitleSize(48);
   //hMC->GetXaxis()->SetTitleFont(43);

   //hData->Scale(1./hData->GetEntries());
   //hMC->Scale(1./hMC->GetEntries());

   hData->Scale(1./hData->Integral("width"));
   hMC->Scale(1./hMC->Integral("width"));

   hMC->Draw();
   hMC->SetXTitle("M_{#mu#mu} (GeV)");
   hData->Draw("same"); 
   //hMC->SetLineColor(2);
   //hMC->SetMarkerColor(2);

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

   //TLatex *pt = new TLatex(0.18,0.82,Form("%d < Centrality < %d",centL,centH));
   TLatex *pt = new TLatex(0.18,0.82,Form("%.0f %%< Centrality < %.0f %%",cent_diff[centL],cent_diff[centH]));
   pt->SetTextFont(42);
   //std::cout<<"TextSize = "<<pt->GetTextSize()<<std::endl;
   pt->SetTextSize(0.03);
   pt->SetNDC(kTRUE);
   pt->Draw();
   TLatex *pt2 = new TLatex(0.18,0.88,Form("%.1f < p_{T} < %.1f",ptL,ptH));
   pt2->SetTextFont(42);
   //std::cout<<"TextSize2 = "<<pt2->GetTextSize()<<std::endl;
   pt2->SetTextSize(0.03);
   pt2->SetNDC(kTRUE);
   pt2->Draw();

   TLatex *ptN = new TLatex(0.6,0.97,Form("N_{MC} = %d, N_{Data} = %d",countM,countD));
   ptN->SetTextFont(42);
   ptN->SetTextSize(0.03);
   ptN->SetNDC(kTRUE);
   ptN->Draw();

   gSystem->Exec(Form("mkdir -p figs/mass/%s",typeofdata));

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   ////style();
   

   hMC_eta->Draw();
   hMC_eta->SetXTitle("#eta_{Z}");
   hData_eta->Draw("same"); 
   leg.Draw();
   pt->Draw();
   pt2->Draw();
   hMC_eta->SetMinimum(0);
   if (ptL==0&&ptH==200) hMC_eta->SetMaximum(0.1); 
   else if(ptL==0&&ptH==20) hMC_eta->SetMaximum(0.1); 
   else if(ptL==80&&ptH==100) hMC_eta->SetMaximum(0.2); 
   else if(ptL==0&&ptH==2000) hMC_eta->SetMaximum(0.1); 
   else hMC_eta->SetMaximum(0.15);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_eta.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();

   ////style();

   hMC_phi->Draw();
   hMC_phi->SetXTitle("#phi_{Z}");

   hData_phi->Draw("same"); 
   leg.Draw();
   pt->Draw();
   pt2->Draw();
   hMC_phi->SetMinimum(0);
   hMC_phi->SetMaximum(0.08);
   if(ptL==80&&ptH==100) hMC_phi->SetMaximum(0.1);
   else if(ptL==0&&ptH==2000) hMC_phi->SetMaximum(0.05);

   ptN->Draw();

   c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f_phi.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();


   if(ptL==0&&(ptH==200||ptH==2000)){
      ////style();
      hMC_pt->Draw();
      hMC_pt->SetXTitle("Z p_{T} (GeV)");
      hMC_pt->GetXaxis()->SetLimits(0,200);
      hData_pt->Draw("same"); 
      hMC_pt->SetMinimum(0);
      hMC_pt->SetMaximum(0.3);
      leg.Draw();
      pt->Draw();
      pt->SetY(0.88);
      //pt2->Draw();
   
      ptN->Draw();

      c->SaveAs(Form("figs/mass/%s/Zmass_%s_%.0f_%.0f_pt.png",typeofdata,typeofdata,cent_diff[centL],cent_diff[centH])); 
      hMC_pt->SetMinimum(0.00001);
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
      TH1D *h = ZmassRoofit_single(hDataMass->GetBinLowEdge(i),hDataMass->GetBinLowEdge(i+1));
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
   hDataWidth->Draw("e");
   hDataWidth->SetMinimum(0.5);
   hDataWidth->SetMaximum(4.5);
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
      //TH1D *h = ZmassRoofit_single(0,200,hDataMass->GetBinLowEdge(i),hDataMass->GetBinLowEdge(i+1));
      TH1D *h = ZmassRoofit_single(0,2000,i-1,i); 
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
   hDataWidth->Draw("e");
   hDataWidth->SetMinimum(-2);
   hDataWidth->SetMaximum(6);
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


void ZmassRoofit(){

   style();

   ZmassRoofit_single();
   loop();
   loopHiBin();
   //ZmassRoofit_single(0,20);

}
