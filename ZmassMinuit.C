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
#include <TSystem.h>
#include <TLatex.h>
#include <TChain.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TMinuit.h>

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

TNtuple* n_Zmass = NULL;

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

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
    double nb     = par[0];
    double ns     = par[2];
    f = 2.*(ns+nb);
    
    //Double_t mass[1];
    double *mass = (double*) n_Zmass->GetArgs();
    //mass[0] = n_Zmass->GetArgs();
    for (int i=0;i<n_Zmass->GetEntries();i++) {
        n_Zmass->GetEntry(i);
        
        double L = bwfun(mass,par);
        
        if (L>0.) f -= 2.*log(L);
        else {f = HUGE; return; }
    }
}

TH1D* ZmassMinuit_single(double ptL=0,double ptH=2000,int centL=0,int centH=4)
{

   int binsize = 40;
   //if((centL==3&&centH==4)||ptL>35) binsize = 30;

   TH1D *hData = new TH1D("hData","",binsize,81.2,101.2);
   TH1D *hMC = new TH1D("hMC","",binsize,81.2,101.2);

   n_Zmass = new TNtuple ( "n_Zmass", "Invariant Mass of #mu^{+}#mu^{-}", "mass");

   TCanvas *c = new TCanvas("c","",800,800);
   
   TChain *tMC = new TChain("t");
   tMC->Add("/eos/cms/store/group/phys_heavyions_ops/pchou/MC/*.root?#t");
   
   TChain *tData = new TChain("t");
   tData->Add(Form("/eos/cms/store/group/phys_heavyions_ops/pchou/%s/*.root?#t",typeofdata));

   vector<double>* tagzMassData = 0;
   vector<double>* tagzPtData = 0;
   Float_t taghiHFData;
   tData->SetBranchAddress("zMass",&tagzMassData);
   tData->SetBranchAddress("zPt", &tagzPtData);
   tData->SetBranchAddress("hiHF",&taghiHFData);

   int no_of_eventsData = (int)tData->GetEntries();
   double TotalEventsData = 0;

   TH1D *h = new TH1D("h","",100,0,1);
   
   tData->Draw("zMass>>hData",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   hData->Sumw2();

   tMC->Draw("zMass>>hMC",Form("zPt>%f&&zPt<%f&&hiHF<=%.4f&&hiHF>%.4f",ptL,ptH,hf_diff[centL],hf_diff[centH]));

   hMC->Sumw2();

   std::cout<<"Start filling Data"<<std::endl;

   for(int idx=0;idx<no_of_eventsData;idx++)
   { 
     tData->GetEntry(idx);  
     int zSize = (int) tagzMassData->size();
     if(zSize==0) continue;
     for(int ii=0; ii<zSize;ii++){
          if((*tagzPtData)[ii]>ptL&&(*tagzPtData)[ii]<ptH&&taghiHFData>hf_diff[centH]&&taghiHFData<=hf_diff[centL]){
              n_Zmass->Fill((*tagzMassData)[ii]);
              TotalEventsData++;
          }
       }
   }

   int countD = hData->GetEntries();
   std::cout<<"Data = "<<countD<<std::endl;
   int countM = hMC->GetEntries();
   std::cout<<"MC = "<<countM<<std::endl;


   hMC->SetMarkerStyle(24);

   hData->SetMarkerColor(kBlack);
   hData->SetLineColor(kBlack);
   
   hMC->SetMarkerColor(kRed);
   hMC->SetLineColor(kRed);

   hMC -> SetLineWidth(2);
   hData->SetLineWidth(2);

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

   hMC->SetMinimum(0);
   hData->SetMinimum(0);

   
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

   TMinuit *gMinuit = new TMinuit(7);
   gMinuit->SetFCN(fcn);

   gMinuit->DefineParameter(0, "Mean",   91., 1.,  80.,  120.);
   gMinuit->DefineParameter(1, "Width",  4.6, 1.,  0.1,   10.);
   gMinuit->DefineParameter(2, "Area",    3., 1.,  0.1,  300.);
   gMinuit->DefineParameter(3, "GSigma", 0.1, 1.,   0.,   30.);
   gMinuit->DefineParameter(4, "BkgArea",  0, 1., -100,  100.);
   gMinuit->DefineParameter(5, "ExpCont",  0, 1., -100,  100.);
   gMinuit->DefineParameter(6, "ExpShift", 0, 1., -100,  100.);
   
   gMinuit->Command("MINI");
   gMinuit->Command("MINI");
   gMinuit->Command("MINOS");  

   double par[7],err[7];
   for(int i=0;i<7;i++) gMinuit->GetParameter(i,par[i],err[i]);

/*
   hData->Fit("f","LL");
   hData->Fit("f","");
   hData->Fit("f","LL m");
   hData->Fit("f","");
   hData->Fit("f","LL m");
   hData->Fit("f","");
   hData->Fit("f","");
   hData->Fit("f","");
   hData->Fit("f","");
 */

   f->SetParameters(par[0],par[1],par[2],par[3],par[4],par[5],par[6],0);
   //f->SetParErrors(err[0],err[1],err[2],err[3],err[4],err[5],err[6],1);
   f->SetParErrors(err);


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

   gSystem->Exec(Form("mkdir -p figs/mass/minuit/%s",typeofdata));

   c->SaveAs(Form("figs/mass/minuit/%s/Zmass_%s_%.0f_%.0f_%.0f_%.0f.png",typeofdata,typeofdata,ptL,ptH,cent_diff[centL],cent_diff[centH])); 
   c->Clear();


   //hData->Reset();
   //hMC->Reset();

   n_Zmass=NULL;
   std::cout<<"n_Zmass=NULL;"<<std::endl;

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
      TH1D *h = ZmassMinuit_single(hDataMass->GetBinLowEdge(i),hDataMass->GetBinLowEdge(i+1));
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

   c1->SaveAs(Form("figs/mass/minuit/%s/Zmass_%s_loop.png",typeofdata,typeofdata)); 
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
   hDataWidth->SetMinimum(0);
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

   c1->SaveAs(Form("figs/mass/minuit/%s/ZmassWidth_%s_loop.png",typeofdata,typeofdata)); 
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
      //TH1D *h = ZmassMinuit_single(0,200,hDataMass->GetBinLowEdge(i),hDataMass->GetBinLowEdge(i+1));
      TH1D *h = ZmassMinuit_single(0,2000,i-1,i); 
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

   c2->SaveAs(Form("figs/mass/minuit/%s/Zmass_%s_loopHiBin.png",typeofdata,typeofdata)); 
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
   hDataWidth->SetMinimum(0);
   hDataWidth->SetMaximum(8);
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

   c2->SaveAs(Form("figs/mass/minuit/%s/ZmassWidth_%s_loopHiBin.png",typeofdata,typeofdata)); 
   c2->Clear();
}  


void ZmassMinuit(){

   style();

   ZmassMinuit_single();
   //loop();
   //loopHiBin();
   //ZmassMinuit_single(0,20);

}
