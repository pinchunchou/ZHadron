#include "HiForest.h"
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TMath.h>
#include <TLorentzVector.h>

#define M_MU 0.1056583755

double dphi(double phi1,double phi2)
{
 double a=phi1-phi2;
 while (a<-TMath::Pi()) a+=2*TMath::Pi();
 while (a>TMath::Pi()) a-=2*TMath::Pi();
   
// if (a<-TMath::Pi()/2) a=2*TMath::Pi()+a;
 return a;
}


class ZhadronData
{
   public:
   vector<double> zMass;
   vector<double> zEta;
   vector<double> zPhi;
   vector<double> zPt;
   vector<double> genZMass;
   vector<double> genZEta;
   vector<double> genZPhi;
   vector<double> genZPt;
   vector<double> trackDphi;
   vector<double> trackPt;
   vector<double> trackDeta;
   int hiBin;
   
   ZhadronData(){};
   ~ZhadronData(){};
   void setBranch(TTree *t){
      t->Branch("zMass",&zMass);
      t->Branch("zEta",&zEta);
      t->Branch("zPhi",&zPhi);
      t->Branch("zPt",&zPt);
      t->Branch("genZMass",&genZMass);
      t->Branch("genZEta",&genZEta);
      t->Branch("genZPhi",&genZPhi);
      t->Branch("genZPt",&genZPt);
      t->Branch("trackDphi",&trackDphi);
      t->Branch("trackDeta",&trackDeta);
      t->Branch("trackPt",&trackPt);
      t->Branch("hiBin",&hiBin);
   };
   
   void clear() {
      zMass.clear();
      zEta.clear();
      zPhi.clear();
      zPt.clear();
      genZMass.clear();
      genZEta.clear();
      genZPhi.clear();
      genZPt.clear();
      trackDphi.clear();
      trackPt.clear();
      trackDeta.clear();
   }
   
};

//void Zhadron(string
//infname="DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8_merged.root", string outfname="/eos/cms/store/group/phys_heavyions_ops/pchou/outputMC.root")
//{
void Zhadron(string infname="DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8_merged.root", string outfname="outputMC_hydjet.root", int genornot=1)
{
   std::cout<<"Start running..."<<std::endl;

   //infname = "/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbSpring21MiniAOD-mva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/20221017_ZHadronInitialRunMCDY/221017_113040/0000/HiForestMiniAOD_100.root";
   //infname = "/eos/cms/store/group/phys_heavyions_ops/pchou/Zhadron_dyll.root";


   //infname = "/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/HIDoubleMuon/HIRun2018A-PbPb18_MiniAODv1-v1/MINIAOD/HIDoubleMuon/20221101_ZHadronDoubleMuonWithMuTree/221101_212941/0000/HiForestMiniAOD_108.root";
   HiForest f(infname.c_str());
  // HiForest f("/eos/cms/store/group/phys_heavyions_ops/pchou/Zhadron_hydjet.root");
   //HiForest fBkg("hydjet.root");
   //HiForest fBkg("/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/MinBias_Hydjet_Drum5F_2018_5p02TeV/HINPbPbSpring21MiniAOD-NoPUmva98_112X_upgrade2018_realistic_HI_v9-v1/MINIAODSIM/MinBias_Hydjet_Drum5F_2018_5p02TeV/20221017_ZHadronInitialRunMCMB/221017_113057/0000/HiForestMiniAOD_20.root");
   //HiForest fBkg("/eos/cms/store/group/phys_heavyions_ops/pchou/Zhadron_hydjet.root");

   std::cout<<"f..."<<std::endl;
   f.doGenParticle=genornot;
   f.doPbPbTracks=1;
   f.doMuTree=1;
   f.doHiTree=1;
   f.Init();
/*
   std::cout<<"fBkg..."<<std::endl;

   fBkg.doGenParticle=0;
   fBkg.doPbPbTracks=1;
   fBkg.doMuTree=0;
   fBkg.doHiTree=0;
   fBkg.Init();
   */
   TFile *outfile = new TFile(outfname.c_str(),"recreate");
   
   TNtuple *nt = new TNtuple("ntZ","Z tree","mass:pt:eta");
   TTree *t = new TTree("t","tree");
   
   TH2D *h2D = new TH2D("h2D","",100,-6,6,100,-3.14159,3.14159);
   TH2D *h2Dmix = new TH2D("h2Dmix","",100,-6,6,100,-3.14159,3.14159);
   
   ZhadronData data;
   data.setBranch(t);
   data.clear();
  
   
   TLorentzVector LgenMu1;
   TLorentzVector LgenMu2;
   TLorentzVector LgenZ;
   
   std::cout<<"Start looping..."<<std::endl;

   for (int i=0;i<f.GetEntries();i++)
   {
     f.GetEntry(i);

     // display the progress
     if (i%1000==0) cout <<i<<"/"<<f.GetEntries()<<endl;
     // Loop over Gen information (single muons)
     if (f.muTree.Gen_nptl>1) {
        for (int igen1=0;igen1<f.muTree.Gen_nptl;igen1++) {
           if (f.muTree.Gen_mom[igen1]==23) {
              LgenMu1.SetPtEtaPhiM(f.muTree.Gen_pt[igen1],
	                          f.muTree.Gen_eta[igen1],
	                          f.muTree.Gen_phi[igen1],
	                          M_MU);
	      for (int igen2=igen1+1;igen2<f.muTree.Gen_nptl;igen2++) {
                 if (f.muTree.Gen_mom[igen2]==23) {
                    LgenMu2.SetPtEtaPhiM(f.muTree.Gen_pt[igen2],
	                                f.muTree.Gen_eta[igen2],
	                                f.muTree.Gen_phi[igen2],
	                                M_MU);
		    
		    LgenZ=LgenMu1+LgenMu2;
		    data.genZMass.push_back(LgenZ.M());
		    data.genZPt.push_back  (LgenZ.Pt());
		    data.genZPhi.push_back (LgenZ.Phi());
		    data.genZEta.push_back (LgenZ.Eta());
		 }
              }
	   }   
        }
     }

     // loop over RECO information (dimuons)
     for (int ipair=0;ipair<f.muTree.Di_npair;ipair++) {
        //cout <<ipair<<" "<<f.muTree.Di_mass[ipair]<<endl;
	if (f.muTree.Di_charge1[ipair]==f.muTree.Di_charge2[ipair]) continue;
	if (fabs(f.muTree.Di_eta1[ipair])>2.4) continue;
	if (fabs(f.muTree.Di_eta2[ipair])>2.4) continue;
	if (fabs(f.muTree.Di_pt1[ipair])<20) continue;
	if (fabs(f.muTree.Di_pt2[ipair])<20) continue;
	data.zMass.push_back(f.muTree.Di_mass[ipair]);
	data.zEta.push_back(f.muTree.Di_eta[ipair]);
	data.zPhi.push_back(f.muTree.Di_phi[ipair]);
	data.zPt.push_back(f.muTree.Di_pt[ipair]);
	data.hiBin=f.hi.hiBin;
	nt->Fill(f.muTree.Di_mass[ipair],f.muTree.Di_pt[ipair], f.muTree.Di_eta[ipair], f.muTree.Di_phi[ipair]);
     }
     
     if (data.zMass.size()>0&&data.zPt.at(0)>30) {
        for (int itrack=0;itrack<f.tracks.trkPt->size();itrack++) {
	   if (!f.tracks.highPurity) continue;
           double deltaPhi = dphi(data.zPhi.at(0),f.tracks.trkPhi->at(itrack)-TMath::Pi());
	   double deltaEta = fabs(data.zEta.at(0)-f.tracks.trkEta->at(itrack));
	   h2D->Fill(deltaEta,deltaPhi,0.25);
	   h2D->Fill(-deltaEta,deltaPhi,0.25);
	   h2D->Fill(-deltaEta,-deltaPhi,0.25);
	   h2D->Fill(deltaEta,-deltaPhi,0.25);
	   data.trackDphi.push_back(deltaPhi);
	   data.trackDeta.push_back(deltaEta);
	   data.trackPt.push_back(f.tracks.trkPt->at(itrack));
	}
   /*
	fBkg.GetEntry(i+1%fBkg.GetEntries());
	if (i+1==fBkg.GetEntries()) fBkg.GetEntry(0);
        for (int itrack=0;itrack<fBkg.tracks.trkPt->size();itrack++) {
	   if (!fBkg.tracks.highPurity) continue;
           double deltaPhi = dphi(data.zPhi.at(0),fBkg.tracks.trkPhi->at(itrack)-TMath::Pi());
	   double deltaEta = fabs(data.zEta.at(0)-fBkg.tracks.trkEta->at(itrack));
	   h2Dmix->Fill(deltaEta,deltaPhi,0.25);
	   h2Dmix->Fill(-deltaEta,deltaPhi,0.25);
	   h2Dmix->Fill(-deltaEta,-deltaPhi,0.25);
	   h2Dmix->Fill(deltaEta,-deltaPhi,0.25);
//	   data.trackDphi.push_back(deltaPhi);
//	   data.trackDeta.push_back(deltaEta);
//	   data.trackPt.push_back(f.tracks.trkPt->at(itrack));
	}
   */
     }
     t->Fill();
     data.clear();
   }

   h2D->Write();   
   nt->Write();
   outfile->Write();
   outfile->Close();      

}

