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
   vector<double> trackPDFId;
   vector<double> trackEta;
   vector<double> trackPhi;

   vector<double> muEta1;
   vector<double> muEta2;
   vector<double> muPhi1;
   vector<double> muPhi2;
   vector<double> muPt1;
   vector<double> muPt2;   

   vector<double> genMuPt1;
   vector<double> genMuPt2;
   vector<double> genMuEta1;
   vector<double> genMuEta2;
   vector<double> genMuPhi1;
   vector<double> genMuPhi2;
   int hiBin;
   float hiHF;

   
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
      t->Branch("trackPDFId",&trackPDFId);
      t->Branch("trackPhi",&trackPhi);
      t->Branch("trackEta",&trackEta);
      t->Branch("hiBin",&hiBin);
      t->Branch("hiHF",&hiHF);

      t->Branch("muEta1",&muEta1);
      t->Branch("muEta2",&muEta2);
      t->Branch("muPhi1",&muPhi1);
      t->Branch("muPhi2",&muPhi2);
      t->Branch("muPt1",&muPt1);
      t->Branch("muPt2",&muPt2);

      t->Branch("genMuPt1",&genMuPt1);
      t->Branch("genMuPt2",&genMuPt2);
      t->Branch("genMuEta1",&genMuEta1);
      t->Branch("genMuEta2",&genMuEta2);
      t->Branch("genMuPhi1",&genMuPhi1);
      t->Branch("genMuPhi2",&genMuPhi2);
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
      trackPDFId.clear();
      trackPhi.clear();
      trackEta.clear();

      muEta1.clear();
      muEta2.clear();
      muPhi1.clear();
      muPhi2.clear();
      muPt1.clear();
      muPt2.clear();

      genMuPt1.clear();
      genMuPt2.clear();
      genMuEta1.clear();
      genMuEta2.clear();
      genMuPhi1.clear();
      genMuPhi2.clear();

   }
   
};

void Zhadron_singleFile(const char *dirname, TString fname, int genornot, TNtuple *nt, TTree *t, TH2D *h2D ){

   HiForest f(Form("%s/%s",dirname,fname.Data()));
   f.doGenParticle=genornot;
   f.doPbPbTracks=1;
   f.doMuTree=1;
   f.doHiTree=1;
   f.Init();

   ZhadronData data;
   data.setBranch(t);
   data.clear();


   TLorentzVector LgenMu1;
   TLorentzVector LgenMu2;
   TLorentzVector LgenZ;

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
                      data.genZPt.push_back(LgenZ.Pt());
                      data.genZPhi.push_back(LgenZ.Phi());
                      data.genZEta.push_back(LgenZ.Eta());

                      data.genMuPt1.push_back(f.muTree.Gen_pt[igen1]);
                      data.genMuPt2.push_back(f.muTree.Gen_pt[igen2]);
                      data.genMuEta1.push_back(f.muTree.Gen_eta[igen1]);
                      data.genMuEta2.push_back(f.muTree.Gen_eta[igen2]);
                      data.genMuPhi1.push_back(f.muTree.Gen_phi[igen1]);
                      data.genMuPhi2.push_back(f.muTree.Gen_phi[igen2]);
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

   data.muEta1.push_back(f.muTree.Di_eta1[ipair]);
   data.muEta2.push_back(f.muTree.Di_eta2[ipair]);
   data.muPhi1.push_back(f.muTree.Di_phi1[ipair]);
   data.muPhi2.push_back(f.muTree.Di_phi2[ipair]);

   data.muPt1.push_back(f.muTree.Di_pt1[ipair]);
   data.muPt2.push_back(f.muTree.Di_pt2[ipair]);

   data.hiBin=f.hi.hiBin;
   data.hiHF=f.hi.hiHF;
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
      data.trackPDFId.push_back(f.tracks.trkPDFId->at(itrack));
      data.trackPhi.push_back(f.tracks.trkPhi->at(itrack));
      data.trackEta.push_back(f.tracks.trkEta->at(itrack));
     }

     }
     t->Fill();
     data.clear();
   }

}

//void Zhadron(string
//infname="DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8_merged.root", string outfname="/eos/cms/store/group/phys_heavyions_ops/pchou/outputMC.root")
//{
void Zhadron_loop_part(string infname, string outfname="outputMC.root", int genornot=0, int partnum=0,int numofpart=50)
{
   std::cout<<"Start running..."<<std::endl;
/*
   HiForest f(infname.c_str());

   f.doGenParticle=genornot;
   f.doPbPbTracks=1;
   f.doMuTree=1;
   f.doHiTree=1;
   f.Init();
  */ 
   TFile *outfile = new TFile(outfname.c_str(),"recreate");
   
   TNtuple *nt = new TNtuple("ntZ","Z tree","mass:pt:eta:phi");
   TTree *t = new TTree("t","tree");
   
   TH2D *h2D = new TH2D("h2D","",100,-6,6,100,-3.14159,3.14159);
   
   
   
   std::cout<<"Start looping..."<<std::endl;

   //const char *dirname="/eos/cms/store/group/phys_heavyions/chenyi/PbPb2018/Forest/HIDoubleMuon/HIRun2018A-PbPb18_MiniAODv1-v1/MINIAOD/HIDoubleMuon/20221101_ZHadronDoubleMuonWithMuTree/221101_212941/0000/";
   const char *dirname=infname.c_str();
   const char *ext=".root";

   int count=0;
   
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         std::cout << fname.Data() << std::endl;
         if (!file->IsDirectory() && fname.EndsWith(ext)&&count%numofpart==partnum%numofpart) {

            ////HiForest f(fname.Data());
            //HiForest f(Form("%s/%s",dirname,fname.Data()));
            //f.doGenParticle=genornot;
            //f.doPbPbTracks=1;
            //f.doMuTree=1;
            //f.doHiTree=1;
            //f.Init();
            ////cout << fname.Data() << endl;
            Zhadron_singleFile(dirname, fname, genornot, nt, t, h2D );
            
         }
         count++;
      }
   }

   

   h2D->Write();   
   nt->Write();
   outfile->Write();
   outfile->Close();      

}

