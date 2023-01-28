#include <TFile.h>
#include <TTree.h>
#include <iostream>

#include "GenParticle.C"
#include "MuTree.C"
#include "PbPbTracks.C"
#include "HiTree.C"


class HiForest
{
   public:
   TFile *hiForestFile;
   long nEntries; 
  
   bool doGenParticle;
   bool doMuTree;
   bool doPbPbTracks;
   bool doHiTree;


   TTree *tGenParticle;
   TTree *tMuTree;
   TTree *tPbPbTracks;
   TTree *tHiTree;
   
   GenParticle genParticle;
   MuTree muTree;
   PbPbTracks tracks;
   HiTree hi;
   
   
   HiForest(string fileName="hiForest.root"){
       hiForestFile=new TFile(fileName.c_str());
       cout <<"Opened "<<fileName.c_str()<<endl;
       cout <<hiForestFile<<endl;
       nEntries=0;
   };
   
   
   ~HiForest(){};

   void Init();
   int GetEntry(int i);
   long GetEntries() {
      return nEntries;
   };
};

void HiForest::Init()
{
   if (doGenParticle) {
       tGenParticle = (TTree*) hiForestFile->Get("HiGenParticleAna/hi");
       if (tGenParticle==0) cout <<"Error! GenParticle"<<endl;
       else {       
          genParticle.Init(tGenParticle);
          nEntries=tGenParticle->GetEntries();
       }
   }   
   if (doMuTree) {
       tMuTree = (TTree*) hiForestFile->Get("hltMuTree/HLTMuTree");
       if (tMuTree==0) cout <<"Error! MuTree"<<endl;
       else {       
          muTree.Init(tMuTree);
          nEntries=tMuTree->GetEntries();
       }
   }   
   if (doPbPbTracks) {
       tPbPbTracks = (TTree*) hiForestFile->Get("PbPbTracks/trackTree");
       if (tPbPbTracks==0) cout <<"Error! PbPbTrack"<<endl;
       else {       
          tracks.Init(tPbPbTracks);
          nEntries=tPbPbTracks->GetEntries();
       }
   }   
   if (doHiTree) {
       tHiTree = (TTree*) hiForestFile->Get("hiEvtAnalyzer/HiTree");
       if (tHiTree==0) cout <<"Error! HiTree"<<endl;
       else {       
          hi.Init(tHiTree);
          nEntries=tHiTree->GetEntries();
       }
   }   
}

int HiForest::GetEntry(int i)
{
   if (doGenParticle) {
       tGenParticle->GetEntry(i);
   }
   if (doMuTree) {
       tMuTree->GetEntry(i);
   }
   if (doPbPbTracks) {
       tPbPbTracks->GetEntry(i);
   }
   if (doHiTree) {
       tHiTree->GetEntry(i);
   }
   return 1;
}
