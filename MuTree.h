//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  1 08:08:05 2022 by ROOT version 6.24/06
// from TTree HLTMuTree/HLTMuTree
// found on file: DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8_merged.root
//////////////////////////////////////////////////////////

#ifndef MuTree_h
#define MuTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MuTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Event;
   Int_t           Lumi;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Int_t           Gen_nptl;
   Int_t           Gen_pid[17];   //[Gen_nptl]
   Int_t           Gen_mom[17];   //[Gen_nptl]
   Int_t           Gen_status[17];   //[Gen_nptl]
   Float_t         Gen_p[17];   //[Gen_nptl]
   Float_t         Gen_pt[17];   //[Gen_nptl]
   Float_t         Gen_eta[17];   //[Gen_nptl]
   Float_t         Gen_phi[17];   //[Gen_nptl]
   Int_t           Glb_nptl;
   Int_t           Glb_charge[7];   //[Glb_nptl]
   Float_t         Glb_p[7];   //[Glb_nptl]
   Float_t         Glb_pt[7];   //[Glb_nptl]
   Float_t         Glb_eta[7];   //[Glb_nptl]
   Float_t         Glb_phi[7];   //[Glb_nptl]
   Float_t         Glb_dxy[7];   //[Glb_nptl]
   Float_t         Glb_dz[7];   //[Glb_nptl]
   Int_t           Glb_nValMuHits[7];   //[Glb_nptl]
   Int_t           Glb_nValTrkHits[7];   //[Glb_nptl]
   Int_t           Glb_nValPixHits[7];   //[Glb_nptl]
   Int_t           Glb_trkLayerWMeas[7];   //[Glb_nptl]
   Int_t           Glb_nMatchedStations[7];   //[Glb_nptl]
   Int_t           Glb_nTrkFound[7];   //[Glb_nptl]
   Float_t         Glb_glbChi2_ndof[7];   //[Glb_nptl]
   Float_t         Glb_trkChi2_ndof[7];   //[Glb_nptl]
   Int_t           Glb_pixLayerWMeas[7];   //[Glb_nptl]
   Float_t         Glb_trkDxy[7];   //[Glb_nptl]
   Float_t         Glb_trkDz[7];   //[Glb_nptl]
   Int_t           Sta_nptl;
   Int_t           Sta_charge[9];   //[Sta_nptl]
   Float_t         Sta_p[9];   //[Sta_nptl]
   Float_t         Sta_pt[9];   //[Sta_nptl]
   Float_t         Sta_eta[9];   //[Sta_nptl]
   Float_t         Sta_phi[9];   //[Sta_nptl]
   Float_t         Sta_dxy[9];   //[Sta_nptl]
   Float_t         Sta_dz[9];   //[Sta_nptl]
   Int_t           Glb_isArbitrated[7];   //[Glb_nptl]
   Int_t           Di_npair;
   Float_t         Di_vProb[21];   //[Di_npair]
   Float_t         Di_mass[21];   //[Di_npair]
   Float_t         Di_e[21];   //[Di_npair]
   Float_t         Di_pt[21];   //[Di_npair]
   Float_t         Di_pt1[21];   //[Di_npair]
   Float_t         Di_pt2[21];   //[Di_npair]
   Float_t         Di_eta[21];   //[Di_npair]
   Float_t         Di_eta1[21];   //[Di_npair]
   Float_t         Di_eta2[21];   //[Di_npair]
   Float_t         Di_rapidity[21];   //[Di_npair]
   Float_t         Di_phi[21];   //[Di_npair]
   Float_t         Di_phi1[21];   //[Di_npair]
   Float_t         Di_phi2[21];   //[Di_npair]
   Int_t           Di_charge[21];   //[Di_npair]
   Int_t           Di_charge1[21];   //[Di_npair]
   Int_t           Di_charge2[21];   //[Di_npair]
   Int_t           Di_isArb1[21];   //[Di_npair]
   Int_t           Di_isArb2[21];   //[Di_npair]
   Int_t           Di_nTrkHit1[21];   //[Di_npair]
   Int_t           Di_nTrkHit2[21];   //[Di_npair]
   Int_t           Di_nMuHit1[21];   //[Di_npair]
   Int_t           Di_nMuHit2[21];   //[Di_npair]
   Int_t           Di_nTrkLayers1[21];   //[Di_npair]
   Int_t           Di_nTrkLayers2[21];   //[Di_npair]
   Int_t           Di_nPixHit1[21];   //[Di_npair]
   Int_t           Di_nPixHit2[21];   //[Di_npair]
   Int_t           Di_nMatchedStations1[21];   //[Di_npair]
   Int_t           Di_nMatchedStations2[21];   //[Di_npair]
   Float_t         Di_trkChi2_1[21];   //[Di_npair]
   Float_t         Di_trkChi2_2[21];   //[Di_npair]
   Float_t         Di_glbChi2_1[21];   //[Di_npair]
   Float_t         Di_glbChi2_2[21];   //[Di_npair]
   Float_t         Di_dxy1[21];   //[Di_npair]
   Float_t         Di_dxy2[21];   //[Di_npair]
   Float_t         Di_dz1[21];   //[Di_npair]
   Float_t         Di_dz2[21];   //[Di_npair]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_Gen_nptl;   //!
   TBranch        *b_Gen_pid;   //!
   TBranch        *b_Gen_mom;   //!
   TBranch        *b_Gen_status;   //!
   TBranch        *b_Gen_p;   //!
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_phi;   //!
   TBranch        *b_Glb_nptl;   //!
   TBranch        *b_Glb_charge;   //!
   TBranch        *b_Glb_p;   //!
   TBranch        *b_Glb_pt;   //!
   TBranch        *b_Glb_eta;   //!
   TBranch        *b_Glb_phi;   //!
   TBranch        *b_Glb_dxy;   //!
   TBranch        *b_Glb_dz;   //!
   TBranch        *b_Glb_nValMuHits;   //!
   TBranch        *b_Glb_nValTrkHits;   //!
   TBranch        *b_Glb_nValPixHits;   //!
   TBranch        *b_Glb_trkLayerWMeas;   //!
   TBranch        *b_Glb_nMatchedStations;   //!
   TBranch        *b_Glb_nTrkFound;   //!
   TBranch        *b_Glb_glbChi2_ndof;   //!
   TBranch        *b_Glb_trkChi2_ndof;   //!
   TBranch        *b_Glb_pixLayerWMeas;   //!
   TBranch        *b_Glb_trkDxy;   //!
   TBranch        *b_Glb_trkDz;   //!
   TBranch        *b_Sta_nptl;   //!
   TBranch        *b_Sta_charge;   //!
   TBranch        *b_Sta_p;   //!
   TBranch        *b_Sta_pt;   //!
   TBranch        *b_Sta_eta;   //!
   TBranch        *b_Sta_phi;   //!
   TBranch        *b_Sta_dxy;   //!
   TBranch        *b_Sta_dz;   //!
   TBranch        *b_Glb_isArbitrated;   //!
   TBranch        *b_Di_npair;   //!
   TBranch        *b_Di_vProb;   //!
   TBranch        *b_Di_mass;   //!
   TBranch        *b_Di_e;   //!
   TBranch        *b_Di_pt;   //!
   TBranch        *b_Di_pt1;   //!
   TBranch        *b_Di_pt2;   //!
   TBranch        *b_Di_eta;   //!
   TBranch        *b_Di_eta1;   //!
   TBranch        *b_Di_eta2;   //!
   TBranch        *b_Di_rapidity;   //!
   TBranch        *b_Di_phi;   //!
   TBranch        *b_Di_phi1;   //!
   TBranch        *b_Di_phi2;   //!
   TBranch        *b_Di_charge;   //!
   TBranch        *b_Di_charge1;   //!
   TBranch        *b_Di_charge2;   //!
   TBranch        *b_Di_isArb1;   //!
   TBranch        *b_Di_isArb2;   //!
   TBranch        *b_Di_nTrkHit1;   //!
   TBranch        *b_Di_nTrkHit2;   //!
   TBranch        *b_Di_nMuHit1;   //!
   TBranch        *b_Di_nMuHit2;   //!
   TBranch        *b_Di_nTrkLayers1;   //!
   TBranch        *b_Di_nTrkLayers2;   //!
   TBranch        *b_Di_nPixHit1;   //!
   TBranch        *b_Di_nPixHit2;   //!
   TBranch        *b_Di_nMatchedStations1;   //!
   TBranch        *b_Di_nMatchedStations2;   //!
   TBranch        *b_Di_trkChi2_1;   //!
   TBranch        *b_Di_trkChi2_2;   //!
   TBranch        *b_Di_glbChi2_1;   //!
   TBranch        *b_Di_glbChi2_2;   //!
   TBranch        *b_Di_dxy1;   //!
   TBranch        *b_Di_dxy2;   //!
   TBranch        *b_Di_dz1;   //!
   TBranch        *b_Di_dz2;   //!

   MuTree(TTree *tree=0);
   virtual ~MuTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MuTree_cxx
MuTree::MuTree(TTree *tree) : fChain(0) 
{

}

MuTree::~MuTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MuTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_run);
   fChain->SetBranchAddress("Event", &Event, &b_event);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_lumi);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("Gen_nptl", &Gen_nptl, &b_Gen_nptl);
   fChain->SetBranchAddress("Gen_pid", Gen_pid, &b_Gen_pid);
   fChain->SetBranchAddress("Gen_mom", Gen_mom, &b_Gen_mom);
   fChain->SetBranchAddress("Gen_status", Gen_status, &b_Gen_status);
   fChain->SetBranchAddress("Gen_p", Gen_p, &b_Gen_p);
   fChain->SetBranchAddress("Gen_pt", Gen_pt, &b_Gen_pt);
   fChain->SetBranchAddress("Gen_eta", Gen_eta, &b_Gen_eta);
   fChain->SetBranchAddress("Gen_phi", Gen_phi, &b_Gen_phi);
   fChain->SetBranchAddress("Glb_nptl", &Glb_nptl, &b_Glb_nptl);
   fChain->SetBranchAddress("Glb_charge", Glb_charge, &b_Glb_charge);
   fChain->SetBranchAddress("Glb_p", Glb_p, &b_Glb_p);
   fChain->SetBranchAddress("Glb_pt", Glb_pt, &b_Glb_pt);
   fChain->SetBranchAddress("Glb_eta", Glb_eta, &b_Glb_eta);
   fChain->SetBranchAddress("Glb_phi", Glb_phi, &b_Glb_phi);
   fChain->SetBranchAddress("Glb_dxy", Glb_dxy, &b_Glb_dxy);
   fChain->SetBranchAddress("Glb_dz", Glb_dz, &b_Glb_dz);
   fChain->SetBranchAddress("Glb_nValMuHits", Glb_nValMuHits, &b_Glb_nValMuHits);
   fChain->SetBranchAddress("Glb_nValTrkHits", Glb_nValTrkHits, &b_Glb_nValTrkHits);
   fChain->SetBranchAddress("Glb_nValPixHits", Glb_nValPixHits, &b_Glb_nValPixHits);
   fChain->SetBranchAddress("Glb_trkLayerWMeas", Glb_trkLayerWMeas, &b_Glb_trkLayerWMeas);
   fChain->SetBranchAddress("Glb_nMatchedStations", Glb_nMatchedStations, &b_Glb_nMatchedStations);
   fChain->SetBranchAddress("Glb_nTrkFound", Glb_nTrkFound, &b_Glb_nTrkFound);
   fChain->SetBranchAddress("Glb_glbChi2_ndof", Glb_glbChi2_ndof, &b_Glb_glbChi2_ndof);
   fChain->SetBranchAddress("Glb_trkChi2_ndof", Glb_trkChi2_ndof, &b_Glb_trkChi2_ndof);
   fChain->SetBranchAddress("Glb_pixLayerWMeas", Glb_pixLayerWMeas, &b_Glb_pixLayerWMeas);
   fChain->SetBranchAddress("Glb_trkDxy", Glb_trkDxy, &b_Glb_trkDxy);
   fChain->SetBranchAddress("Glb_trkDz", Glb_trkDz, &b_Glb_trkDz);
   fChain->SetBranchAddress("Sta_nptl", &Sta_nptl, &b_Sta_nptl);
   fChain->SetBranchAddress("Sta_charge", Sta_charge, &b_Sta_charge);
   fChain->SetBranchAddress("Sta_p", Sta_p, &b_Sta_p);
   fChain->SetBranchAddress("Sta_pt", Sta_pt, &b_Sta_pt);
   fChain->SetBranchAddress("Sta_eta", Sta_eta, &b_Sta_eta);
   fChain->SetBranchAddress("Sta_phi", Sta_phi, &b_Sta_phi);
   fChain->SetBranchAddress("Sta_dxy", Sta_dxy, &b_Sta_dxy);
   fChain->SetBranchAddress("Sta_dz", Sta_dz, &b_Sta_dz);
   fChain->SetBranchAddress("Glb_isArbitrated", Glb_isArbitrated, &b_Glb_isArbitrated);
   fChain->SetBranchAddress("Di_npair", &Di_npair, &b_Di_npair);
   fChain->SetBranchAddress("Di_vProb", Di_vProb, &b_Di_vProb);
   fChain->SetBranchAddress("Di_mass", Di_mass, &b_Di_mass);
   fChain->SetBranchAddress("Di_e", Di_e, &b_Di_e);
   fChain->SetBranchAddress("Di_pt", Di_pt, &b_Di_pt);
   fChain->SetBranchAddress("Di_pt1", Di_pt1, &b_Di_pt1);
   fChain->SetBranchAddress("Di_pt2", Di_pt2, &b_Di_pt2);
   fChain->SetBranchAddress("Di_eta", Di_eta, &b_Di_eta);
   fChain->SetBranchAddress("Di_eta1", Di_eta1, &b_Di_eta1);
   fChain->SetBranchAddress("Di_eta2", Di_eta2, &b_Di_eta2);
   fChain->SetBranchAddress("Di_rapidity", Di_rapidity, &b_Di_rapidity);
   fChain->SetBranchAddress("Di_phi", Di_phi, &b_Di_phi);
   fChain->SetBranchAddress("Di_phi1", Di_phi1, &b_Di_phi1);
   fChain->SetBranchAddress("Di_phi2", Di_phi2, &b_Di_phi2);
   fChain->SetBranchAddress("Di_charge", Di_charge, &b_Di_charge);
   fChain->SetBranchAddress("Di_charge1", Di_charge1, &b_Di_charge1);
   fChain->SetBranchAddress("Di_charge2", Di_charge2, &b_Di_charge2);
   fChain->SetBranchAddress("Di_isArb1", Di_isArb1, &b_Di_isArb1);
   fChain->SetBranchAddress("Di_isArb2", Di_isArb2, &b_Di_isArb2);
   fChain->SetBranchAddress("Di_nTrkHit1", Di_nTrkHit1, &b_Di_nTrkHit1);
   fChain->SetBranchAddress("Di_nTrkHit2", Di_nTrkHit2, &b_Di_nTrkHit2);
   fChain->SetBranchAddress("Di_nMuHit1", Di_nMuHit1, &b_Di_nMuHit1);
   fChain->SetBranchAddress("Di_nMuHit2", Di_nMuHit2, &b_Di_nMuHit2);
   fChain->SetBranchAddress("Di_nTrkLayers1", Di_nTrkLayers1, &b_Di_nTrkLayers1);
   fChain->SetBranchAddress("Di_nTrkLayers2", Di_nTrkLayers2, &b_Di_nTrkLayers2);
   fChain->SetBranchAddress("Di_nPixHit1", Di_nPixHit1, &b_Di_nPixHit1);
   fChain->SetBranchAddress("Di_nPixHit2", Di_nPixHit2, &b_Di_nPixHit2);
   fChain->SetBranchAddress("Di_nMatchedStations1", Di_nMatchedStations1, &b_Di_nMatchedStations1);
   fChain->SetBranchAddress("Di_nMatchedStations2", Di_nMatchedStations2, &b_Di_nMatchedStations2);
   fChain->SetBranchAddress("Di_trkChi2_1", Di_trkChi2_1, &b_Di_trkChi2_1);
   fChain->SetBranchAddress("Di_trkChi2_2", Di_trkChi2_2, &b_Di_trkChi2_2);
   fChain->SetBranchAddress("Di_glbChi2_1", Di_glbChi2_1, &b_Di_glbChi2_1);
   fChain->SetBranchAddress("Di_glbChi2_2", Di_glbChi2_2, &b_Di_glbChi2_2);
   fChain->SetBranchAddress("Di_dxy1", Di_dxy1, &b_Di_dxy1);
   fChain->SetBranchAddress("Di_dxy2", Di_dxy2, &b_Di_dxy2);
   fChain->SetBranchAddress("Di_dz1", Di_dz1, &b_Di_dz1);
   fChain->SetBranchAddress("Di_dz2", Di_dz2, &b_Di_dz2);
   Notify();
}

Bool_t MuTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuTree_cxx
