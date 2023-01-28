//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  1 15:36:27 2022 by ROOT version 6.24/06
// from TTree HiTree/
// found on file: DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8_merged.root
//////////////////////////////////////////////////////////

#ifndef HiTree_h
#define HiTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "utility"
#include "utility"
#include "vector"
#include "vector"

class HiTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       evt;
   UInt_t          lumi;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         Npart;
   Float_t         Ncoll;
   Float_t         Nhard;
   Float_t         phi0;
   Float_t         b;
   Int_t           Ncharged;
   Int_t           NchargedMR;
   Float_t         MeanPt;
   Float_t         MeanPtMR;
   Float_t         EtMR;
   Int_t           NchargedPtCut;
   Int_t           NchargedPtCutMR;
   Int_t           ProcessID;
   Float_t         pthat;
   Float_t         weight;
   Float_t         alphaQCD;
   Float_t         alphaQED;
   Float_t         qScale;
   Int_t           nMEPartons;
   Int_t           nMEPartonsFiltered;
   vector<float>   *ttbar_w;
   vector<int>     *npus;
   vector<float>   *tnpus;
   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiHFECut;
   Float_t         hiHFECutPlus;
   Float_t         hiHFECutMinus;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
   Float_t         hiZDC;
   Float_t         hiZDCplus;
   Float_t         hiZDCminus;
   Float_t         hiHFhit;
   Float_t         hiHFhitPlus;
   Float_t         hiHFhitMinus;
   Float_t         hiET;
   Float_t         hiEE;
   Float_t         hiEB;
   Float_t         hiEEplus;
   Float_t         hiEEminus;
   Int_t           hiNpix;
   Int_t           hiNpixPlus;
   Int_t           hiNpixMinus;
   Int_t           hiNpixelTracks;
   Int_t           hiNpixelTracksPlus;
   Int_t           hiNpixelTracksMinus;
   Int_t           hiNtracks;
   Int_t           hiNtracksPtCut;
   Int_t           hiNtracksEtaCut;
   Int_t           hiNtracksEtaPtCut;
   Int_t           hiNevtPlane;
   Float_t         hiEvtPlanes[24];   //[hiNevtPlane]
   Int_t           numMinHFTower2;
   Int_t           numMinHFTower3;
   Int_t           numMinHFTower4;
   Int_t           numMinHFTower5;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_Npart;   //!
   TBranch        *b_Ncoll;   //!
   TBranch        *b_Nhard;   //!
   TBranch        *b_NPhi0;   //!
   TBranch        *b_b;   //!
   TBranch        *b_Ncharged;   //!
   TBranch        *b_NchargedMR;   //!
   TBranch        *b_MeanPt;   //!
   TBranch        *b_MeanPtMR;   //!
   TBranch        *b_EtMR;   //!
   TBranch        *b_NchargedPtCut;   //!
   TBranch        *b_NchargedPtCutMR;   //!
   TBranch        *b_ProcessID;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_alphaQCD;   //!
   TBranch        *b_alphaQED;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_nMEPartons;   //!
   TBranch        *b_nMEPartonsFiltered;   //!
   TBranch        *b_ttbar_w;   //!
   TBranch        *b_npus;   //!
   TBranch        *b_tnpus;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_hiHFplus;   //!
   TBranch        *b_hiHFminus;   //!
   TBranch        *b_hiHFECut;   //!
   TBranch        *b_hiHFECutPlus;   //!
   TBranch        *b_hiHFECutMinus;   //!
   TBranch        *b_hiHFplusEta4;   //!
   TBranch        *b_hiHFminusEta4;   //!
   TBranch        *b_hiZDC;   //!
   TBranch        *b_hiZDCplus;   //!
   TBranch        *b_hiZDCminus;   //!
   TBranch        *b_hiHFhit;   //!
   TBranch        *b_hiHFhitPlus;   //!
   TBranch        *b_hiHFhitMinus;   //!
   TBranch        *b_hiET;   //!
   TBranch        *b_hiEE;   //!
   TBranch        *b_hiEB;   //!
   TBranch        *b_hiEEplus;   //!
   TBranch        *b_hiEEminus;   //!
   TBranch        *b_hiNpix;   //!
   TBranch        *b_hiNpixPlus;   //!
   TBranch        *b_hiNpixMinus;   //!
   TBranch        *b_hiNpixelTracks;   //!
   TBranch        *b_hiNpixelTracksPlus;   //!
   TBranch        *b_hiNpixelTracksMinus;   //!
   TBranch        *b_hiNtracks;   //!
   TBranch        *b_hiNtracksPtCut;   //!
   TBranch        *b_hiNtracksEtaCut;   //!
   TBranch        *b_hiNtracksEtaPtCut;   //!
   TBranch        *b_hiNevtPlane;   //!
   TBranch        *b_hiEvtPlanes;   //!
   TBranch        *b_numMinHFTower2;   //!
   TBranch        *b_numMinHFTower3;   //!
   TBranch        *b_numMinHFTower4;   //!
   TBranch        *b_numMinHFTower5;   //!

   HiTree(TTree *tree=0);
   virtual ~HiTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HiTree_cxx
HiTree::HiTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
}

HiTree::~HiTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HiTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HiTree::LoadTree(Long64_t entry)
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

void HiTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ttbar_w = 0;
   npus = 0;
   tnpus = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("Npart", &Npart, &b_Npart);
   fChain->SetBranchAddress("Ncoll", &Ncoll, &b_Ncoll);
   fChain->SetBranchAddress("Nhard", &Nhard, &b_Nhard);
   fChain->SetBranchAddress("phi0", &phi0, &b_NPhi0);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("Ncharged", &Ncharged, &b_Ncharged);
   fChain->SetBranchAddress("NchargedMR", &NchargedMR, &b_NchargedMR);
   fChain->SetBranchAddress("MeanPt", &MeanPt, &b_MeanPt);
   fChain->SetBranchAddress("MeanPtMR", &MeanPtMR, &b_MeanPtMR);
   fChain->SetBranchAddress("EtMR", &EtMR, &b_EtMR);
   fChain->SetBranchAddress("NchargedPtCut", &NchargedPtCut, &b_NchargedPtCut);
   fChain->SetBranchAddress("NchargedPtCutMR", &NchargedPtCutMR, &b_NchargedPtCutMR);
   fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("alphaQCD", &alphaQCD, &b_alphaQCD);
   fChain->SetBranchAddress("alphaQED", &alphaQED, &b_alphaQED);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("nMEPartons", &nMEPartons, &b_nMEPartons);
   fChain->SetBranchAddress("nMEPartonsFiltered", &nMEPartonsFiltered, &b_nMEPartonsFiltered);
   fChain->SetBranchAddress("ttbar_w", &ttbar_w, &b_ttbar_w);
   fChain->SetBranchAddress("npus", &npus, &b_npus);
   fChain->SetBranchAddress("tnpus", &tnpus, &b_tnpus);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   fChain->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
   fChain->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
   fChain->SetBranchAddress("hiHFECut", &hiHFECut, &b_hiHFECut);
   fChain->SetBranchAddress("hiHFECutPlus", &hiHFECutPlus, &b_hiHFECutPlus);
   fChain->SetBranchAddress("hiHFECutMinus", &hiHFECutMinus, &b_hiHFECutMinus);
   fChain->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
   fChain->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);
   fChain->SetBranchAddress("hiZDC", &hiZDC, &b_hiZDC);
   fChain->SetBranchAddress("hiZDCplus", &hiZDCplus, &b_hiZDCplus);
   fChain->SetBranchAddress("hiZDCminus", &hiZDCminus, &b_hiZDCminus);
   fChain->SetBranchAddress("hiHFhit", &hiHFhit, &b_hiHFhit);
   fChain->SetBranchAddress("hiHFhitPlus", &hiHFhitPlus, &b_hiHFhitPlus);
   fChain->SetBranchAddress("hiHFhitMinus", &hiHFhitMinus, &b_hiHFhitMinus);
   fChain->SetBranchAddress("hiET", &hiET, &b_hiET);
   fChain->SetBranchAddress("hiEE", &hiEE, &b_hiEE);
   fChain->SetBranchAddress("hiEB", &hiEB, &b_hiEB);
   fChain->SetBranchAddress("hiEEplus", &hiEEplus, &b_hiEEplus);
   fChain->SetBranchAddress("hiEEminus", &hiEEminus, &b_hiEEminus);
   fChain->SetBranchAddress("hiNpix", &hiNpix, &b_hiNpix);
   fChain->SetBranchAddress("hiNpixPlus", &hiNpixPlus, &b_hiNpixPlus);
   fChain->SetBranchAddress("hiNpixMinus", &hiNpixMinus, &b_hiNpixMinus);
   fChain->SetBranchAddress("hiNpixelTracks", &hiNpixelTracks, &b_hiNpixelTracks);
   fChain->SetBranchAddress("hiNpixelTracksPlus", &hiNpixelTracksPlus, &b_hiNpixelTracksPlus);
   fChain->SetBranchAddress("hiNpixelTracksMinus", &hiNpixelTracksMinus, &b_hiNpixelTracksMinus);
   fChain->SetBranchAddress("hiNtracks", &hiNtracks, &b_hiNtracks);
   fChain->SetBranchAddress("hiNtracksPtCut", &hiNtracksPtCut, &b_hiNtracksPtCut);
   fChain->SetBranchAddress("hiNtracksEtaCut", &hiNtracksEtaCut, &b_hiNtracksEtaCut);
   fChain->SetBranchAddress("hiNtracksEtaPtCut", &hiNtracksEtaPtCut, &b_hiNtracksEtaPtCut);
   fChain->SetBranchAddress("hiNevtPlane", &hiNevtPlane, &b_hiNevtPlane);
   fChain->SetBranchAddress("hiEvtPlanes", hiEvtPlanes, &b_hiEvtPlanes);
   fChain->SetBranchAddress("numMinHFTower2", &numMinHFTower2, &b_numMinHFTower2);
   fChain->SetBranchAddress("numMinHFTower3", &numMinHFTower3, &b_numMinHFTower3);
   fChain->SetBranchAddress("numMinHFTower4", &numMinHFTower4, &b_numMinHFTower4);
   fChain->SetBranchAddress("numMinHFTower5", &numMinHFTower5, &b_numMinHFTower5);
   Notify();
}

Bool_t HiTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HiTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HiTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HiTree_cxx
