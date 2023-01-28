//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  1 08:09:06 2022 by ROOT version 6.24/06
// from TTree trackTree/v1
// found on file: DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8_merged.root
//////////////////////////////////////////////////////////

#ifndef PbPbTracks_h
#define PbPbTracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class PbPbTracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nRun;
   Int_t           nEv;
   Int_t           nLumi;
   vector<float>   *xVtx;
   vector<float>   *yVtx;
   vector<float>   *zVtx;
   vector<float>   *xErrVtx;
   vector<float>   *yErrVtx;
   vector<float>   *zErrVtx;
   vector<float>   *chi2Vtx;
   vector<float>   *ndofVtx;
   vector<bool>    *isFakeVtx;
   vector<int>     *nTracksVtx;
   vector<float>   *ptSumVtx;
   vector<float>   *trkPt;
   vector<float>   *trkPtError;
   vector<float>   *trkEta;
   vector<float>   *trkPhi;
   vector<char>    *trkCharge;
   vector<int>     *trkPDFId;
   vector<char>    *trkNHits;
   vector<char>    *trkNPixHits;
   vector<char>    *trkNLayers;
   vector<float>   *trkNormChi2;
   vector<bool>    *highPurity;
   vector<int>     *trkAssociatedVtxIndx;
   vector<int>     *trkAssociatedVtxQuality;
   vector<float>   *trkDzAssociatedVtx;
   vector<float>   *trkDzErrAssociatedVtx;
   vector<float>   *trkDxyAssociatedVtx;
   vector<float>   *trkDxyErrAssociatedVtx;
   vector<int>     *trkFirstVtxQuality;
   vector<float>   *trkDzFirstVtx;
   vector<float>   *trkDzErrFirstVtx;
   vector<float>   *trkDxyFirstVtx;
   vector<float>   *trkDxyErrFirstVtx;

   // List of branches
   TBranch        *b_nRun;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_xErrVtx;   //!
   TBranch        *b_yErrVtx;   //!
   TBranch        *b_zErrVtx;   //!
   TBranch        *b_chi2Vtx;   //!
   TBranch        *b_ndofVtx;   //!
   TBranch        *b_isFakeVtx;   //!
   TBranch        *b_nTracksVtx;   //!
   TBranch        *b_ptSumVtx;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkPDFId;   //!
   TBranch        *b_trkNHits;   //!
   TBranch        *b_trkNPixHits;   //!
   TBranch        *b_trkNLayers;   //!
   TBranch        *b_trkNormChi2;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_trkAssociatedVtxIndx;   //!
   TBranch        *b_trkAssociatedVtxQuality;   //!
   TBranch        *b_trkDzAssociatedVtx;   //!
   TBranch        *b_trkDzErrAssociatedVtx;   //!
   TBranch        *b_trkDxyAssociatedVtx;   //!
   TBranch        *b_trkDxyErrAssociatedVtx;   //!
   TBranch        *b_trkFirstVtxQuality;   //!
   TBranch        *b_trkDzFirstVtx;   //!
   TBranch        *b_trkDzErrFirstVtx;   //!
   TBranch        *b_trkDxyFirstVtx;   //!
   TBranch        *b_trkDxyErrFirstVtx;   //!

   PbPbTracks(TTree *tree=0);
   virtual ~PbPbTracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PbPbTracks_cxx
PbPbTracks::PbPbTracks(TTree *tree) : fChain(0) 
{

}

PbPbTracks::~PbPbTracks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PbPbTracks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PbPbTracks::LoadTree(Long64_t entry)
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

void PbPbTracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   xVtx = 0;
   yVtx = 0;
   zVtx = 0;
   xErrVtx = 0;
   yErrVtx = 0;
   zErrVtx = 0;
   chi2Vtx = 0;
   ndofVtx = 0;
   isFakeVtx = 0;
   nTracksVtx = 0;
   ptSumVtx = 0;
   trkPt = 0;
   trkPtError = 0;
   trkEta = 0;
   trkPhi = 0;
   trkCharge = 0;
   trkPDFId = 0;
   trkNHits = 0;
   trkNPixHits = 0;
   trkNLayers = 0;
   trkNormChi2 = 0;
   highPurity = 0;
   trkAssociatedVtxIndx = 0;
   trkAssociatedVtxQuality = 0;
   trkDzAssociatedVtx = 0;
   trkDzErrAssociatedVtx = 0;
   trkDxyAssociatedVtx = 0;
   trkDxyErrAssociatedVtx = 0;
   trkFirstVtxQuality = 0;
   trkDzFirstVtx = 0;
   trkDzErrFirstVtx = 0;
   trkDxyFirstVtx = 0;
   trkDxyErrFirstVtx = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nEv", &nEv, &b_nRun);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("xVtx", &xVtx, &b_xVtx);
   fChain->SetBranchAddress("yVtx", &yVtx, &b_yVtx);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   fChain->SetBranchAddress("xErrVtx", &xErrVtx, &b_xErrVtx);
   fChain->SetBranchAddress("yErrVtx", &yErrVtx, &b_yErrVtx);
   fChain->SetBranchAddress("zErrVtx", &zErrVtx, &b_zErrVtx);
   fChain->SetBranchAddress("chi2Vtx", &chi2Vtx, &b_chi2Vtx);
   fChain->SetBranchAddress("ndofVtx", &ndofVtx, &b_ndofVtx);
   fChain->SetBranchAddress("isFakeVtx", &isFakeVtx, &b_isFakeVtx);
   fChain->SetBranchAddress("nTracksVtx", &nTracksVtx, &b_nTracksVtx);
   fChain->SetBranchAddress("ptSumVtx", &ptSumVtx, &b_ptSumVtx);
   fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkPtError", &trkPtError, &b_trkPtError);
   fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", &trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkCharge", &trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("trkPDFId", &trkPDFId, &b_trkPDFId);
   fChain->SetBranchAddress("trkNHits", &trkNHits, &b_trkNHits);
   fChain->SetBranchAddress("trkNPixHits", &trkNPixHits, &b_trkNPixHits);
   fChain->SetBranchAddress("trkNLayers", &trkNLayers, &b_trkNLayers);
   fChain->SetBranchAddress("trkNormChi2", &trkNormChi2, &b_trkNormChi2);
   fChain->SetBranchAddress("highPurity", &highPurity, &b_highPurity);
   fChain->SetBranchAddress("trkAssociatedVtxIndx", &trkAssociatedVtxIndx, &b_trkAssociatedVtxIndx);
   fChain->SetBranchAddress("trkAssociatedVtxQuality", &trkAssociatedVtxQuality, &b_trkAssociatedVtxQuality);
   fChain->SetBranchAddress("trkDzAssociatedVtx", &trkDzAssociatedVtx, &b_trkDzAssociatedVtx);
   fChain->SetBranchAddress("trkDzErrAssociatedVtx", &trkDzErrAssociatedVtx, &b_trkDzErrAssociatedVtx);
   fChain->SetBranchAddress("trkDxyAssociatedVtx", &trkDxyAssociatedVtx, &b_trkDxyAssociatedVtx);
   fChain->SetBranchAddress("trkDxyErrAssociatedVtx", &trkDxyErrAssociatedVtx, &b_trkDxyErrAssociatedVtx);
   fChain->SetBranchAddress("trkFirstVtxQuality", &trkFirstVtxQuality, &b_trkFirstVtxQuality);
   fChain->SetBranchAddress("trkDzFirstVtx", &trkDzFirstVtx, &b_trkDzFirstVtx);
   fChain->SetBranchAddress("trkDzErrFirstVtx", &trkDzErrFirstVtx, &b_trkDzErrFirstVtx);
   fChain->SetBranchAddress("trkDxyFirstVtx", &trkDxyFirstVtx, &b_trkDxyFirstVtx);
   fChain->SetBranchAddress("trkDxyErrFirstVtx", &trkDxyErrFirstVtx, &b_trkDxyErrFirstVtx);
   Notify();
}

Bool_t PbPbTracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PbPbTracks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PbPbTracks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PbPbTracks_cxx
