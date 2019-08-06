//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug  1 12:02:43 2019 by ROOT version 6.18/00
// from TTree layers/layers
// found on file: MergeLayers.root
//////////////////////////////////////////////////////////

#ifndef Layers_compar_h
#define Layers_compar_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Layers_compar {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           acqNumber;
   Int_t           chipid[9][16];
   Int_t           nColumns[9][16];
   Int_t           bcid[9][16][15];
   Int_t           corrected_bcid[9][16][15];
   Int_t           badbcid[9][16][15];
   Int_t           nhits[9][16][15];
   Int_t           charge_lowGain[9][16][15][64];
   Int_t           charge_hiGain[9][16][15][64];
   Int_t           gain_hit_low[9][16][15][64];
   Int_t           gain_hit_high[9][16][15][64];

   // List of branches
   TBranch        *b_acqNumber;   //!
   TBranch        *b_chipid;   //!
   TBranch        *b_nColumns;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_corrected_bcid;   //!
   TBranch        *b_badbcid;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_lowGain;   //!
   TBranch        *b_highGain;   //!
   TBranch        *b_gain_hit_low;   //!
   TBranch        *b_gain_hit_high;   //!

   Layers_compar(TTree *tree=0);
   virtual ~Layers_compar();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Layers_compar_cxx
Layers_compar::Layers_compar(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../../EvBuilt/run_32015_merge.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../../EvBuilt/run_32015_merge.root");
      }
      f->GetObject("ecal_raw",tree);

   }
   Init(tree);
}

Layers_compar::~Layers_compar()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Layers_compar::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Layers_compar::LoadTree(Long64_t entry)
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

void Layers_compar::Init(TTree *tree)
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

   fChain->SetBranchAddress("acqNumber", &acqNumber, &b_acqNumber);
   fChain->SetBranchAddress("chipid", chipid, &b_chipid);
   fChain->SetBranchAddress("nColumns", nColumns, &b_nColumns);
   fChain->SetBranchAddress("bcid", bcid, &b_bcid);
   fChain->SetBranchAddress("corrected_bcid", corrected_bcid, &b_corrected_bcid);
   fChain->SetBranchAddress("badbcid", badbcid, &b_badbcid);
   fChain->SetBranchAddress("nhits", nhits, &b_nhits);
   fChain->SetBranchAddress("charge_lowGain", charge_lowGain, &b_lowGain);
   fChain->SetBranchAddress("charge_hiGain", charge_hiGain, &b_highGain);
   fChain->SetBranchAddress("gain_hit_low", gain_hit_low, &b_gain_hit_low);
   fChain->SetBranchAddress("gain_hit_high", gain_hit_high, &b_gain_hit_high);
   Notify();
}

Bool_t Layers_compar::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Layers_compar::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Layers_compar::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Layers_compar_cxx
