//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 18 11:06:32 2017 by ROOT version 5.34/34
// from TTree fev10/fev10
// found on file: cosmics_dif_1_1_commissioning.raw.root
//////////////////////////////////////////////////////////

#ifndef example_h
#define example_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
//#include "TIter.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class example {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t NSLABS=9;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           spill;
   Int_t           bcid;
   Int_t           prev_bcid;
   Int_t           next_bcid;
   Int_t           nhit_slab;
   Int_t           nhit_chip;
   Int_t           nhit_chan;
   Float_t         sum_hg;
   Float_t         sum_energy;
   Int_t           hit_slab[8000];   //[nhit_chan]
   Int_t           hit_chip[8000];   //[nhit_chan]
   Int_t           hit_chan[8000];   //[nhit_chan]
   Int_t           hit_sca[8000];   //[nhit_chan]
   Float_t         hit_x[8000];   //[nhit_chan]
   Float_t         hit_y[8000];   //[nhit_chan]
   Float_t         hit_z[8000];   //[nhit_chan]
   Float_t         hit_x0[8000];   //[nhit_chan]
   Float_t         hit_hg[8000];   //[nhit_chan]
   Float_t         hit_lg[8000];   //[nhit_chan]
   Float_t         hit_energy[8000];   //[nhit_chan]
   Int_t           hit_isHit[8000];   //[nhit_chan]
   Int_t           hit_isMasked[8000];   //[nhit_chan]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_spill;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_prev_bcid;   //!
   TBranch        *b_next_bcid;   //!
   TBranch        *b_nhit_slab;   //!
   TBranch        *b_nhit_chip;   //!
   TBranch        *b_nhit_chan;   //!
   TBranch        *b_sum_hg;   //!
   TBranch        *b_sum_energy;   //!
   TBranch        *b_hit_slab;   //!
   TBranch        *b_hit_chip;   //!
   TBranch        *b_hit_chan;   //!
   TBranch        *b_hit_sca;   //!
   TBranch        *b_hit_x;   //!
   TBranch        *b_hit_y;   //!
   TBranch        *b_hit_z;   //!
   TBranch        *b_hit_x0;   //!
   TBranch        *b_hit_hg;   //!
   TBranch        *b_hit_lg;   //!
   TBranch        *b_hit_energy;   //!
   TBranch        *b_hit_isHit;   //!
   TBranch        *b_hit_isMasked;   //!

   //example(TTree *tree=0);
   example(TString tree_s);
   example(TList *f=0);

   virtual ~example();
   // int     main(int argc, char* argv[2]);
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);

   //Service functions
   virtual bool     IsHit(int ihit);
   virtual bool     IsPedestal(int ihit);
   virtual bool     TrackBasicSelection(int);

   //AnalysisFunctions
   virtual void SimpleDistributionsTrack(TString );
   virtual void SimpleEvDisplayTrack(TString );
   
   virtual TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);

   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

private :

   

};

#endif

#ifdef example_cxx
example::example(TString tree_s) : fChain(0) 
{

  
  TFile *f = new TFile(tree_s);
  TTree *tree = (TTree*)f->Get("ecal");
  //  tree->Print();
  Init(tree);
  
}

example::example(TList *f) : fChain(0) 
{
// if parameter tree is not specified (or zero), use a list of of files provided as input

  TIter next(f);
  TSystemFile *file;
  TString fname;
  while((file = (TSystemFile*)next())){
      fname = file->GetName();
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
      TTree *tree=0;
      f->GetObject("ecal",tree);
      Init(tree);
  }


}

example::~example()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t example::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t example::LoadTree(Long64_t entry)
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




void example::Init(TTree *tree)
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("spill", &spill, &b_spill);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("prev_bcid", &prev_bcid, &b_prev_bcid);
   fChain->SetBranchAddress("next_bcid", &next_bcid, &b_next_bcid);
   fChain->SetBranchAddress("nhit_slab", &nhit_slab, &b_nhit_slab);
   fChain->SetBranchAddress("nhit_chip", &nhit_chip, &b_nhit_chip);
   fChain->SetBranchAddress("nhit_chan", &nhit_chan, &b_nhit_chan);
   fChain->SetBranchAddress("sum_hg", &sum_hg, &b_sum_hg);
   fChain->SetBranchAddress("sum_energy", &sum_energy, &b_sum_energy);
   fChain->SetBranchAddress("hit_slab", hit_slab, &b_hit_slab);
   fChain->SetBranchAddress("hit_chip", hit_chip, &b_hit_chip);
   fChain->SetBranchAddress("hit_chan", hit_chan, &b_hit_chan);
   fChain->SetBranchAddress("hit_sca", hit_sca, &b_hit_sca);
   fChain->SetBranchAddress("hit_x", hit_x, &b_hit_x);
   fChain->SetBranchAddress("hit_y", hit_y, &b_hit_y);
   fChain->SetBranchAddress("hit_z", hit_z, &b_hit_z);
   fChain->SetBranchAddress("hit_x0", hit_x0, &b_hit_x0);
   fChain->SetBranchAddress("hit_hg", hit_hg, &b_hit_hg);
   fChain->SetBranchAddress("hit_lg", hit_lg, &b_hit_lg);
   fChain->SetBranchAddress("hit_energy", hit_energy, &b_hit_energy);
   fChain->SetBranchAddress("hit_isHit", hit_isHit, &b_hit_isHit);
   fChain->SetBranchAddress("hit_isMasked", hit_isMasked, &b_hit_isMasked);
   Notify();
}

Bool_t example::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void example::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t example::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef example_cxx
