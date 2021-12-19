//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 18 11:06:32 2017 by ROOT version 5.34/34
// from TTree slboard/slboard
// found on file: cosmics_dif_1_1_commissioning.raw.root
//////////////////////////////////////////////////////////

#ifndef DecodedSLBAnalysis_h
#define DecodedSLBAnalysis_h
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
#include "../include/utils.h"

//#include "TIter.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class DecodedSLBAnalysis {
public :

   
 //Tree - MyClass
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TString         datadir;


  DecodedSLBAnalysis(TString tree_s, TString treename);
  //  DecodedSLBAnalysis(TString tree_s);
  //DecodedSLBAnalysis(TList *f=0);
  virtual ~DecodedSLBAnalysis();
  // int     main(int argc, char* argv[2]);
  virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);

  //Proto Analysis
  virtual void     HitMonitoring(TString, int, int);
  virtual int     NSlabsAnalysis(TString, int, int);
  

   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

  
private :


};

#endif

#ifdef DecodedSLBAnalysis_cxx
DecodedSLBAnalysis::DecodedSLBAnalysis(TString tree_s, TString treename="siwecaldecoded") : fChain(0)
{
  cout<<tree_s<<endl;
  TFile *f = new TFile(tree_s);
  TTree *tree = (TTree*)f->Get(treename);
  //  tree->Print();                                                                                                                                                                                       
  Init(tree);

}
/*
DecodedSLBAnalysis::DecodedSLBAnalysis(TList *f) : fChain(0) 
{
  // if parameter tree is not specified (or zero), use a list of of files provided as input

  TIter next(f);
  TSystemFile *file;
  TString fname;
  while((file = (TSystemFile*)next())){
    fname = file->GetName();
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
    TTree *tree=0;
    f->GetObject("siwecaldecoded",tree);
    Init(tree);
  }

}
*/
DecodedSLBAnalysis::~DecodedSLBAnalysis()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


Int_t DecodedSLBAnalysis::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t DecodedSLBAnalysis::LoadTree(Long64_t entry)
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



void DecodedSLBAnalysis::Init(TTree *tree)
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
   SetBranchAddressFunction(fChain);
   Notify();
}

Bool_t DecodedSLBAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DecodedSLBAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DecodedSLBAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DecodedSLBAnalysis_cxx
