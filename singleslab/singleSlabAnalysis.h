//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 18 11:06:32 2017 by ROOT version 5.34/34
// from TTree fev10/fev10
// found on file: cosmics_dif_1_1_commissioning.raw.root
//////////////////////////////////////////////////////////

#ifndef singleSlabAnalysis_h
#define singleSlabAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
//#include "TIter.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class singleSlabAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TString         datadir;

   // Declaration of leaf types
   Int_t           event;
   Int_t           thr;
   Int_t           acqNumber;
   Int_t           chipid[16];
   Int_t           nColumns[16];
   Int_t           bcid[16][15];
   Int_t           corrected_bcid[16][15];
   Int_t           badbcid[16][15];
   Int_t           nhits[16][15];
   Int_t           charge_lowGain[16][15][64];
   Int_t           charge_hiGain[16][15][64];
   Int_t           gain_hit_low[16][15][64];
   Int_t           gain_hit_high[16][15][64];

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_thr;   //!
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

   //singleSlabAnalysis(TTree *tree=0);
   singleSlabAnalysis(TString tree_s);
   singleSlabAnalysis(TList *f=0);

   virtual ~singleSlabAnalysis();
   // int     main(int argc, char* argv[2]);
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   // write file with masked channels
   virtual void     FindMasked(TString dif);
   // analysis of pedestal and writting of the file with pedestals per chi/channel/sca
   virtual void     PedestalAnalysis(TString dif,TString grid,TString map_filename);
   virtual void     PedestalAnalysis_gridpoints(TString dif,TString grid,TString map_filename);
   virtual void     PedestalAnalysis_bcid(TString dif,TString grid,TString map_filename);
   // read pedestals, maps, masked channels
   virtual void     ReadMap(TString filename);
   virtual void     ReadMasked(TString filename);
   virtual void     ReadPedestals(TString filename);
   //signal analysis: MIP fitt and signal/noise 
   virtual void     SignalAnalysis(TString dif, TString outputname, bool readpedestal, TString map_filename);
   virtual TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
   //   virtual Double_t langaufun(Double_t *x, Double_t *par);
   // bcid correlations of retriggers (good vs bad events)
   virtual void     BcidCorrelations(TString filename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

private :

  Float_t map_pointX[16][64];
  Float_t map_pointY[16][64];
  Int_t masked[16][64];
  std::vector<std::vector<std::vector<Double_t> > > ped_mean;
  std::vector<std::vector<std::vector<Double_t> > > ped_error;
  std::vector<std::vector<std::vector<Double_t> > > ped_width;



};

#endif

#ifdef singleSlabAnalysis_cxx
singleSlabAnalysis::singleSlabAnalysis(TString tree_s) : fChain(0) 
{

  
  TFile *f = new TFile(tree_s);
  TTree *tree = (TTree*)f->Get("fev10");
  //  tree->Print();
  Init(tree);
  
}

singleSlabAnalysis::singleSlabAnalysis(TList *f) : fChain(0) 
{
// if parameter tree is not specified (or zero), use a list of of files provided as input

  TIter next(f);
  TSystemFile *file;
  TString fname;
  while((file = (TSystemFile*)next())){
      fname = file->GetName();
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
      TTree *tree=0;
      f->GetObject("fev10",tree);
      Init(tree);
  }


}

singleSlabAnalysis::~singleSlabAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t singleSlabAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t singleSlabAnalysis::LoadTree(Long64_t entry)
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



void singleSlabAnalysis::Init(TTree *tree)
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
   fChain->SetBranchAddress("thr", &thr, &b_thr);
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

Bool_t singleSlabAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
  return kTRUE;
}

void singleSlabAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t singleSlabAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef singleSlabAnalysis_cxx
