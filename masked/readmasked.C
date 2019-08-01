#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
#include "TString.h"

void readmasked() 
{

  TString filename="masked_SLB_2.txt";//test.txt";
  std::ifstream reading_file(filename);
  cout<<filename<<endl;
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }

  int total_mask[16];
  for(int i=0; i<16; i++) {
    total_mask[i]=0;
  }
  Int_t tmp_chip = 0,tmp_channel = 0;
  Int_t tmp_mask = 0;

  //  #list of masked channels, per slboard: _SLB_0
  // #chip channel mask (0=not masked, 1=masked)
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_channel >> tmp_mask ;
    if(tmp_mask==1) total_mask[tmp_chip] ++ ;
  }

  for(int i=0; i<16; i++) {
    std::cout<< "chip="<<i<< " total mask=" <<100.*total_mask[i]/64.<<std::endl;
  }

  
}
