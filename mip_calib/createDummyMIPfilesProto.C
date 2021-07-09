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


void createDummyMIPfilesProto() 
{

  //beam line
  //9 layers, w.r.t. beam point of view
  //dif_1_1_1 , FEV13Jp P1, 4x650um wafers
  //dif_1_1_2 , FEV13Jp P2, 4x650um wafers
  //dif_1_1_3 , FEV13Jp P3, 4x320um wafers
  //dif_1_1_4 , FEV13Jp K1, 4x650um wafers
  //dif_1_1_5 , FEV13Jp K1, 4x650um wafers
  //SLB_2 = COB_a, 500um wafer on chips 0-3
  //SLB_1 = FEV12, 500um wafer on chips 0-3
  //SLB_3 = FEV12, 500um wafer on chips 0-3
  //SLB_0 = COB_c, 500um wafer on chips 0-3

  double mip_value[15]={1};

  ofstream fout_mip("MIP_PROTO15_dummy.txt",ios::out);
  fout_mip<<"#mip results PROTO15"<<endl;
  fout_mip<<"#layer chip channel mpv empv widthmpv chi2ndf nentries"<<endl;

  for(int i=0; i<15; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	fout_mip<<i<<" "<<j<<" "<<k<<" "<<mip_value[i]<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;
      }
    }
  }

}

