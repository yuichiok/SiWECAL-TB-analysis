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


void createDummyMIPfiles() 
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

  
  TString filename[9];
  filename[0]="dif_1_1_1";
  filename[1]="dif_1_1_2";
  filename[2]="dif_1_1_3";
  filename[3]="dif_1_1_4";
  filename[4]="dif_1_1_5";
  filename[5]="SLB_2";
  filename[6]="SLB_1";
  filename[7]="SLB_3";
  filename[8]="SLB_0";

  double mip_value[9];
  mip_value[0]=650./500.;
  mip_value[1]=650./500.;
  mip_value[2]=650./500.;
  mip_value[3]=320./500.;
  mip_value[4]=650./500.;
  mip_value[5]=500./500.;
  mip_value[6]=500./500.;
  mip_value[7]=500./500.;
  mip_value[8]=500./500.;

  for(int i=0; i<9; i++) {
    ofstream fout_mip("MIP_"+filename[i]+"_dummy.txt",ios::out);
    fout_mip<<"#mip results "<<filename[i]<<endl;
    fout_mip<<"#chip channel mpv empv widthmpv chi2ndf nentries"<<endl;

    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
	if(j<4) fout_mip<<j<<" "<<k<<" "<<mip_value[i]<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;
	else {
	  if(i<5)  fout_mip<<j<<" "<<k<<" "<<mip_value[i]<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;
	  else fout_mip<<j<<" "<<k<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<endl;
	}
      }
    }
  }
}
