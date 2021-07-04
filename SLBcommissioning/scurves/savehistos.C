//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include <iostream>
#include <fstream>

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;

Float_t value[2][15][300][64][100]; //daughter, slab, asic, chn, DAC
Float_t error_value[2][15][300][64][100]; //dauchter, slab, asic, chn, DAC
Float_t x[100];
Int_t nsteps;
Int_t nslabs[2];
int idslab[2][15];

void ReadScurves(TString filename)
  
{
  /*
 === RATE Vs THRESHOLD HISTO SAVED WITH ILC SL-SOFTWARE VERSION: V2.10  == DATE OF SCAN: UnixTime = 1576752746.986 date = 2019.12.19 ===
=== SL_COREMODULE_INTERFACE SerNum 2.41A FPGA Version V2.1.11 == NB OF CORE DAUGHTERS: 1 ===
====== CORE DAUGHTER 0 == FPGA Version V1.2.6 == NB OF CONNECTED SLABs 1 ======
========= DAUGHTER 0 == SLAB 0 == SL BOARD ADD 3 == FPGA Version V1.3.4 == NB OF CONNECTED ASUs 1 =========
=== Hit Rate Vs Threshold Scan,for All Channels modulo 64 at once, AcqWindow 2 ms, MinThreshold = 170, MaxThreshold= 240, Nb Of Steps: 7, Nb Of Cycles per Step= 5 === 
  */


 
 
  cout<<"Reading scurve file: "<<filename<<endl;
  std::ifstream reading_file(filename);
  if(!reading_file){
    std::cout<<"NO FILE!!"<<std::endl;
  }

  for(int k=0; k<100; k++) {
    x[k] = 0;
    for(int icore=0; icore<2; icore++) {
      for(int islab=0; islab<15; islab++) {
	for(int i=0; i<300; i++) {
	  for(int j=0; j<64; j++) {
	    value[icore][islab][i][j][k] = 0.;
	    error_value[icore][islab][i][j][k] = 0.;
	  }
	}
     }
    }
  }
  
 
  TString tmpst;
  //  int nslabs[2];
  int ndaughter=0;
  int nasu[2][15];

  for(int i=0; i<2; i++) {
    nslabs[i]=0;
    for(int j=0; j<15; j++) {
      nasu[i][j]=0;
      idslab[i][j]=0;
    }
  }
      
  TString nsteps_tmp;

  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst>> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> ndaughter >> tmpst;
  for(int k=0; k<ndaughter; k++) {
    reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> nslabs[k]>> tmpst ;
    for(int i=0; i<nslabs[k]; i++) {
      reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> idslab[k][i] >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst>> tmpst >> tmpst >> tmpst >> tmpst >> nasu[k][i] >> tmpst ;
    }
  }
   
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst>> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >>  tmpst >> nsteps_tmp >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;

  std::cout<< "-----" << ndaughter<< " " <<nslabs[0] << " " << nslabs[1] << " " << std::endl;
  /*  for(int i=0; i<2; i++) {
    for(int j=0; j<15; j++) {
      std::cout<< nasu[i][j]<< " "<< idslab[i][j] << std::endl;
    }
    }*/
  nsteps_tmp.Remove(nsteps_tmp.Length()-1);
  nsteps=nsteps_tmp.Atoi();

  /*  == X AXIS : Threshold[7 values] ==
170 180 190 200 210 220 230 
== Y AXIS for CORE 0 SLAB 0 ASU 0 SKIROC 0 CH 0 ==*/

  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  for(int i=0; i<nsteps; i++) {
    reading_file >> x[i] ;
  }

  while(reading_file){
    int core =-1, asu=-1, asic=-1, chn=-1, slab=-1;
    reading_file >> tmpst >> tmpst  >> tmpst >> tmpst >> tmpst >> core >> tmpst >> slab >> tmpst >> asu >> tmpst >> asic >> tmpst >> chn >> tmpst;
    asic=asic+16*asu;
    if(asu >-1 && asic>-1 && chn >-1 && slab >-1) {
      for(int i=0; i<nsteps; i++) {
	reading_file >> value[core][slab][asic][chn][i];
	if(value[core][slab][asic][chn][i]>0) error_value[core][slab][asic][chn][i]=sqrt(value[core][slab][asic][chn][i]);
      }
    }
  }
   

  

}

void savehistos(TString filename = ""){//TString filename, int slabadd){
  

  TString filename_input = "../"+filename+"/RateVsThresholdScan_"+filename+"_SLBoard.txt";
  TGraphErrors *scurve[15][16][64];

  cout<<filename_input<<endl;
 ReadScurves(filename_input);

 for(int i=0; i<nslabs[0]; i++) {
   for(int iasic=0; iasic<16; iasic++) {
     for(int ichn=0; ichn<64; ichn++) {
       scurve[i][iasic][ichn]= new TGraphErrors(nsteps,x,value[0][i][iasic][ichn],0,error_value[0][i][iasic][ichn]);
     }
   }
 }


 TCanvas *canvas[15];
  TCanvas *canvas_asic[15][16];

  TFile *file_summary = new TFile("histos/scurves_"+filename+".root" , "RECREATE");
  file_summary->cd();

  for(int i=0; i<nslabs[0]; i++) {
    canvas[i]= new TCanvas(TString::Format("canvas_slboard%i",idslab[0][i]),TString::Format("canvas_slboard%i",idslab[0][i]),1600,1000);
    canvas[i]->Divide(4,4);
    
    //  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);

    for(int iasic=0; iasic<16; iasic++) {
      canvas[i]->cd(1+iasic);
      for(int ichn=0; ichn<64; ichn++) {
	if(ichn>0) scurve[i][iasic][ichn]->Draw("l");
	else {
	  scurve[i][iasic][ichn]->GetXaxis()->SetTitle("DAC");
	  scurve[i][iasic][ichn]->GetYaxis()->SetTitle("Nhits");
	  scurve[i][iasic][ichn]->SetTitle(TString::Format("ASIC %i",iasic));
	  //	  if(iasic==0) leg->AddEntry(scurve[i][iasic][ichn],"HV=100V","l");
	  scurve[i][iasic][ichn]->Draw("al");
	}
	scurve[i][iasic][ichn]->SetName(TString::Format("scurve_slboard%i_chip%i_chn%i",idslab[0][i],iasic,ichn));
	file_summary->cd();
	scurve[i][iasic][ichn]->Write();
	
      }
      //  leg->Draw();
    }
    
  }
  
  
}





