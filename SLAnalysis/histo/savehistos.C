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

Float_t value[2][15][300][64][100]; //slab, asic, chn, DAC
Float_t error_value[2][15][300][64][100]; //slab, asic, chn, DAC
Float_t x[100];
Int_t nsteps;

void ReadScurves(TString filename, int slabadd, int coreadd)
  
{
  /*
 === RATE Vs THRESHOLD HISTO SAVED WITH ILC SL-SOFTWARE VERSION: V2.10  == DATE OF SCAN: UnixTime = 1576752746.986 date = 2019.12.19 ===
=== SL_COREMODULE_INTERFACE SerNum 2.41A FPGA Version V2.1.11 == NB OF CORE DAUGHTERS: 1 ===
====== CORE DAUGHTER 0 == FPGA Version V1.2.6 == NB OF CONNECTED SLABs 1 ======
========= DAUGHTER 0 == SLAB 0 == SL BOARD ADD 3 == FPGA Version V1.3.4 == NB OF CONNECTED ASUs 1 =========
=== Hit Rate Vs Threshold Scan,for All Channels modulo 64 at once, AcqWindow 2 ms, MinThreshold = 170, MaxThreshold= 240, Nb Of Steps: 7, Nb Of Cycles per Step= 5 === 
  */

  filename="/home/calice/TB2020/commissioning/data_slboard/Histos/"+filename;
  
  std::ifstream reading_file(filename);
  if(!reading_file){
    std::cout<<" dameyo - damedame"<<std::endl;
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
  int nslabs[2];
  int ndaughter=0;
  int nasu[2][15];
  int idslab[2][15];

  for(int i=0; i<2; i++) {
    for(int j=0; j<15; j++) {
      nslabs[i]=0;
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

  std::cout<< ndaughter<< " " <<nslabs[0] << " " << nslabs[1] << " " << std::endl;
  for(int i=0; i<2; i++) {
    for(int j=0; j<15; j++) {
      std::cout<< nasu[i][j]<< " "<< idslab[i][j] << std::endl;
    }
  }
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
    if(core == coreadd && asu >-1 && asic>-1 && chn >-1 && slab >-1 && idslab[core][slab]==slabadd) {
      for(int i=0; i<nsteps; i++) {
	reading_file >> value[core][slab][asic][chn][i];
	if(value[core][slab][asic][chn][i]>0) error_value[core][slab][asic][chn][i]=sqrt(value[core][slab][asic][chn][i]);
      }
    }
  }
   

  

}

void savehistos(){//TString filename, int slabadd){
  

  TString filename = "RateVsThresholdScan_01022020_SLBoard_slab14_HV100.txt";
  int slabadd=3;
  ReadScurves(filename,slabadd,0);

  TGraphErrors *scurve_1[16][64];
  for(int iasic=0; iasic<16; iasic++) {
      for(int ichn=0; ichn<64; ichn++) {
	scurve_1[iasic][ichn]= new TGraphErrors(nsteps,x,value[0][0][iasic][ichn],0,error_value[0][0][iasic][ichn]);
      }
  }

  filename = "RateVsThresholdScan_01022020_SLBoard_slab14_HV180.txt";
  ReadScurves(filename,slabadd,0);

  TGraphErrors *scurve_2[16][64];
  for(int iasic=0; iasic<16; iasic++) {
      for(int ichn=0; ichn<64; ichn++) {
	scurve_2[iasic][ichn]= new TGraphErrors(nsteps,x,value[0][0][iasic][ichn],0,error_value[0][0][iasic][ichn]);
      }
  }

  TCanvas *canvas = new TCanvas("canvas","canvas",1600,1000);
  canvas->Divide(4,4);

  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);

  for(int iasic=0; iasic<16; iasic++) {
    canvas->cd(1+iasic);
    for(int ichn=0; ichn<64; ichn++) {
      if(ichn>0) scurve_1[iasic][ichn]->Draw("l");
      else {
	scurve_1[iasic][ichn]->GetXaxis()->SetTitle("DAC");
	scurve_1[iasic][ichn]->GetYaxis()->SetTitle("Nhits");
	scurve_1[iasic][ichn]->SetTitle(TString::Format("ASIC %i",iasic));
	if(iasic==0) leg->AddEntry(scurve_1[iasic][ichn],"HV=100V","l");
	scurve_1[iasic][ichn]->Draw("al");
      }
      scurve_2[iasic][ichn]->SetLineColor(2);
      scurve_2[iasic][ichn]->Draw("l");
      if(iasic==0 && ichn==0) leg->AddEntry(scurve_2[iasic][ichn],"HV=180V","l");
      
    }
    leg->Draw();
  }


  TCanvas *canvas_asic[16];

  for(int iasic=0; iasic<16; iasic++) {
    
    canvas_asic[iasic]= new TCanvas(TString::Format("asic_%i",iasic),TString::Format("asic_%i",iasic),1600,1000);
    canvas_asic[iasic]->Divide(8,8);
    
    TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
    
    for(int ichn=0; ichn<64; ichn++) {
      canvas_asic[iasic]->cd(1+ichn);
      scurve_1[iasic][ichn]->GetXaxis()->SetTitle("DAC");
      scurve_1[iasic][ichn]->GetYaxis()->SetTitle("Nhits");
      scurve_1[iasic][ichn]->SetTitle(TString::Format("ASIC %i",iasic));
      if(ichn==0) leg2->AddEntry(scurve_1[iasic][ichn],"HV=100V","l");
      scurve_1[iasic][ichn]->Draw("al");
      
      scurve_2[iasic][ichn]->SetLineColor(2);
      scurve_2[iasic][ichn]->Draw("l");
      if(ichn==0) leg2->AddEntry(scurve_2[iasic][ichn],"HV=180V","l");
      leg2->Draw();
    }
  }
  
  
}





