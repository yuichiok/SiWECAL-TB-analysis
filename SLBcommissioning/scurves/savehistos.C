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
#include "../conf_struct.h"

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
  
 
  // TString tmpst;
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

  std::string line;
  bool isCore = false;

  // check USB or Core
  getline(reading_file, line);
  getline(reading_file, line);


  std::size_t found_core = line.find("SL_COREMODULE_INTERFACE");
  std::size_t found_usb  = line.find("SL_DIRECT_INTERFACE");

  if(found_core!=std::string::npos){
    isCore = true;
  }else if(found_usb!=std::string::npos){
    isCore = false;
  }else{
    exit (EXIT_FAILURE);
  }



  std::cout << line << std::endl;

  if(!isCore){
    ndaughter = 1;
    nslabs[0] = 1;
    idslab[0][0] = 0;

    search_string_nocolon(line,"NB OF CONNECTED ASUS",nasu[0][0]);

  }else{

    search_string(line,"NB OF CORE DAUGHTERS:",ndaughter);

    for(int k=0; k<ndaughter; k++) {

      getline(reading_file,line);
      search_string_nocolon(line,"NB OF CONNECTED SLABs",nslabs[k]);

      for(int i=0; i<nslabs[k]; i++) {

        getline(reading_file,line);
        search_string_nocolon(line,"SL BOARD ADD",idslab[k][i]);
        search_string_nocolon(line,"NB OF CONNECTED ASUs",nasu[k][i]);

      }

    }

  }

  getline(reading_file,line);

  std::string nsteps_tmp;
  search_string(line,"Nb Of Steps:",nsteps_tmp);
  nsteps_tmp = nsteps_tmp.substr(0,nsteps_tmp.size()-1); // Nb Of Steps: 15, (rm ,)
  nsteps=std::stoi(nsteps_tmp);


  /*  
  == X AXIS : Threshold[7 values] ==
  170 180 190 200 210 220 230 
  */

  getline(reading_file,line);
  getline(reading_file,line);

  std::string space_delimiter = " ";
  size_t pos = 0;
  int i = 0;
  while ((pos = line.find(space_delimiter)) != string::npos) {
      if(i>nsteps) break;
      x[i] = std::stoi(line.substr(0, pos));
      line.erase(0, pos + space_delimiter.length());
      i++;
  }



  while( getline(reading_file,line) ) {
    // == Y AXIS for CORE 0 SLAB 0 ASU 0 SKIROC 0 CH 0 ==
  
    int core =-1, asu=-1, asic=-1, chn=-1, slab=-1;
    search_string_nocolon(line,"CORE",core);
    search_string_nocolon(line,"SLAB",slab);
    search_string_nocolon(line,"ASU",asu);
    search_string_nocolon(line,"SKIROC",asic);
    search_string_nocolon(line,"CH",chn);
    asic = asic + 16*asu;
  
    getline(reading_file,line);
    // 103 119 118 121 125 123 122 75 11 3 0 0 0 0 0
  
    size_t pos2 = 0;
    int j = 0;
    while ((pos2 = line.find(space_delimiter)) != string::npos) {
        if(j>nsteps) break;
        value[core][slab][asic][chn][j] = std::stoi(line.substr(0, pos2));  
        if(value[core][slab][asic][chn][j]>0) error_value[core][slab][asic][chn][j]=sqrt(value[core][slab][asic][chn][j]);
        line.erase(0, pos2 + space_delimiter.length());
        j++;
    }
  
  } // next == Y AXIS ...
  

}

void savehistos(TString filename = "", TString date=""){//TString filename, int slabadd){
  
  if(date=="") date=filename;
  TString filename_input = "../"+date+"/RateVsThresholdScan_"+filename+"_SLBoard.txt";
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
    canvas[i]= new TCanvas(TString::Format("canvas_slbAdd%i",idslab[0][i]),TString::Format("canvas_slbAdd%i",idslab[0][i]),1600,1000);
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
	scurve[i][iasic][ichn]->SetName(TString::Format("scurve_slbAdd%i_chip%i_chn%i",idslab[0][i],iasic,ichn));
	file_summary->cd();
	scurve[i][iasic][ichn]->Write();
	
      }
      //  leg->Draw();
    }
    
  }
  
  
}





