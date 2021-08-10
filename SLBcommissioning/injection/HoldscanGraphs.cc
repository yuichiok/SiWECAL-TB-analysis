//# Copyright 2020  Adrián Irles

#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"
#include "TGraphErrors.h"
#include "TLegend.h"

void HoldscanGraphs(TString filename_in="../../converter_SLB/convertedfiles/", TString run="07282021_dac1.2V_small", int nslboards=15){
  
  cout<<" Holdscan file: "<<filename_in<<endl;

  float y[15][16][64][30];
  float x[30];
  float ex[30];
  float ey[15][16][64][30];

  for(int l=0; l<30; l++) {
    x[l]=0;
    ex[l]=0;
    for(int i=0; i<15; i++) {
      for(int j=0; j<16; j++) {
  	for(int k=0; k<64; k++) {
  	  y[i][j][k][l]=0;
  	  ey[i][j][k][l]=0;
  	}
      }
    }
  }


  
  int board_id[15];
 
  for(int i=0; i<15; i++) {
    board_id[i]=-1;
  }
  
  int count=0;
  for(int i=20; i<170; i=i+20) {
    cout<<" i: "<<i<<endl;

    TString filename=filename_in+TString::Format("%s/holdscan_hold%i.root",run.Data(),i);
    cout<<" Holdscan file: "<<filename_in<<endl;
    DecodedSLBAnalysis ss(filename);
    std::vector<std::array<float,6>> holdvalues =ss.HoldscanAnalysis(-1);
    
    x[count]=i;
    
    for(int i=0; i<holdvalues.size(); i++) {
      if(board_id[int(holdvalues.at(i)[0])]==-1) board_id[int(holdvalues.at(i)[0])]=int(holdvalues.at(i)[1]);

      y[int(holdvalues.at(i)[0])][int(holdvalues.at(i)[2])][int(holdvalues.at(i)[3])][count]=holdvalues.at(i)[4];
      ey[int(holdvalues.at(i)[0])][int(holdvalues.at(i)[2])][int(holdvalues.at(i)[3])][count]=holdvalues.at(i)[5];
    }
    count++;
    
  }


  TGraphErrors* holdscan[15][16][64];

  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
  	holdscan[i][j][k]= new TGraphErrors(count,x,y[i][j][k],ex,ey[i][j][k]);
      }
    }
  }

  TFile *file_summary = new TFile("results/graphs_holdscan_"+run+".root" , "RECREATE");
  file_summary->cd();

  TH1F *h1 = new TH1F("layer_slboard_relation","layer_slboard_relation",nslboards,-0.5,nslboards-0.5);
  for(int i=0; i<nslboards; i++) {
    h1->Fill(i,board_id[i]);
  }
  h1->Write();
  
  for(int i=0; i<nslboards; i++) {
    for(int j=0; j<16; j++) {
      for(int k=0; k<64; k++) {
        holdscan[i][j][k]->SetName(TString::Format("holdscan_layer%i_slboard%i_chip%i_chn%i",i,board_id[i],j,k));
        holdscan[i][j][k]->Write();
      }
    }
  }
  
  file_summary->Close();

}
  

