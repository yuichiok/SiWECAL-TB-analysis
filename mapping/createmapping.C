#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "TLatex.h"
#include "../include/utils.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;

void createmapping(TString run="PedestalMIP_3GeVMIPscan", TString gain="low"){
  TString name[2]={"fev10","fev11_cob_rotate"};

  for(int iname=0; iname<2; iname++) {

    for(int iflipx=0; iflipx<2; iflipx++) {
      for(int iflipy=0; iflipy<2; iflipy++) {

	float flipx=1;
	if(iflipx==1) flipx=-1;
    	float flipy=1;
	if(iflipy==1) flipy=-1;
	
  // TString mapstring="../../mapping/"+name[iname]+"_chip_channel_x_y_mapping.txt";
  TString mapstring=name[iname]+"_chip_channel_x_y_mapping.txt";
  ReadMap(mapstring,0);
  TH2F* map =new TH2F("map","map; chip; chn",16,-0.5,15.5,64,-0.5,63.5);
  TH2F* map_chip =new TH2F("map_chip","map; chip; chn",16,-0.5,15.5,1,-0.5,63.5);
  TH2F* mapxy =new TH2F("mapxy","map-xy; x; y",32,-90,90,32,-90,90);
  TH2F* mapIJ =new TH2F("mapIJ","map-IJ; I; J",32,-0.5,31.5,32,-0.5,31.5);
  TH2F* mapxy_chip =new TH2F("mapxy_chip","map-xy; x; y",32,-90,90,32,-90,90);
  TH2F* mapIJ_chip =new TH2F("mapIJ_chip","map-IJ; I; J",32,-0.5,31.5,32,-0.5,31.5);

  for(int i=0;i<16;i++){
    map_chip->Fill(i,0.,i+1);
    for(int j=0; j<64; j++) {
      mapxy_chip->Fill(-flipx*map_pointX[0][i][j],-flipy*map_pointY[0][i][j],i+1);
      mapxy->Fill(-flipx*map_pointX[0][i][j],-flipy*map_pointY[0][i][j],j);
      map->Fill(i,j);
    }
  }

  ofstream fout_ped(TString::Format("IJmapping/IJmapping_type_%s_flipx%i_flipy%i.txt",name[iname].Data(),iflipx,iflipy).Data(),ios::out);
  fout_ped<<"### IJmapping " <<endl;
  fout_ped<<"## TYPE: "<<name[iname]<< " FLIPX: "<<iflipx<< " FLIPY: "<<iflipy <<endl;
  fout_ped<<"#chip channel I J"<<endl;

  double mappingarray_i[32][64];
  double mappingarray_j[32][64];

  for(int i=0; i<32; i++) {
    for(int j=0; j<32; j++) {
      mapIJ->Fill(i,j,mapxy->GetBinContent(i+1,j+1));
      //if(mapIJ_chip->GetBinContent(i+1,j+1)==0)
      mapIJ_chip->Fill(i,j,mapxy_chip->GetBinContent(i+1,j+1));
      mappingarray_i[int(mapxy_chip->GetBinContent(i+1,j+1)-1)][int(mapxy->GetBinContent(i+1,j+1))]=i;
      mappingarray_j[int(mapxy_chip->GetBinContent(i+1,j+1)-1)][int(mapxy->GetBinContent(i+1,j+1))]=j;

    }
  }
  
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      fout_ped<<i<<" "<<j<<" "<<mappingarray_i[i][j]<<" "<<mappingarray_j[i][j]<<endl;
    }
  }
  


  gStyle->SetOptStat(0);
  TCanvas* canvas2= new TCanvas(TString::Format("IJmapping/IJmapping_type_%s_flipx%i_flipy%i.txt",name[iname].Data(),iflipx,iflipy),TString::Format("IJmapping_type_%s_flipx%i_flipy%i.txt",name[iname].Data(),iflipx,iflipy),1800,400);   
  canvas2->Divide(3,1);
  canvas2->cd(1);
  map_chip->Draw("colz");
  // map->Draw("textsame");

  canvas2->cd(2);
  mapxy_chip->Draw("colz");
  //  mapxy->SetMarkerColor(0);
  mapxy->Draw("textsame");
  canvas2->cd(2);

  canvas2->cd(3);
  mapIJ_chip->Draw("colz");
  mapIJ->Draw("textsame");
  canvas2->cd(2);
  canvas2->Print(TString::Format("IJmapping/IJmapping_type_%s_flipx%i_flipy%i.pdf",name[iname].Data(),iflipx,iflipy));
      }
    }
  }

}

