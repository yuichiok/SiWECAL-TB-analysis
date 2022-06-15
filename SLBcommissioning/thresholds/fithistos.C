#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "TF1.h"

#include "../conf_struct.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;

int nslabs=15;
int nchips=16;
int nchannels=64;
int nchannels_fit=6;


Float_t charge_low[16][64]; //daughter, slab, asic, chn, DAC
Float_t charge_high[16][64]; //daughter, slab, asic, chn, DAC
void ReadCharge(TString filename)
  
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
  for(int i=0; i<nchips; i++) {
    for(int j=0; j<nchannels; j++) {
      charge_low[i][j] = 0.;
      charge_high[i][j] = 0.;
    }
  }
  
  std::string line;

  // check Firstline
  getline(reading_file, line);

  int ichn=-1, ski=-1;
  int gain=-1;
  while( getline(reading_file,line) ) {
    search_string_nocolon(line,"Skiroc",ski);
    search_string_nocolon(line,"Ch",ichn);

    if(ski>-1 && ichn>-1) {
      //Low gain
      getline(reading_file,line);
      search_string_nocolon(line,"Gain",gain);
      getline(reading_file,line);
      std::string space_delimiter = " ";
      size_t pos = 0;
      int i = 0;
      int scas=0;
      while ((pos = line.find(space_delimiter)) != string::npos) {
	if(i>14) break;
	float tmp=0;
	try {
	  tmp  = std::stof(line.substr(0, pos));
	} catch(std::invalid_argument e){ tmp=0;};
	  
	if(tmp>0) {
	  charge_low[ski][ichn]+=tmp;
	  scas++;
	}
	line.erase(0, pos + space_delimiter.length());
	i++;
      }
      //High gain
      //cout<<"HG "<<line<<endl;

      getline(reading_file,line);
      search_string_nocolon(line,"Gain",gain);
      getline(reading_file,line);
      space_delimiter = " ";
      pos = 0;
      i = 0;
      scas=0;
      //	cout<<line<<endl;
      while ((pos = line.find(space_delimiter)) != string::npos) {
	if(i>15) break;
	float tmp=0;
	try {
	  tmp  = std::stof(line.substr(0, pos));
	} catch(std::invalid_argument e){ tmp=0;};
	if(tmp>0) {
	  charge_high[ski][ichn]+=tmp;
	  scas++;
	}
	line.erase(0, pos + space_delimiter.length());
	i++;
      }
      if(scas>0) charge_high[ski][ichn]/=scas;
      //      cout<<scas<<" "<<charge_high[ski][ichn]<<endl;
    }
  }


}

float FirstZero(TGraphErrors *gr)
{


  int imax = TMath::LocMax(gr->GetN(),gr->GetY()); 
  Double_t xmax,ymax;
  gr->GetPoint(imax, xmax, ymax);

  //if(imax==0) return 0;

  for(int i=imax; i<gr->GetN(); i++) {
    gr->GetPoint(i, xmax, ymax);
    if(ymax<0.05) return xmax;
  }

  return 0;
  
}

TF1 *FitScurve(TGraphErrors *gr, float xmin_=200, float xmax_=350, float x0=200)
{
   
  double par1=0, par2=0, par3=0;

  int imax = TMath::LocMax(gr->GetN(),gr->GetY());
  int n=gr->GetN();
  double* x = gr->GetX();
  double* y = gr->GetY();

  if(n==0 || y[imax]==0) return NULL;

  double xmin = FirstZero(gr)-15;//max(x[imax],x[n/3])-5;
  double xmax = FirstZero(gr)+15;

  double ymax = y[imax];

  //  std::cout<<"  ----------- 1 "<<endl;

  TF1 *fit1 = new TF1("fit1","[0]*TMath::Erfc((x-[1])/(sqrt(2)*[2]))",xmin,xmax);

  fit1->SetParLimits(0,ymax*0.1,ymax*10);
  fit1->SetParLimits(1,xmin_,xmax_);

  fit1->SetParameter(0,ymax);
  // fit1->SetParameter(1,250);
  fit1->SetParameter(2,3); 
  gr->Fit("fit1","QR");

  par1=fit1->GetParameter(0); 
  par2=fit1->GetParameter(1); 
  par3=fit1->GetParameter(2);


  // std::cout<<par1<<" "<<par2<<" "<<par3<<endl;

  //  std::cout<<"  ----------- 2 "<<endl;

  //if(par2<185) par2=190;
  //  if(par3<1) par3=1.5;
  xmin=max(par2-8*par3,x[n/3]);
  xmax=par2+20*par3;
  TF1 *fit2 = new TF1("fit2","[0]*TMath::Erfc((x-[1])/[2]+[3])",xmin,xmax);

  fit2->SetParLimits(0,par1*0.3,par1+par1*0.5);
  fit2->SetParLimits(1,xmin_,xmax_);
  fit2->SetParLimits(2,1.5,10);

  /* fit2->SetParameter(0,0.5*ymax);
     fit2->SetParameter(1,200);
     fit2->SetParameter(2,10); */
  fit2->SetParameter(0,par1);
  fit2->SetParameter(1,par2);
  fit2->SetParameter(2,par3);
  gr->Fit("fit2","QR");
    
  return fit2;
}

void fithistos(){


  read_configuration_file("Run_Settings_090185.txt",false);
  double dac_value[15][16]={0};
  double dac_value_p10[15][16]={0};
  double dac_value_p20[15][16]={0};
  for(int islab=0; islab<nslabs; islab++) {
    for(int ichip=0; ichip<nchips; ichip++) {
      double nchn=0;
      double th=0;
      for(int ichn=0; ichn<nchannels; ichn++) {
	if(detector.slab[0][islab].asu[0].skiroc[ichip].mask[ichn]==0 ){
	  th+=(detector.slab[0][islab].asu[0].skiroc[ichip].threshold_dac - detector.slab[0][islab].asu[0].skiroc[ichip].chn_threshold_adj[ichn]);
	  nchn++;
	}
      }
      if(nchn!=0) th/=nchn;
      dac_value[islab][ichip]=th;
      dac_value_p10[islab][ichip]=th+10.;
      dac_value_p20[islab][ichip]=th+20.;
    }
  }

  int slabsorting[16]={4,8,3,7,2,6,1,5,12,16,11,15,10,14,9,13};

  double xtest2[2]={0.5,0.5};
  double ytest2[2]={0,400};
  TGraph* g_test2= new TGraph(2,xtest2,ytest2);

  double xtest3[2]={0.75,0.75};
  double ytest3[2]={0,400};
  TGraph* g_test3= new TGraph(2,xtest3,ytest3);

  double xtest4[2]={1,1};
  double ytest4[2]={0,400};
  TGraph* g_test4= new TGraph(2,xtest4,ytest4);
  
  TCanvas * c_dac_vs_mip[15];
  for(int i=0; i<nslabs; i++) {
    c_dac_vs_mip[i] = new TCanvas(TString::Format("slab_%i",i),TString::Format("slab_%i",i),1200,900);
    c_dac_vs_mip[i]->Divide(4,4);
  }

  TF1 *scurvefit[10][15][16][4];

    
  for(int ifile=0; ifile<4; ifile++) {
    TString filename="injection_3p0high_chn0to8";
    if(ifile==1) filename="injection_3p05high_chn0to8";
    if(ifile==2) filename="injection_3p1high_chn0to8";
    if(ifile==3) filename="injection_3p2high_chn0to8";
    
    
    TGraphErrors *scurve[15][16][4];
    
    TFile *file_summary = new TFile("histos/scurves_"+filename+".root" , "READ");
    
    file_summary->cd();
    for(int i=0; i<nslabs; i++) {
      for(int iasic=0; iasic<nchips; iasic++) {
	for(int ichn=0; ichn<nchannels_fit; ichn++) {
	  
	  scurve[i][iasic][ichn]=(TGraphErrors*)file_summary->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",i,iasic,ichn));
	  scurvefit[ifile][i][iasic][ichn]=FitScurve(scurve[i][iasic][ichn]);
	}
      }
    }
  }
  cout<<" AAA "<<endl;  
  TFile *file_summary2 = new TFile("../../mip_calib/MIPSummary_pedestalsubmode1_raw_siwecal_90021to90070_highgain.root" , "READ");
  for(int i=0; i<nslabs; i++) {

    for(int iasic=0; iasic<nchips; iasic++) {

      TH1F *h_temp=(TH1F*)file_summary2->Get(TString::Format("layer_%i/mip_high_layer%i_chip%i",i,i,iasic));
      TF1 * f_temp=h_temp->GetFunction(TString::Format("Fitfcn_mip_high_layer%i_chip%i",i,iasic));
      float mpv=9999999.;
      float empv=0.;

      if(f_temp!=NULL) {
	mpv=f_temp->GetParameter(1);
	if(mpv<2) mpv=999999.;
	empv=f_temp->GetParameter(0);
      }

      c_dac_vs_mip[i]->cd(slabsorting[iasic]);
      TGraphErrors * dac_vs_mip[8];
      for(int ichn=0; ichn<nchannels_fit; ichn++) {
	double x[10]={0}, y[10]={0}, ex[10]={0}, ey[10]={0};

	int nfiles=0;
	for(int ifile=0; ifile<4; ifile++) {
	  TString filename="injection_3p0high_chn0to8";
	  if(ifile==1) filename="injection_3p05high_chn0to8";
	  if(ifile==2) filename="injection_3p1high_chn0to8";
	  if(ifile==3) filename="injection_3p2high_chn0to8";
	  TString filename_input =filename+TString::Format("/ChargeMeanValues_Core0_SlabAdd%i_Asu0.txt",i);
	  ReadCharge(filename_input);
	  if(scurvefit[ifile][i][iasic][ichn]!=NULL) {
	    y[ifile]=scurvefit[ifile][i][iasic][ichn]->GetParameter(1);
	    ey[ifile]=scurvefit[ifile][i][iasic][ichn]->GetParError(1);
	    
	    x[ifile]=charge_high[iasic][ichn]/mpv;
	    ex[ifile]=empv/mpv;//chargefit[ifile][i][iasic][ichn]->GetParameter(2);
	    nfiles++;
	  } 
	}
                
	dac_vs_mip[ichn]= new TGraphErrors(nfiles,x,y,ex,ey);
	dac_vs_mip[ichn]->SetLineColor(ichn+1);
	dac_vs_mip[ichn]->SetMarkerColor(ichn+1);
	dac_vs_mip[ichn]->SetMarkerStyle(ichn+20);
	if(ichn==0) {
	  double xtest[2]={0,2.5};
	  double ytest[2]={180,300};
	  TGraph* g_test= new TGraph(2,xtest,ytest);
	  g_test->SetLineColor(0);
	  g_test->SetTitle(TString::Format("Layer:%i, ASIC:%i",i,iasic));
	  g_test->GetYaxis()->SetTitle("DAC");
	  g_test->GetYaxis()->SetTitleOffset(1);
	  g_test->GetYaxis()->SetTitleSize(0.05);
	  g_test->GetYaxis()->SetLabelSize(0.05);
	  g_test->GetXaxis()->SetTitle("MIP");
	  g_test->GetXaxis()->SetTitleSize(0.05);
	  g_test->GetXaxis()->SetLabelSize(0.05);

	  double xdac[2]={0,2.5};
          double ydac[2]={dac_value[i][iasic],dac_value[i][iasic]};
          TGraph* g_dac= new TGraph(2,xdac,ydac);

	  double xdac2[2]={0,2.5};
          double ydac2[2]={dac_value_p10[i][iasic],dac_value_p10[i][iasic]};
          TGraph* g_dac2= new TGraph(2,xdac2,ydac2);

	  double xdac3[2]={0,2.5};
          double ydac3[2]={dac_value_p20[i][iasic],dac_value_p20[i][iasic]};
          TGraph* g_dac3= new TGraph(2,xdac3,ydac3);

	  g_test->Draw("al");
	  g_dac->SetLineColor(2);
          g_dac->SetLineStyle(1);
          g_dac->SetLineWidth(3);
          g_dac->Draw("l");

	  g_dac2->SetLineColor(2);
          g_dac2->SetLineStyle(1);
          g_dac2->SetLineWidth(2);
          g_dac2->Draw("l");

	  g_dac3->SetLineColor(2);
          g_dac3->SetLineStyle(1);
          g_dac3->SetLineWidth(1);
          g_dac3->Draw("l");

	  g_test2->SetLineColor(1);
	  g_test2->SetLineStyle(2);
	  g_test2->SetLineWidth(3);
	  g_test2->Draw("l");
	  g_test3->SetLineColor(1);
          g_test3->SetLineStyle(2);
          g_test3->SetLineWidth(2);
          g_test3->Draw("l");
	  g_test4->SetLineColor(1);
	  g_test4->SetLineStyle(2);
	  g_test4->SetLineWidth(1);
	  g_test4->Draw("l");
	}
	dac_vs_mip[ichn]->Draw("lp");
      }
    }
  }

  
  for(int i=0; i<nslabs; i++) {
    c_dac_vs_mip[i]->Print(TString::Format("slab_%i.pdf",i));
  }
}


