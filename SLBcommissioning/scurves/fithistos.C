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
// #include "../conf_struct.h"

//#include "../../style/Style.C"
//#include "../../style/Labels.C"

using namespace std;


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

TF1 *FitScurve(TGraphErrors *gr, float xmin_=170, float xmax_=250, float x0=200)
{
   
  double par1=0, par2=0, par3=0;

  int imax = TMath::LocMax(gr->GetN(),gr->GetY());
  int n=gr->GetN();
  double* x = gr->GetX();
  double* y = gr->GetY();

  if(n==0 || y[imax]==0) return NULL;

  double xmin = max(x[imax],x[n/3])-3;
  double xmax = FirstZero(gr)+10;

  double ymax = y[imax];

  //  std::cout<<"  ----------- 1 "<<endl;

  TF1 *fit1 = new TF1("fit1","[0]*TMath::Erfc((x-[1])/(sqrt(2)*[2]))",xmin,xmax);

  fit1->SetParLimits(0,ymax*0.1,ymax*10);
  fit1->SetParLimits(1,xmin_,xmax_);

  fit1->SetParameter(0,ymax);
  fit1->SetParameter(1,x0);
  fit1->SetParameter(2,3); 
  gr->Fit("fit1","QR");

  par1=fit1->GetParameter(0); 
  par2=fit1->GetParameter(1); 
  par3=fit1->GetParameter(2);


  // std::cout<<par1<<" "<<par2<<" "<<par3<<endl;

  //  std::cout<<"  ----------- 2 "<<endl;

  if(par2<185) par2=190;
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


std::vector<double> MeanSigma(TGraph *scurve_threshold, TString type="SK2a",int islab=-1, int iasic=-1) {

  std::vector<double> result;

  if(islab==-1 || iasic==-1) {
    cout<<"ERROR in MeanSigma.... what slboard, iasic is this?"<<endl;
    return result;
  }
  double mean=0;
  double sigma=0;
  int max=0;
  int min=10000000;
  double nch=0;
  for(int ichn=0; ichn<64; ichn++) {
    if(detector.slab[0][islab].asu[0].skiroc[iasic].mask[ichn]==1) continue;
    Double_t x,y;
    scurve_threshold->GetPoint(ichn, x, y);
    if(y<400 && y>100) {
      mean +=  y;
      nch++;
    }
    if(y>max && y<400) max=y;
    if(y<min && y>100) min=y;
  }
  mean/=nch;
  for(int ichn=0; ichn<64; ichn++) {
    if(detector.slab[0][islab].asu[0].skiroc[iasic].mask[ichn]==1) continue;
    Double_t x,y;
    scurve_threshold->GetPoint(ichn, x, y);
    if(y<400 && y>100) sigma +=  (y-mean)*(y-mean);
  }
  sigma = sqrt( sigma/float(nch-1));


  cout<<" max:" <<max<<" optimal:"<<mean<<" std:"<<sigma<<" max-min:"<<max-min<<" type:"<<type;

  if( (max-min)>15)  max=mean+15;//3*sigma;
  if(type=="SK2") max=TMath::Max(mean,235.0);
  if(type=="SK2a") max=TMath::Max(max,225);
  result.push_back(max);

  int th[64];
  for(int ichn=0; ichn<64; ichn++) {
    th[ichn]=0;
    Double_t x,y;
    scurve_threshold->GetPoint(ichn, x, y);
    if(y<max)th[ichn]=int(max-y);
    if(th[ichn]>15) th[ichn]=15;
    if(type=="SK2") th[ichn]=0;
    if(detector.slab[0][islab].asu[0].skiroc[iasic].mask[ichn]==1) th[ichn]=0;
    result.push_back(th[ichn]);

  }



  return result;

}


void fithistos(TString filename, std::vector<int> slboards , int iteration=0, bool draw=true, bool write=true){

  if(slboards.size()==0 || slboards.size()>15) {
    cout<<"ERROR: Size of vector of ids of the slboards = "<<slboards.size()<<endl;
  }

  filename=filename;
  if(write==true) {
    read_configuration_file(TString::Format("../"+filename+"/Run_Settings_it%i.txt",iteration),false);
    cout<< " Read FILE ..... " << TString::Format("../"+filename+"/Run_Settings_it%i.txt",iteration) <<endl;
  }
  TGraphErrors *scurve[15][16][64];
  TF1 *scurvefit[15][16][64];

  TFile *file_summary = new TFile("histos/scurves_"+filename+".root" , "READ");
  file_summary->cd();

  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {
      for(int ichn=0; ichn<64; ichn++) {

       scurve[i][iasic][ichn]=(TGraphErrors*)file_summary->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
       scurvefit[i][iasic][ichn]=FitScurve(scurve[i][iasic][ichn]);
     }
   }
 }

 std::cout << scurve[1][0][5]->GetErrorX(2) << std::endl;


 TGraph *scurve_mean[15][16];
 TGraph *scurve_width[15][16];
 TGraph *scurve_threshold[15][16];
 TGraph *scurve_offset[15][16];


 for(int i=0; i<slboards.size(); i++) {
  for(int iasic=0; iasic<16; iasic++) {
    double mean[64];
    double width[64];
    double threshold[64];
    double offset[64];
    double chn[64];

    for(int ichn=0; ichn<64; ichn++) {
     mean[ichn]=-1;
     width[ichn]=-1;
     threshold[ichn]=-1;
     offset[ichn]=-1;

     chn[ichn]=ichn;

     if(scurvefit[i][iasic][ichn]!=NULL) {
       mean[ichn]=scurvefit[i][iasic][ichn]->GetParameter(1);
       width[ichn]=scurvefit[i][iasic][ichn]->GetParameter(2);
       threshold[ichn]=scurvefit[i][iasic][ichn]->GetParameter(1)+5*scurvefit[i][iasic][ichn]->GetParameter(2);
       offset[ichn]=scurvefit[i][iasic][ichn]->GetParameter(3);
     }
	  //if(ichn==37)
	  //cout<<" Result  " <<i<<  " "<<iasic <<" "<<ichn<<" "<< mean[ichn] << " " <<width[ichn] <<" "<<threshold[ichn]<<" "<<offset[ichn]<<endl;
   }
   scurve_mean[i][iasic]= new TGraph(64,chn,mean);
   scurve_width[i][iasic]= new TGraph(64,chn,width);
   scurve_threshold[i][iasic]= new TGraph(64,chn,threshold);
   scurve_offset[i][iasic]= new TGraph(64,chn,offset);

 }
}
if(write ==true) {
  for(int i=0; i<slboards.size(); i++) {
    TString type="SK2";
    // if(detector.slab[0][i].add==4 || detector.slab[0][i].add==5 || detector.slab[0][i].add==6 || detector.slab[0][i].add==7 ){type = "SK2a";}
    if(detector.slab[0][i].add==4 || detector.slab[0][i].add==5 || detector.slab[0][i].add==6 || detector.slab[0][i].add==7 || detector.slab[0][i].add==14 ){type = "SK2a";}
    //if(detector.slab[0][i].add==0 ){ type = "SK2a"; }

   for(int iasic=0; iasic<16; iasic++) {
     cout<<" Thresholds for layer:" <<i<<" skiroc:"<<iasic <<" "<<type<<" ";
     std::vector<double> mean_sigma = MeanSigma(scurve_threshold[i][iasic],type,slboards.at(i),iasic);
     cout<<endl;
     double mean=mean_sigma.at(0);
     double fine_tuning[64];
	//std::cout<<" -------------------------------------- "<< endl;
	////      std::cout<<"slab: "<<i<<" asic:"<<iasic<<"    TH:"<<mean<< endl;
     detector.slab[0][i].asu[0].skiroc[iasic].threshold_dac=mean;
     for(int ichn=0; ichn<64; ichn++) {
       fine_tuning[ichn]=mean_sigma.at(1+ichn);
	  //	if(type=="SK2a") fine_tuning[ichn]=15;//mean_sigma.at(1+ichn);
	  //else fine_tuning[ichn]=0;//mean_sigma.at(1+ichn);
	  //std::cout<<" chn:"<<ichn<<" TH : "<<fine_tuning[ichn]<<endl;
       detector.slab[0][i].asu[0].skiroc[iasic].chn_threshold_adj[ichn]=fine_tuning[ichn];
     }
   }
 }

 write_configuration_file(TString::Format("../"+filename+"/Run_Settings_comm_it%i.txt",iteration+1));

}

if(draw==true) {
  TCanvas *canvas_asic[15];
  
  // for(int i=0; i<15; i++) {
  for(int i=0; i<slboards.size(); i++) {
   canvas_asic[i]= new TCanvas(TString::Format("canvas_slbAdd%i",i),TString::Format("canvas_slbAdd%i",i),1600,1600);
   canvas_asic[i]->Divide(4,4);

   for(int iasic=0; iasic<16; iasic++) {
     canvas_asic[i]->cd(iasic+1);
	  //scurve_threshold[i][iasic]->Draw("alp");
	  // scurve_offset[i][iasic]->Draw("alp");

     for(int ichn=0; ichn<64; ichn++) {
	    //	    canvas_asic[i][iasic]->cd(1+ichn);
       scurve[i][iasic][ichn]->GetXaxis()->SetTitle("DAC");
       scurve[i][iasic][ichn]->GetYaxis()->SetTitle("Nhits");
       scurve[i][iasic][ichn]->SetTitle(TString::Format("chn %i",ichn));
       if(i==0) scurve[i][iasic][ichn]->Draw("alp");
       if(i>0) scurve[i][iasic][ichn]->Draw("lp");
     }
   }
 }

}

}



void fithistos_injection(TString filename, std::vector<int> slboards , int iteration=0, bool draw=true, bool write=true){

  if(slboards.size()==0 || slboards.size()>15) {
    cout<<"ERROR: Size of vector of ids of the slboards = "<<slboards.size()<<endl;
  }

  filename=filename;
  if(write==true) {
    read_configuration_file(TString::Format("../"+filename+"/Run_Settings_it%i.txt",iteration),false);
    cout<< " Read FILE ..... " << TString::Format("../"+filename+"/Run_Settings_it%i.txt",iteration) <<endl;
  }

  //----------------------
  TGraphErrors *scurve_1[15][16][64]={NULL};
  TF1 *scurvefit_1[15][16][64]={NULL};

  TFile *file_1_diag1 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC85_diag1.root" , "READ");
  TFile *file_1_diag2 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC85_diag2.root" , "READ");
  TFile *file_1_diag3 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC85_diag3.root" , "READ");
  TFile *file_1_diag4 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC85_diag4.root" , "READ");

  cout<<"histos/scurves_"+filename+"_PROTO15_ADC85_diag1.root"<<endl;
  cout<<"histos/scurves_"+filename+"_PROTO15_ADC85_diag2.root"<<endl;
  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {
      for(int ichn=0; ichn<64; ichn++) {
	file_1_diag1->cd();
	scurve_1[i][iasic][ichn]=(TGraphErrors*)file_1_diag1->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	
	if(scurve_1[i][iasic][ichn]->Integral()<1) {
	  file_1_diag2->cd();
	  scurve_1[i][iasic][ichn]=(TGraphErrors*)file_1_diag2->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	}
	if(scurve_1[i][iasic][ichn]->Integral()<1) {
	  file_1_diag3->cd();
	  scurve_1[i][iasic][ichn]=(TGraphErrors*)file_1_diag3->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	}
	if(scurve_1[i][iasic][ichn]->Integral()<1) {
	  file_1_diag4->cd();
	  scurve_1[i][iasic][ichn]=(TGraphErrors*)file_1_diag4->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	}
	if(scurve_1[i][iasic][ichn]->Integral()>2) scurvefit_1[i][iasic][ichn]=FitScurve(scurve_1[i][iasic][ichn],300,400,350);
      }
    }
  }
  
  //----------------------
  TGraphErrors *scurve_2[15][16][64]={NULL};
  TF1 *scurvefit_2[15][16][64]={NULL};

  TFile *file_2_diag1 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC120_diag1.root" , "READ");
  TFile *file_2_diag2 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC120_diag2.root" , "READ");
  TFile *file_2_diag3 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC120_diag3.root" , "READ");
  TFile *file_2_diag4 = new TFile("histos/scurves_"+filename+"_PROTO15_ADC120_diag4.root" , "READ");

  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {
      for(int ichn=0; ichn<64; ichn++) {
	file_2_diag1->cd();
	scurve_2[i][iasic][ichn]=(TGraphErrors*)file_2_diag1->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	
	if(scurve_2[i][iasic][ichn]->Integral()<1) {
	  file_2_diag2->cd();
	  scurve_2[i][iasic][ichn]=(TGraphErrors*)file_2_diag2->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	}
	if(scurve_2[i][iasic][ichn]->Integral()<1) {
	  file_2_diag3->cd();
	  scurve_2[i][iasic][ichn]=(TGraphErrors*)file_2_diag3->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	}
	if(scurve_2[i][iasic][ichn]->Integral()<1) {
	  file_2_diag4->cd();
	  scurve_2[i][iasic][ichn]=(TGraphErrors*)file_2_diag4->Get(TString::Format("scurve_slbAdd%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	}
	if(scurve_2[i][iasic][ichn]->Integral()>2) scurvefit_2[i][iasic][ichn]=FitScurve(scurve_2[i][iasic][ichn],420,500,440);
      }
    }
  }
  


  double mean_1[15][16]={0};
  double n_1[15][16]={0};
  double mean_2[15][16]={0};
  double n_2[15][16]={0};
  
  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {   
	for(int ichn=0; ichn<64; ichn++) {  
	  if(scurvefit_1[i][iasic][ichn]!=NULL) {
	    if(scurvefit_1[i][iasic][ichn]->GetParError(1)<20) {
	      mean_1[i][iasic]+=scurvefit_1[i][iasic][ichn]->GetParameter(1);
	      n_1[i][iasic]++;
	    }
	  }
	  if(scurvefit_2[i][iasic][ichn]!=NULL) {
	    if(scurvefit_2[i][iasic][ichn]->GetParError(1)<20) {
	      mean_2[i][iasic]+=scurvefit_2[i][iasic][ichn]->GetParameter(1);
	      n_2[i][iasic]++;
	    }
	  }
	}
    }
  }

  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {
      if(n_1[i][iasic]!=0) mean_1[i][iasic]/=n_1[i][iasic];
      else mean_1[i][iasic]=0;
      if(n_2[i][iasic]!=0) mean_2[i][iasic]/=n_2[i][iasic];
      else mean_2[i][iasic]=0;
    }
  }

  // HARD CODED ADC-to-MIP dependence.. based on the cosmic run_050010
  // 320 um ->  1MIP = 58ADC
  // 500 um ->  1MIP = 70ADC
  cout<<" FIRST CASE"<<endl;
  double MIP[15][16]={0};
  double injected1=85;
  double injected2=120;

  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {
      MIP[i][iasic]=58;
      if(i>3 && i<8) MIP[i][iasic]=70;
      if(i==2 && iasic<4) MIP[i][iasic]=70;
    }
  }

  float dac[15]={0};
  float nfits[15]={0};
  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {

      double x1=injected1/MIP[i][iasic];
      double x2=injected2/MIP[i][iasic];
      float optimal = 0; 
      if(n_1[i][iasic]>1 && n_2[i][iasic]>1) {
	dac[i] += mean_1[i][iasic] - ( (x1-0.5)/ (x2-x1) ) * ( mean_2[i][iasic] - mean_1[i][iasic]);
	nfits[i]++;
      }
    }
    dac[i]/=nfits[i];
  }

  float dac_chip[15][16]={0};
  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {
      double x1=injected1/MIP[i][iasic];
      double x2=injected2/MIP[i][iasic];
      
      if(n_1[i][iasic]>1 && n_2[i][iasic]>1) {
	dac_chip[i][iasic] = mean_1[i][iasic] - ( (x1-0.5)/ (x2-x1) ) * ( mean_2[i][iasic] - mean_1[i][iasic]);
      } else dac_chip[i][iasic]=dac[i];
    }
  }

  if(write ==true) {
    for(int i=0; i<slboards.size(); i++) {	       
      for(int iasic=0; iasic<16; iasic++) {

	cout<<i<<" "<<iasic<<" "<<detector.slab[0][i].asu[0].skiroc[iasic].threshold_dac<<endl;
	if(detector.slab[0][i].add==4 || detector.slab[0][i].add==5 || detector.slab[0][i].add==6 || detector.slab[0][i].add==7 ) {
	  detector.slab[0][i].asu[0].skiroc[iasic].threshold_dac=dac_chip[i][iasic]+15;
	  for(int ichn=0; ichn<64; ichn++) detector.slab[0][i].asu[0].skiroc[iasic].chn_threshold_adj[ichn]=15;
	} else {
	  detector.slab[0][i].asu[0].skiroc[iasic].threshold_dac=dac_chip[i][iasic];
	  for(int ichn=0; ichn<64; ichn++) detector.slab[0][i].asu[0].skiroc[iasic].chn_threshold_adj[ichn]=0;
	}
	cout<<i<<" "<<iasic<<" "<<detector.slab[0][i].asu[0].skiroc[iasic].threshold_dac<<endl;
      }
    }
      
    write_configuration_file(TString::Format("../"+filename+"/Run_Settings_comm_it%i.txt",iteration+1));
    cout<<TString::Format("../"+filename+"/Run_Settings_comm_it%i.txt",iteration+1)<<endl;
  }



}



