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

TF1 *FitScurve(TGraphErrors *gr)
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
  fit1->SetParLimits(1,170,250);

  fit1->SetParameter(0,ymax);
  fit1->SetParameter(1,200);
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
  fit2->SetParLimits(1,170,240);
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

void fithistos(TString filename, std::vector<int> slboards , int iteration=0){

  if(slboards.size()==0 || slboards.size()>15) {
    cout<<"ERROR: Size of vector of ids of the slboards = "<<slboards.size()<<endl;
  }

  filename=filename;
  read_configuration_file(TString::Format("../"+filename+"/Run_Settings_it%i.txt",iteration),false);
  cout<< " Read FILE ..... " << TString::Format("../"+filename+"/Run_Settings_it%i.txt",iteration) <<endl;
  TGraphErrors *scurve[15][16][64];
  TF1 *scurvefit[15][16][64];

  TFile *file_summary = new TFile("histos/scurves_"+filename+".root" , "READ");
  file_summary->cd();

  for(int i=0; i<slboards.size(); i++) {
    for(int iasic=0; iasic<16; iasic++) {
      for(int ichn=0; ichn<64; ichn++) {
	scurve[i][iasic][ichn]=(TGraphErrors*)file_summary->Get(TString::Format("scurve_slboard%i_chip%i_chn%i",slboards.at(i),iasic,ichn));
	scurvefit[i][iasic][ichn]=FitScurve(scurve[i][iasic][ichn]);
      }
    }
  }

  
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

  for(int i=0; i<slboards.size(); i++) {
    TString type="SK2";
    if(slboards.at(i)==3 || slboards.at(i)==7 || slboards.at(i)==8|| slboards.at(i)==12) type="SK2a";			       
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
	
  /*  TCanvas *canvas_asic[15];
  
  for(int i=0; i<nslabs; i++) {
    canvas_asic[i]= new TCanvas(TString::Format("canvas_slboard%i",i),TString::Format("canvas_slboard%i_chip%i",i),1600,1600);
    canvas_asic[i]->Divide(4,4);
    
    for(int iasic=0; iasic<16; iasic++) {
      canvas_asic[i]->cd(iasic+1);
      scurve_threshold[i][iasic]->Draw("alp");
      //scurve_offset[i][iasic]->Draw("alp");

   //      for(int ichn=0; ichn<64; ichn++) {
//	canvas_asic[i][iasic]->cd(1+ichn);
//	scurve[i][iasic][ichn]->GetXaxis()->SetTitle("DAC");
//	scurve[i][iasic][ichn]->GetYaxis()->SetTitle("Nhits");
//	scurve[i][iasic][ichn]->SetTitle(TString::Format("chn %i",ichn));
//	scurve[i][iasic][ichn]->Draw("al");
//	}
    }
  }
*/
    
  
}





