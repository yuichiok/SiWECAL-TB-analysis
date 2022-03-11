#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2F.h"
#include "TString.h"
#include "TImage.h"
#include "TMinuit.h"
#include "../../../include/utils.h"

double cov_data[64][64]={0};
double n_data[64][64]={0};
double sigma_i[64]={0};
double sigma_c1[64]={0};
double sigma_c2[64]={0};
double sigma_c3[64]={0};
double sigma_c4[64]={0};
double esigma_i[64]={0};
double esigma_c1[64]={0};
double esigma_c2[64]={0};
double esigma_c3[64]={0};
double esigma_c4[64]={0};

double pedestal[64]={0};
double epedestal[64]={0};


//______________________________________________________________________________


Double_t func(int i,int j, Double_t *par)
{
  Double_t value;
  //if(i==j) value=(par[i]*par[j]+par[i+64]*par[j+64]);
  //  else value=(par[i+64]*par[j+64]);

  //if(i==j) value=(par[i]*par[j]+par[i+64]*par[j+64]+par[i+128]*par[j+128]);
  //else value=(par[i+64]*par[j+64]+par[i+128]*par[j+128]);
  // if(i==j) value=par[i]*par[j]+par[i+64]*par[j+64]+par[i+128]*par[j+128];
  // else value=par[i+64]*par[j+64]+par[i+128]*par[j+128];
  if(i==j) value=(par[i]*par[j]+par[i+64]*par[j+64]+par[i+128]*par[j+128]);//+par[i+192]*par[j+192]+par[i+256]*par[j+256]);
  else value=(par[i+64]*par[j+64]+par[i+128]*par[j+128]);//+par[i+192]*par[j+192]+par[i+256]*par[j+256]);
  return value;
 
}


//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;

  for (int i=0;i<64; i++) {
    for (int j=0;j<64; j++) {
      delta = pow(func(i,j,par)-cov_data[i][j],2);
      chisq += delta;
    }
  }
  f = chisq;

}


TH2F* OpenCovMatrix(TString filename="3GeVMIPscan_scagt0",TString gain="highgain", int layer=5, int chip=0)
{

  TFile * file = TFile::Open("../../results_noise/NoiseCovariance_"+filename+".root");
  cout<<"../../results_noise/NoiseCovariance_"+filename+".root"<<endl;
  TH2F* cov_temp=(TH2F*)file->Get(TString::Format("layer_%i/cov_unnorm_layer%i_chip%i",layer,layer,chip));
  cout<<TString::Format("layer_%i/cov_unnorm_layer%i_chip%i",layer,layer,chip)<<endl;
  TH2F* norm_temp=(TH2F*)file->Get(TString::Format("layer_%i/nevents_layer%i_chip%i",layer,layer,chip));
  cout<<TString::Format("layer_%i/nevents_layer%i_chip%i",layer,layer,chip)<<endl;
  TH2F* cov= new TH2F(TString::Format("%s_%s_layer%i_chip%i",filename.Data(),gain.Data(),layer,chip), TString::Format("%s_%s_layer%i_chip%i",filename.Data(),gain.Data(),layer,chip),64,-0.5,63.5,64,-0.5,63.5);
  
  
  if(cov_temp==NULL || norm_temp==NULL){//|| norm_temp==NULL) {
    cout<< "cov_temp==NULL || norm_temp==NULL"<<endl;
    file->Close();
    return NULL;
  }

  double integral=0;
  for (int i=0; i<64; i++) {
    integral+=norm_temp->GetBinContent(i+1,i+1);
  }
  if ( integral/64. <100) {
    delete file;
    cout<< "integral/64. <500, integral="<< integral/64.<<endl;
    return NULL;
  }
  
  for (int i=0; i<64; i++) {
    TH1F* htemp=(TH1F*)file->Get(TString::Format("layer_%i/h_ped_layer%i_chip%i_chn%i",layer,layer,chip,i));
    if(htemp==NULL) {
      cout<<"htemp == NUL"<<endl;
      continue;
    }
    if(htemp->GetEntries()>100) {
      htemp->GetXaxis()->SetRangeUser(htemp->GetMean()-20,htemp->GetMean()+20);
      pedestal[i]=htemp->GetMean();
      epedestal[i]=htemp->GetRMS();
    
      for (int j=0; j<64; j++) {
	if(norm_temp->GetBinContent(i+1,j+1)>0  ) {
	  cov->SetBinContent(i+1,j+1,cov_temp->GetBinContent(i+1,j+1)/norm_temp->GetBinContent(i+1,j+1));
	  n_data[i][j]=norm_temp->GetBinContent(i+1,j+1);
	}
	else {
	  cov->SetBinContent(i+1,j+1,0);
	  n_data[i][j]=0;
	}
	cov_data[i][j]=cov->GetBinContent(i+1,j+1);
      }
    }
  }
  delete file;
  return cov;
}



int Fit(TString name="3GeVMIPscan", TString gain="highgain", int layer=0, int chip=0){
  TH2F* cov_matrix = OpenCovMatrix(name,gain,layer,chip);
  if(cov_matrix==NULL) return 0;
  // RooFit::Offset(true);
  //for(int i=0; i<4; i++) std::cout<<z[i]<<" "<<rho[i]<<"    ....     "<<std::endl;
  Int_t nparameters=192;
  TMinuit *gMinuit = new TMinuit(nparameters);  //initialize TMinuit with a maximm of 5 params
  gMinuit->SetFCN(fcn);

  //gMinuit->SetFCN(fcn_cov);

  gMinuit->mnhelp("*");
  gMinuit->mnhelp("MIGRAD");
  //m.mnhelp("MIGrad")
  Double_t arglist[nparameters];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,0,ierflg);

  // Set starting values and step sizes for parameters
  Double_t vstart[nparameters];
  Double_t step[nparameters];
  for(int i=0; i<nparameters; i++) {
    if(i<64) {
      vstart[i]=3.0;
      step[i]=0.001;
      gMinuit->mnparm(i, TString::Format("I_%i",i), vstart[i], step[i], 1.0,100,ierflg);
    } else if(i<128) {
      vstart[i]=0.5;
      step[i]=0.001;
      gMinuit->mnparm(i, TString::Format("C1_%i",i-64), vstart[i], step[i], 0.0,5,ierflg);
    } else if(i<192) {
      vstart[i]=0.1;
      step[i]=0.001;
      gMinuit->mnparm(i, TString::Format("C2_%i",i-128), vstart[i], step[i], 0.0,05,ierflg);
    } else if(i<256)  {
      vstart[i]=0.05;
      step[i]=0.0001;
      gMinuit->mnparm(i, TString::Format("C3_%i",i-192), vstart[i], step[i], 0.,5,ierflg);
    } else {
      vstart[i]=0.005;
      step[i]=0.001;
      gMinuit->mnparm(i, TString::Format("C4_%i",i-256), vstart[i], step[i], 0,1,ierflg);
    }
  }

  // Now data for minimization step
  arglist[0] = 20000;//1000;//10000; //2000 for the HG
  arglist[1] = 1;
  gMinuit->SetErrorDef(1);
  gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);

  //    Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat); 
  gMinuit->mnprin(3,amin); 

  Int_t npars = gMinuit->GetNumPars();
  Double_t *covar = new Double_t[npars*npars];
  gMinuit->mnemat(covar,npars); 

  //  cout<<covar[0][0]<<endl;

  for(int i=0; i<64; i++) {
    sigma_i[i]=0;
    sigma_c1[i]=0;
    sigma_c2[i]=0;
    gMinuit->GetParameter(i,sigma_i[i],esigma_i[i]);
    //    double temp_c[2]={0};
    gMinuit->GetParameter(i+64,sigma_c1[i],esigma_c1[i]);
    gMinuit->GetParameter(i+128,sigma_c2[i],esigma_c2[i]);
    // gMinuit->GetParameter(i+192,sigma_c3[i],esigma_c3[i]);
    //gMinuit->GetParameter(i+256,temp_c[3],esigma_c4[i]);
    //int n = sizeof(temp_c)/sizeof(temp_c[0]); 
    //sort(temp_c, temp_c+n);
    // cout<<sigma_c1[i]<<endl;
    //sigma_c1[i]=temp_c[1];
    //sigma_c2[i]=temp_c[0];
    // sigma_c3[i]=temp_c[1];
    // sigma_c4[i]=temp_c[0];
    
    
  }

  return 1;
}





