#define savePedestal_cxx
#include "savePedestal.h"
#include <TPaveStats.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <TFitResult.h>
#include <TSpectrum.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;



void savePedestal::FindMasked(TString dif)
{

  //function that reads a root file and check which channels have or not signal (hit bit ==1).
  //should be used for root files that contain a full position scan, in order to provide meaninful list of masked channels.

  ofstream fout_masked("masked_"+dif+".log",ios::out);

  fout_masked<<"#list of masked channels, per dif: "<<dif<<endl;
  fout_masked<<"#chip channel mask (0=not masked, 1=masked)"<<endl;

  
  bool global = true;
  int maxnhit=32;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  std::vector<std::vector<TH1F*> > h_charge_channel;

  for(int ichip=0; ichip<16; ichip++) {
    
    std::vector<TH1F*> chargetemp_sca;
    for(int ichn=0; ichn<64; ichn++) {
      TH1F *tmp_charge_chn = new TH1F(TString::Format("tmp_charge_chip%i_chn%i",ichip,ichn),TString::Format("tmp_charge_chip%i_chn%i",ichip,ichn),4096,0.5,4096.5);
      chargetemp_sca.push_back(tmp_charge_chn);
    }
    h_charge_channel.push_back(chargetemp_sca);
  }


  
  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    for(int ichip=0; ichip<16; ichip++) {
      for(int isca=0; isca<15; isca++) {
	bool gooddata=true;
	if(global == true) {
	  int ntaggedasbad = 0;
	  for(int ichn=0; ichn<64; ichn++) {
	    if(charge_hiGain[ichip][isca][ichn]<10 && charge_hiGain[ichip][isca][ichn]>-1 ) 
	      ntaggedasbad++;
	  }//ichn 
	  if ( ntaggedasbad > 0) gooddata=false;
	}
	
	for(int ichn=0; ichn<64; ichn++) {

	  bool selection=false;
	  if(charge_hiGain[ichip][isca][ichn]>10 && badbcid[ichip][isca]==0 && nhits[ichip][isca]<maxnhit+1) selection=true;    
	  if(gain_hit_high[ichip][isca][ichn]==1 && selection==true && gooddata==true)
	    h_charge_channel.at(ichip).at(ichn)->Fill(charge_hiGain[ichip][isca][ichn]);
	}//ichn
           
      }//isca
    }//ichip 
  }  // end first loop analysis to fill pedestal historgrams


  //------------------------------------------------------------------
  // do signal analysis (chip/channel/sca based) analysis
  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      
      fout_masked << ichip <<" " <<ichn<< " "; 
      //minimum of 100 entries per SCA
      if(h_charge_channel.at(ichip).at(ichn)->GetEntries()>100) fout_masked<<"0"<<endl;
      else fout_masked<<"1"<<endl;
      
    }
  }


}

void savePedestal::ReadMap(TString filename) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      map_pointX[i][j] = -1000.;
      map_pointY[i][j] = -1000.;
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Float_t tmp_x0 = 0 ,tmp_y0 = 0 , tmp_x = 0 , tmp_y = 0 ;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y ;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x ;
    map_pointY[tmp_chip][tmp_channel] = -tmp_y ;
   }

}

void savePedestal::ReadMasked(TString filename) 
{

  std::ifstream reading_file(filename);
  if(!reading_file){
    cout<<" dameyo - damedame"<<endl;
  }
  for(int i=0; i<16; i++) {
    for(int j=0; j<64; j++) {
      masked[i][j] = 0;
    }
  }

  Int_t tmp_chip = 0,tmp_channel = 0;
  Int_t tmp_masked = 0;
  TString tmpst;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst ;
  reading_file >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst >> tmpst  ;

  while(reading_file){
    reading_file >> tmp_chip >> tmp_channel >> tmp_masked ;
    masked[tmp_chip][tmp_channel] = tmp_masked ;
   }

}

void savePedestal::PedestalAnalysis(TString dif,TString map_filename="../fev10_chip_channel_x_y_mapping.txt")
{

  //Read the channel/chip -- x/y mapping
  ReadMap(map_filename);
  //Read the list of masked channels
  ReadMasked("masked_"+dif+".log");

  ofstream fout_ped("Pedestal_"+dif+".log",ios::out);

  fout_ped<<"#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis : "<<dif<<endl;
  fout_ped<<"#chip channel ped0 eped1 ped1 eped1 ... ped14 ped14 (all SCA)"<<endl;

  
  bool global = true;
  int maxnhit=32;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  std::vector<std::vector<std::vector<float> > > ped_mean;
  std::vector<std::vector<std::vector<float> > > ped_rms;

  for(int ichip=0; ichip<16; ichip++) {
    ped_mean.push_back(std::vector<std::vector<float> >() );
    ped_rms.push_back(std::vector<std::vector<float> >() );

    for(int ichn=0; ichn<64; ichn++) {
      ped_mean.at(ichip).push_back(std::vector<float> () );
      ped_rms.at(ichip).push_back(std::vector<float> () );
    }
  }
  //



  TH2F* pedestal_map[15];
  TH2F* pedestal_width_map[15];
  TH2F* pedestal_error_map[15];
  TH2F* pedestal_chi2ndf_map[15];
  TH2F* pedestal_npeaks_map[15];
  TH2F* pedestal_entries_map[15];

  for(int isca=0; isca<15; isca++) {
    pedestal_map[isca]= new TH2F(TString::Format("pedestal_map_sca%i",isca),TString::Format("pedestal_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_width_map[isca]= new TH2F(TString::Format("pedestal_width_map_sca%i",isca),TString::Format("pedestal_width_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_error_map[isca]= new TH2F(TString::Format("pedestal_error_map_sca%i",isca),TString::Format("pedestal_error_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_chi2ndf_map[isca]= new TH2F(TString::Format("pedestal_chi2ndf_map_sca%i",isca),TString::Format("pedestal_chi2ndf_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_npeaks_map[isca]= new TH2F(TString::Format("pedestal_npeaks_map_sca%i",isca),TString::Format("pedestal_npeaks_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_entries_map[isca]= new TH2F(TString::Format("pedestal_entries_map_sca%i",isca),TString::Format("pedestal_entries_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);

  }

  // --------------------
  // all sca
  //  std::vector<TCanvas*> chip;
  std::vector<std::vector<std::vector<TH1F*> > > ped_sca ;
  std::vector<TH1F*> pedestal_chip ;

  for(int ichip=0; ichip<16; ichip++) {
    TH1F *ped_chip = new TH1F(TString::Format("ped_chip%i",ichip),TString::Format("ped_chip%i",ichip),500,0.5,500.5);
    pedestal_chip.push_back(ped_chip);
    
    std::vector<std::vector<TH1F*> >pedtemp_sca;
    for(int ichn=0; ichn<64; ichn++) {
      std::vector<TH1F*> pedtemp_sca2;
      for(int isca=0; isca<15; isca++) {
	TH1F *ped_sca2 = new TH1F(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca),TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca),500,0.5,500.5);
	pedtemp_sca2.push_back(ped_sca2);
      }
      pedtemp_sca.push_back(pedtemp_sca2);
    }
    ped_sca.push_back(pedtemp_sca);
  }


  
  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    for(int ichip=0; ichip<16; ichip++) {

      for(int isca=0; isca<15; isca++) {

	bool gooddata=true;
	if(global == true) {
	  int ntaggedasbad = 0;
	  for(int ichn=0; ichn<64; ichn++) {
	    if(charge_hiGain[ichip][isca][ichn]<10 && charge_hiGain[ichip][isca][ichn]>-1 ) 
	      ntaggedasbad++;
	  }//ichn 
	  if ( ntaggedasbad > 0) gooddata=false;
	}
	
	for(int ichn=0; ichn<64; ichn++) {

	  bool selection=false;
	  if(charge_hiGain[ichip][isca][ichn]>10 && (badbcid[ichip][isca]==0 || badbcid[ichip][isca]==0 ) && nhits[ichip][isca]<maxnhit+1) selection=true;

	  if(masked[ichip][ichn]==1) selection=false;
 
	  if(gain_hit_high[ichip][isca][ichn]==0 && selection==true && gooddata==true)
	    ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_hiGain[ichip][isca][ichn]);

	  /* if(charge_hiGain[ichip][isca][ichn]>10 && gain_hit_high[ichip][isca][ichn]==1 && badbcid[ichip][isca]==1 && nhits[ichip][isca]<maxnhit+1 && gooddata==true) {
	    std::cout<<ichip<< " " << ichn<< " "<< isca << " " <<charge_hiGain[ichip][isca][ichn]<< " " << gain_hit_high[ichip][isca][ichn] << " " ;
	    if(isca>0) std::cout<< charge_hiGain[ichip][isca-1][ichn]<< " " << gain_hit_high[ichip][isca-1][ichn] << " " << badbcid[ichip][isca-1]<< " "<< bcid[ichip][isca-1]<< " ";
	    if(isca<15) std::cout<< charge_hiGain[ichip][isca+1][ichn]<< " "<< gain_hit_high[ichip][isca+1][ichn] << " " << badbcid[ichip][isca+1]<< " "<< bcid[ichip][isca+1]<< " ";
	    if(isca<14) std::cout<< charge_hiGain[ichip][isca+2][ichn]<< " "<< gain_hit_high[ichip][isca+2][ichn] << " " << badbcid[ichip][isca+2]<< " "<< bcid[ichip][isca+2]<< " ";
	    if(isca<13) std::cout<< charge_hiGain[ichip][isca+3][ichn]<< " "<< gain_hit_high[ichip][isca+3][ichn] << " " << badbcid[ichip][isca+3]<< " "<< bcid[ichip][isca+3]<< " ";
	    if(isca<12) std::cout<< charge_hiGain[ichip][isca+4][ichn]<< " "<< gain_hit_high[ichip][isca+4][ichn] << " " << badbcid[ichip][isca+4]<< bcid[ichip][isca+4]<< " ";

	    std::cout<<std::endl;
	    }*/
	}
           
      }//isca

    }//ichip 
   
  }  // end first loop analysis to fill pedestal historgrams


  TFile *pedfile = new TFile("Pedestal_"+dif+".root" , "RECREATE");
  pedfile->cd();

  
  // do pedestal (chip/channel/sca based) analysis
  for(int ichip=0; ichip<16; ichip++) {
    for(int ichn=0; ichn<64; ichn++) {
      
      fout_ped << ichip <<" " <<ichn<< " "; 
      for(int isca=0; isca<15; isca++) {

	ped_sca.at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca.at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca.at(ichip).at(ichn).at(isca)->Write();

	pedestal_entries_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , ped_sca.at(ichip).at(ichn).at(isca)->GetEntries());

	
	if(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries()> 500 ){ //max_entries/2 ) {
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(ichip).at(ichn).at(isca),2,"",0.05); 
	  pedestal_npeaks_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , npeaks);

	  if(npeaks == 1) {
	    Double_t *mean_peak=s->GetPositionX();
	    
	    TF1 *f0 = new TF1("f0","gaus",mean_peak[0]-2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS(),mean_peak[0]+2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS());
	    ped_sca.at(ichip).at(ichn).at(isca)->Fit("f0","RQNOC");

	    TF1 *f1 = new TF1("f1","gaus",f0->GetParameter(1)-2.*f0->GetParameter(2),f0->GetParameter(1)+2.*f0->GetParameter(2));
	    ped_sca.at(ichip).at(ichn).at(isca)->Fit("f1","RQME");
	    //  if(f1->GetParameter(1)>150 && f1->GetParameter(2)>2. ) {
	      ped_mean.at(ichip).at(ichn).push_back(f1->GetParameter(1));
	      ped_rms.at(ichip).at(ichn).push_back(f1->GetParameter(2));
	      fout_ped<<f1->GetParameter(1) << " " << f1->GetParError(1)<<" ";
	      pedestal_chip.at(ichip)->Fill(f1->GetParameter(1));

	      pedestal_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(1));
	      pedestal_width_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(2));
	      pedestal_error_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParError(1));
	      pedestal_chi2ndf_map[isca] -> Fill( map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetChisquare() / f1->GetNDF());
	      
	      /*  } else {
	      ped_mean.at(ichip).at(ichn).push_back(0);
	      ped_rms.at(ichip).at(ichn).push_back(0);
      	      fout_ped<<0<< " " << 0<<" ";
	      }*/
	  } else {
	    ped_mean.at(ichip).at(ichn).push_back(0);
	    ped_rms.at(ichip).at(ichn).push_back(0);
	    fout_ped<<0<< " " << 0<<" ";
	  }
	} else {
	  ped_mean.at(ichip).at(ichn).push_back(0);
          ped_rms.at(ichip).at(ichn).push_back(0);
	  fout_ped<<0<< " " << 0<<" ";
	}
      }
      fout_ped<<endl;
    }
  }

  pedfile->Close();

  TFile *pedfile_summary = new TFile("Pedestal_summary_"+dif+".root" , "RECREATE");
  pedfile_summary->cd();
  
  TCanvas *canvas_map = new TCanvas("pedestal_map_"+dif, "pedestal_map_"+dif,1200,1200);
  canvas_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_map->cd(isca+1);
    pedestal_map[isca]->SetStats(kFALSE);
    pedestal_map[isca]->SetTitle("pedestal_map, "+dif0);
    pedestal_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_map[isca]->GetYaxis()->SetTitle("y");
    pedestal_map[isca]->GetZaxis()->SetRangeUser(200,400);

    //noise.at(ichip)->SetLineColor(2);
    pedestal_map[isca]->Draw("colz");
    pedestal_map[isca]->Write();
   }

  TCanvas *canvas_width_map = new TCanvas("pedestal_width_map_"+dif, "pedestal_width_map_"+dif,1200,1200);
  canvas_width_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_width_map->cd(isca+1);
    pedestal_width_map[isca]->SetStats(kFALSE);
    pedestal_width_map[isca]->SetTitle("pedestal_width_map, "+dif0);
    pedestal_width_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_width_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_width_map[isca]->Draw("colz");
    pedestal_width_map[isca]->Write();
   }

  TCanvas *canvas_error_map = new TCanvas("pedestal_error_map_"+dif, "pedestal_error_map_"+dif,1200,1200);
  canvas_error_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_error_map->cd(isca+1);
    pedestal_error_map[isca]->SetStats(kFALSE);
    pedestal_error_map[isca]->SetTitle("pedestal_error_map, "+dif0);
    pedestal_error_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_error_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_error_map[isca]->Draw("colz");
    pedestal_error_map[isca]->Write();
   }

  TCanvas *canvas_npeaks_map = new TCanvas("pedestal_npeaks_map_"+dif, "pedestal_npeaks_map_"+dif,1200,1200);
  canvas_npeaks_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_npeaks_map->cd(isca+1);
    pedestal_npeaks_map[isca]->SetStats(kFALSE);
    pedestal_npeaks_map[isca]->SetTitle("pedestal_npeaks_map, "+dif0);
    pedestal_npeaks_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_npeaks_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_npeaks_map[isca]->Draw("colz");
    pedestal_npeaks_map[isca]->Write();
   }

  TCanvas *canvas_chi2ndf_map = new TCanvas("pedestal_chi2ndf_map_"+dif, "pedestal_chi2ndf_map_"+dif,1200,1200);
  canvas_chi2ndf_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_chi2ndf_map->cd(isca+1);
    pedestal_chi2ndf_map[isca]->SetStats(kFALSE);
    pedestal_chi2ndf_map[isca]->SetTitle("pedestal_chi2ndf_map, "+dif0);
    pedestal_chi2ndf_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_chi2ndf_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chi2ndf_map[isca]->Draw("colz");
    pedestal_chi2ndf_map[isca]->Write();
   }

  TCanvas *canvas_entries_map = new TCanvas("pedestal_entries_map_"+dif, "pedestal_entries_map_"+dif,1200,1200);
  canvas_entries_map->Divide(4,4);
  for(int isca=0; isca<15; isca ++) {
    TString dif0=dif+TString::Format("_sca%i",isca);
    canvas_entries_map->cd(isca+1);
    pedestal_entries_map[isca]->SetStats(kFALSE);
    pedestal_entries_map[isca]->SetTitle("pedestal_entries_map, "+dif0);
    pedestal_entries_map[isca]->GetXaxis()->SetTitle("x");
    pedestal_entries_map[isca]->GetYaxis()->SetTitle("y");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_entries_map[isca]->Draw("colz");
    pedestal_entries_map[isca]->Write();
   }
  

  canvas_map->Write();
  canvas_width_map->Write();
  canvas_error_map->Write();
  canvas_npeaks_map->Write();
  canvas_entries_map->Write();
  canvas_chi2ndf_map->Write();

   
  TCanvas *canvas_pedestal = new TCanvas("pedestal_average_"+dif, "pedestal_average_"+dif,1200,1200);
  canvas_pedestal->Divide(4,4);
  for(int ichip=0; ichip<16; ichip++) {
    canvas_pedestal->cd(ichip+1);
    //gPad->SetLogy();
    pedestal_chip.at(ichip)->GetXaxis()->SetRangeUser(0,500);
    pedestal_chip.at(ichip)->SetTitle(TString::Format("Average pedestal, chip-%i",ichip));
    pedestal_chip.at(ichip)->GetXaxis()->SetTitle("ADC");
    pedestal_chip.at(ichip)->GetYaxis()->SetTitle("#");
    //noise.at(ichip)->SetLineColor(2);
    pedestal_chip.at(ichip)->Draw("hs");
    //pedestal_chip.at(ichip)->Write();
  }

  canvas_pedestal->Write();
  pedfile_summary->Close();



}
