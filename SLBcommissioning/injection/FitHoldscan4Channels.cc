
#include "TROOT.h"
#include "TFile.h"
#include "DecodedSLBAnalysis.cc"
#include "TGraphErrors.h"
#include "TLegend.h"

void FitHoldscan4Channels(TString run="20200226_dac1.15V_chn0to3", int nslboards=6, int channel1=0, int channel2=8, int channel3=16, int channel4=24){

  gROOT->SetBatch(kTRUE);
  
  cout<<" Holdscan file: " << run << endl;

  TGraphErrors* holdScan[15][16][4];

  TFile* input_file = new TFile("results/graphs_holdscan_" + run + ".root", "READ");
  input_file->cd();

  cout << input_file << endl;
  
  TH1F* slboardRelations = (TH1F*)input_file->Get("layer_slboard_relation");

  for(Int_t iLayer = 0; iLayer < nslboards; iLayer++)
    {
      Int_t slboard = slboardRelations->GetBinContent(iLayer+1);
      for(int j = 0; j < 16; j++) {
	for(int k = 0; k < 4; k++) {
	  int chan=0;
	  if(k==0) chan=channel1;
	  if(k==1) chan=channel2;
	  if(k==2) chan=channel3;
	  if(k==3) chan=channel4;
	  holdScan[iLayer][j][k] = (TGraphErrors*)input_file->Get(TString::Format("holdscan_layer%i_slboard%i_chip%i_chn%i",iLayer,slboard,j,chan));
	}
      }
    }
  
  input_file->Close();

  TFile* outPlotsFile = new TFile("results/Plots_holdscan_" + run + ".root", "RECREATE");
  
  //TGraphErrors* chn0_hold = new TGraphErrors(holdScan[0][0][channel1]->GetN(),)
  
  TFile* outputFile = new TFile("results/Fits_holdScan_" + run + ".root", "RECREATE");
  
  fstream hold_write;
  hold_write.open("results/hold_" + run + ".txt", fstream::out);

  TString header = "#hold scan results: " + run + "\n";
  header += "#layer chip hold error fitfunc\n";

  hold_write << header;
  
  for(Int_t i = 0; i < nslboards; i++)
    {
      outPlotsFile->cd();
      TCanvas* canvas_Allchips_holdscan = new TCanvas(TString::Format("canvas_layer%i_Allchips_holdscan",i),TString::Format("canvas_layer%i_Allchips_holdscan",i), 1200,1000);
      canvas_Allchips_holdscan->Divide(2,2);

      Double_t holds_Allchips[16] = {};
      Double_t holds_err_Allchips[16] = {};
      Double_t chipNr[16] = {};
      
      for(Int_t j= 0; j < 16; j++)
	{
	  chipNr[j] = j;
	  
	  outputFile->cd();
	  TCanvas* canvas = new TCanvas(TString::Format("canvas_layer%i_chip%i",i,j),TString::Format("canvas_layer%i_chip%i",i,j),1200,1000);
	  canvas->Divide(2,2);

	  vector<Double_t> holds;
	  Double_t mean_hold = 0.;
	  Double_t mean_hold_error = 0.;
	  TString fit_func_names;

	  TString option_Allchips = "lp";
	  if(j == 0) option_Allchips = "alp";
	  
	  for(Int_t ichan = 0; ichan < 4; ichan++)
	    {
	      int chan=0;
	      if(ichan==0) chan=channel1;
	      if(ichan==1) chan=channel2;
	      if(ichan==2) chan=channel3;
	      if(ichan==3) chan=channel4;

	      outPlotsFile->cd();
	      canvas_Allchips_holdscan->cd(ichan+1);

	      if(j<9)
		{
		  holdScan[i][j][ichan]->SetLineColor(j+1);
		  holdScan[i][j][ichan]->SetLineStyle(j+1);
		  holdScan[i][j][ichan]->SetMarkerColor(j+1);
		  holdScan[i][j][ichan]->SetMarkerStyle(j+1);
		}
	      else
		{
		  holdScan[i][j][ichan]->SetLineColor(j-9+2);
		  holdScan[i][j][ichan]->SetLineStyle(j-9+1);
		  holdScan[i][j][ichan]->SetMarkerColor(j-9+2);
		  holdScan[i][j][ichan]->SetMarkerStyle(j-9+1);		  
		  }

	      holdScan[i][j][ichan]->GetYaxis()->SetRangeUser(0.,500.);
	      holdScan[i][j][ichan]->SetTitle(TString::Format("Layer %i, channel: %i",i,chan) + ";Hold;Charge");
	      holdScan[i][j][ichan]->DrawClone(option_Allchips);
	      
	      outputFile->cd();
	      canvas->cd(ichan+1);
	      
	      TF1* best_fit_func = nullptr;
	      Double_t best_chi_ndf = 999.;
	      Int_t nParams;

	      TFitResultPtr res_pol3 = holdScan[i][j][ichan]->Fit("pol3","Q0S", "", 0., 220.);
	      TF1* pol3_func = new TF1(*(TF1*)holdScan[i][j][ichan]->GetListOfFunctions()->FindObject("pol3"));
	      Double_t chi_ndf_pol3 = pol3_func->GetChisquare()/pol3_func->GetNDF();
	      if(res_pol3->IsValid())
		{
		  best_fit_func = pol3_func;
		  best_chi_ndf = chi_ndf_pol3;
		  nParams = 4;
		}
	      
	      TFitResultPtr res_pol4 = holdScan[i][j][ichan]->Fit("pol4","Q0S", "", 0., 220.);
	      TF1* pol4_func = new TF1(*(TF1*)holdScan[i][j][ichan]->GetListOfFunctions()->FindObject("pol4"));	      
	      Double_t chi_ndf_pol4 = pol4_func->GetChisquare()/pol4_func->GetNDF();
	      if(chi_ndf_pol4 < best_chi_ndf && res_pol4->IsValid())
		{
		  best_fit_func = pol4_func;
		  best_chi_ndf = chi_ndf_pol4;
		  nParams = 5;
		}
	      
	      TFitResultPtr res_pol5 = holdScan[i][j][ichan]->Fit("pol5","Q0S", "", 0., 220.);
	      TF1* pol5_func = new TF1(*(TF1*)holdScan[i][j][ichan]->GetListOfFunctions()->FindObject("pol5"));
	      Double_t chi_ndf_pol5 = pol5_func->GetChisquare()/pol5_func->GetNDF();
	      if(chi_ndf_pol5 < best_chi_ndf && res_pol5->IsValid())
		{
		  best_fit_func = pol5_func;
		  best_chi_ndf = chi_ndf_pol5;
		  nParams = 6;
		}
	      
	      TFitResultPtr res_pol6 = holdScan[i][j][ichan]->Fit("pol6","Q0S", "", 0., 220.);
	      TF1* pol6_func = new TF1(*(TF1*)holdScan[i][j][ichan]->GetListOfFunctions()->FindObject("pol6"));
	      Double_t chi_ndf_pol6 = pol6_func->GetChisquare()/pol6_func->GetNDF();
	      if(chi_ndf_pol6 < best_chi_ndf && res_pol6->IsValid())
		{
		  best_fit_func = pol6_func;
		  best_chi_ndf = chi_ndf_pol6;
		  nParams = 7;
		}
	      
	      TLegend *leg = new TLegend(0.6,0.5,0.895, 0.895);
	      leg->AddEntry((TObject*)0, TString::Format("pol3 #chi^{2}/NDF: %.3f", chi_ndf_pol3));
	      leg->AddEntry((TObject*)0, TString::Format("pol4 #chi^{2}/NDF: %.3f", chi_ndf_pol4));
	      leg->AddEntry((TObject*)0, TString::Format("pol5 #chi^{2}/NDF: %.3f", chi_ndf_pol5));
	      leg->AddEntry((TObject*)0, TString::Format("pol6 #chi^{2}/NDF: %.3f", chi_ndf_pol6));
	      
	      holdScan[i][j][ichan]->DrawClone();

	      if(best_fit_func != nullptr)
		{
		  best_fit_func->DrawClone("same");
		  Double_t hold = best_fit_func->GetMaximumX(0., 220.);
		  
		  leg->AddEntry((TObject*)0, TString::Format("xMax: %.3f", hold));

		  mean_hold += hold;
		  holds.push_back(hold);
		  
		  fit_func_names += best_fit_func->GetName();
		  fit_func_names += " ";
		}

	      leg->Draw();

	      delete pol3_func;
	      delete pol4_func;
	      delete pol5_func;
	      delete pol6_func;
	    }

	  canvas->Write();
	  delete canvas;

	  mean_hold /= holds.size();
	  
	  for(Int_t iHold = 0; iHold < holds.size(); iHold++)
	    {
	      mean_hold_error += (holds[iHold] - mean_hold)*(holds[iHold] - mean_hold);
	    }

	  mean_hold_error /= holds.size();
	  mean_hold_error = TMath::Sqrt(mean_hold_error);
	  
	  hold_write << i << " " << j << " " << mean_hold << " " << mean_hold_error << " " << fit_func_names << "\n";

	  holds_Allchips[j] = mean_hold;
	  holds_err_Allchips[j] = mean_hold_error;
	  
	}

      outPlotsFile->cd();
      canvas_Allchips_holdscan->Write();
      delete canvas_Allchips_holdscan;

      TCanvas* canvas_Allchips_meanhold = new TCanvas(TString::Format("canvas_layer%i_Allchips_meanhold",i),TString::Format("canvas_layer%i_Allchips_meanhold",i), 1200,1000);

      TGraphErrors* meanhold_Allchips = new TGraphErrors(16,chipNr, holds_Allchips, 0, holds_err_Allchips);     
      meanhold_Allchips->GetYaxis()->SetRangeUser(0.,120.);
      meanhold_Allchips->SetTitle(TString::Format("Mean hold values. Layer %i",i) + ";ChipNr;Hold");

      meanhold_Allchips->Draw();

      canvas_Allchips_meanhold->Write();
      delete canvas_Allchips_meanhold;
      
    }
}

