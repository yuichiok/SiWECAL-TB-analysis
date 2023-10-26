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

// #include "../../style/Style.C"
// #include "../../style/Labels.C"

Double_t langaufun(Double_t *x, Double_t *par)
{
  // Fit parameters:
  // par[0]=Width (scale) parameter of Landau density
  // par[1]=Most Probable (MP, location) parameter of Landau density
  // par[2]=Total area (integral -inf to inf, normalization constant)
  // par[3]=Width (sigma) of convoluted Gaussian function
  //
  // In the Landau distribution (represented by the CERNLIB approximation),
  // the maximum is located at x=-0.22278298 with the location parameter=0.
  // This shift is corrected within this function, so that the actual
  // maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
  Double_t mpshift = -0.22278298;      // Landau maximum location

  // Control constants
  Double_t np = 100.0; // number of convolution steps
  Double_t sc = 5.0;   // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp - xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for (i = 1.0; i <= np / 2; i++)
  {
    xx = xlow + (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);

    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[4]  reasonable start values for the fit
  //   parlimitslo[4]  lower parameter limits
  //   parlimitshi[4]  upper parameter limits
  //   fitparams[4]    returns the final fit parameters
  //   fiterrors[4]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf

  Int_t i;
  Char_t FunName[100];

  sprintf(FunName, "Fitfcn_%s", his->GetName());

  TF1 *ffitold = (TF1 *)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold)
    delete ffitold;

  TF1 *ffit = new TF1(FunName, langaufun, fitrange[0], fitrange[1], 4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width", "MP", "Area", "GSigma");

  for (i = 0; i < 4; i++)
  {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  his->Fit(FunName, "RBQM"); // fit within specified range, use ParLimits, quiet, improve fit results (TMINUIT)

  ffit->GetParameters(fitparams); // obtain fit parameters
  for (i = 0; i < 4; i++)
  {
    fiterrors[i] = ffit->GetParError(i); // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare(); // obtain chi^2
  NDF[0] = ffit->GetNDF();          // obtain ndf

  return (ffit); // return fit function
}

void FitHoldscan_TB2021(TString rootfile = "MIPs_050031", int chip0 = 12)
{

  TFile *outFile = new TFile(TString("HoldScans_Run_") + rootfile(5, 11) + ".root", "RECREATE");
  outFile->mkdir("LangauFits");
  outFile->mkdir("HoldScans");

  int nSLB = 15;
  int holdv[10] = {30, 40, 50, 60, 70, 80, 100, 140};

  double holdsMPV[15][4][8];

  for (int ihold = 0; ihold < 8; ihold++)
  {

    TFile *_file0 = TFile::Open(TString::Format("%s_hold%i.root", rootfile.Data(), holdv[ihold]));
    for (int ilayer = 0; ilayer < nSLB; ilayer++)
    {

      for (int ichip = chip0; ichip < chip0 + 4; ichip++)
      {
        std::cout << "Processing hold: " << holdv[ihold] << " layer: " << ilayer << " chip: " << ichip << std::endl;
        TH1F *htemp = (TH1F *)_file0->Get(TString::Format("layer_%i/charge_layer%i_chip%i", ilayer, ilayer, ichip));

        Double_t fr[2];
        Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
        // Width MP Area GSigma
        pllo[0] = 1.0;
        pllo[1] = 20;
        pllo[2] = 100.0;
        pllo[3] = 1;
        plhi[0] = 100.0;
        plhi[1] = 100.0;
        plhi[2] = 100000000.0;
        plhi[3] = 30.0;
        sv[0] = 3.0;
        sv[1] = 60.0;
        sv[2] = 10000.0;
        sv[3] = 8.0;
        Double_t chisqr;
        Int_t ndf;

        if (htemp != nullptr && htemp->GetEntries() > 100)
        {
          fr[0] = 0;
          fr[1] = 150 + 100 * ilayer;
          TF1 *fitsnr_temp = langaufit(htemp, fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf);

          outFile->cd("LangauFits");
          htemp->SetTitle(htemp->GetTitle() + TString::Format("_hold%i", holdv[ihold]) + ";MP;Entries");
          htemp->SetName(htemp->GetTitle() + TString::Format("_hold%i", holdv[ihold]));
          htemp->Write();

          double mpv = fitsnr_temp->GetParameter(1);
          double empv = fitsnr_temp->GetParError(1);
          double wmpv = fitsnr_temp->GetParameter(0);
          double chi2ndf = 0;

          holdsMPV[ilayer][ichip % 4][ihold] = mpv;
        }
        else
        {
          holdsMPV[ilayer][ichip % 4][ihold] = -1;
        }
      }
    }

    _file0->Close();
  }

  outFile->cd("HoldScans");

  std::ofstream output;
  output.open("FittedHoldValues.txt", ios_base::trunc);
  output << "Layer \t Chip \t MaxHold" << std::endl;

  double holdsDoubles[8] = {30, 40, 50, 60, 70, 80, 100, 140};

  for (int iLayer = 0; iLayer < 15; iLayer++)
  {

    TCanvas *layerCanvas = new TCanvas(TString::Format("HoldsVsMPV_Layer%i", iLayer));
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (int iChip = 0; iChip < 4; iChip++)
    {

      double *y = holdsMPV[iLayer][iChip];

      TGraph *graph = new TGraph();

      int n = 0;
      for (int iHold = 0; iHold < 8; iHold++)
      {

        if (y[iHold] > 0)
        {
          graph->Set(n + 1);
          graph->SetPoint(n, holdsDoubles[iHold], y[iHold]);
          n++;
        }
      }

      graph->SetMarkerColor(iChip + 1);
      graph->SetTitle(TString::Format("MPVVsHold_Layer%i;Hold;MPV", iLayer));
      graph->GetYaxis()->SetRangeUser(0, 150);

      graph->Fit("pol2", "Q");

      TF1 *fitfunc = (TF1 *)graph->GetListOfFunctions()->FindObject("pol2");
      if (fitfunc != nullptr)
      {
        output << iLayer << "\t" << chip0 + iChip << "\t" << fitfunc->GetMaximumX(0, 140) << std::endl;
        fitfunc->SetLineColor(iChip + 1);
      }

      legend->AddEntry(graph, TString::Format("Chip%i", chip0 + iChip));

      TString drawOption = "AP*";
      if (iChip != 0)
        drawOption = "P*SAME";

      graph->Draw(drawOption);
    }

    legend->Draw();
    layerCanvas->Write();
    delete layerCanvas;
  }

  output.close();
  outFile->Close();
}
