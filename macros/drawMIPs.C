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
#include "Styles.hh"
#include "../SLBperformance/TBchecks/analysis.h"

// #include "../SLBco/conf_struct.h"

// #include "../../style/Style.C"
// #include "../../style/Labels.C"

using namespace std;

int nslabs = 15;
// int nchips = 16;
int nchannels = 64;
int nchannels_fit = 6;

Float_t charge_low[16][64];  // daughter, slab, asic, chn, DAC
Float_t charge_high[16][64]; // daughter, slab, asic, chn, DAC

void drawMIPs()
{

  //  TCanvas *Tlva = new TCanvas("Tlva","Tlva",500,500);
  //  Tlva->SetGrid();
  //  Tlva->DrawFrame(0,0,1,1);
  //  const char *longstring = "K_{S}... K^{*0}... #frac{2s}{#pi#alpha^{2}} #frac{d#sigma}{dcos#theta} (e^{+}e^{-} #rightarrow f#bar{f} ) = #left| #frac{1}{1 - #Delta#alpha} #right|^{2} (1+cos^{2}#theta)";
 
  //  TLatex latex;
  //  latex.SetTextSize(0.025);
  //  latex.SetTextAlign(13);  //align at top
  //  latex.DrawLatex(.2,.9,"K_{S}");
  //  latex.DrawLatex(.3,.9,"K^{*0}");
  //  latex.DrawLatex(.2,.8,longstring);
 
  //  Tlva->Draw();







  TFile *file_MIPSummary = new TFile("../mip_calib/MIPSummary_pedestalsubmode1_raw_siwecal_90021to90070_highgain.root", "READ");

  TCanvas * mip_all = new TCanvas("mip_all","mip_all",800,800);
  mip_all->SetLeftMargin(0.12);
  mip_all->SetRightMargin(0.15);
  TH2F * h_mip_all = (TH2F*) file_MIPSummary->Get("mpv_all");
  h_mip_all->SetTitle("");
  h_mip_all->GetXaxis()->SetTitle("Layer #times 20 + Chip");
  h_mip_all->GetYaxis()->SetTitle("Channel");
  h_mip_all->GetZaxis()->SetTitle("MIP value [ADC]");
  h_mip_all->GetZaxis()->SetTitleOffset(1.2);
  h_mip_all->GetZaxis()->SetRangeUser(0, 50);
  h_mip_all->Draw("colz");

  TCanvas * width_all = new TCanvas("width_all","width_all",800,800);
  width_all->SetLeftMargin(0.12);
  width_all->SetRightMargin(0.15);
  TH2F * h_width_all = (TH2F*) file_MIPSummary->Get("width_all");
  h_width_all->SetTitle("");
  h_width_all->GetXaxis()->SetTitle("Layer #times 20 + Chip");
  h_width_all->GetYaxis()->SetTitle("Channel");
  h_width_all->GetZaxis()->SetTitle("MIP width [ADC]");
  h_width_all->GetZaxis()->SetTitleOffset(1.2);
  h_width_all->GetZaxis()->SetRangeUser(0, 10);
  h_width_all->Draw("colz");



  Bool_t drawIt = true;
  TCanvas *c_MIPSummary_layer7_0 = new TCanvas(TString::Format("MIPSummary_layer%i_0", 7), TString::Format("MIPSummary_layer%i_0", 7), 800, 800);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);

  for (int i = 0; i < nslabs; i++)
  {
    TCanvas * c_MIPSummary = (TCanvas*) file_MIPSummary->Get(TString::Format("MIPSummary_layer%i", i));
    TList *primitives = c_MIPSummary->GetListOfPrimitives();
    TIter next(primitives);
    TObject* obj;
    while ((obj = next())) {
        TPad *pad = (TPad*)obj;
        // StylePad(pad,0,0.15,0,0.17);
        pad->SetLeftMargin(0.15);
        pad->SetGrid(1,1);

        TList *primitives2 = ((TPad*)obj)->GetListOfPrimitives();
        TIter next2(primitives2);
        TObject* obj2;
        while ((obj2 = next2())) {
          if (obj2->InheritsFrom("TH1")) {
              TH1F* h = (TH1F*)obj2;
              h->GetXaxis()->SetTitle("Charge [ADC]");
              h->GetYaxis()->SetTitle("Entries");

              if(i==7 && drawIt){
                c_MIPSummary_layer7_0->cd();
                c_MIPSummary_layer7_0->DrawFrame(0,0,1,1);
                TH1F* h_0 = (TH1F*)h->Clone();
                c_MIPSummary_layer7_0->SetGrid();
                c_MIPSummary_layer7_0->SetLeftMargin(0.15);
                std::vector<double> result = resultfit(h_0,"high");
                double MIP     = result.at(0);
                double MIP_err = result.at(1);
                double MIP_w   = result.at(2);
                double chi2ndf = result.at(3);
                double MIP_w_err = result.at(4);

                h_0->GetXaxis()->SetRangeUser(0, 100);
                h_0->Draw();

                TLatex latex;
                latex.SetTextSize(0.025);
                latex.SetTextAlign(13);  //align at top
                // latex.DrawLatex(65,5000,TString::Format("MIP = %.2f #pm %.2f\\\\
                // MIP width = %.2f #pm %.2f\\\\
                // #chi^{2}/ndf = %.2f",
                // MIP, MIP_err, MIP_w, MIP_w_err, chi2ndf));

                latex.DrawLatex(60,4700,TString::Format("MIP = %.2f #pm %.2f", MIP, MIP_err));
                latex.DrawLatex(60,4500,TString::Format("MIP width = %.2f #pm %.2f", MIP_w, MIP_w_err));
                latex.DrawLatex(60,4300,TString::Format("#chi^{2}/ndf = %.2f", chi2ndf));
                c_MIPSummary_layer7_0->Draw();


                drawIt = false;
              }

              // h->Draw();
          }
        }

    }

    c_MIPSummary->cd();
    // c_MIPSummary->Draw();
    // c_MIPSummary->SaveAs(TString::Format("plots/MIPs/MIPSummary_layer%i.pdf", i));
  }

}

