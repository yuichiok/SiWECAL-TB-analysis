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

  // TCanvas *c_MIPSummary[15];

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
                TH1F* h_0 = (TH1F*)h->Clone();
                std::vector<double> result = resultfit(h_0,"high");

                h_0->Draw();

                TLatex latex;
                const char *longstring = "K_{S}... K^{*0}... #frac{2s}{#pi#alpha^{2}} #frac{d#sigma}{dcos#theta} (e^{+}e^{-} #rightarrow f#bar{f} ) = #left| #frac{1}{1 - #Delta#alpha} #right|^{2} (1+cos^{2}#theta)";
                latex.SetTextSize(0.025);
                latex.SetTextAlign(13);  //align at top
                latex.DrawLatex(.2,.9,"K_{S}");
                latex.DrawLatex(.3,.9,"K^{*0}");
                latex.DrawLatex(.2,.8,longstring);
                c_MIPSummary_layer7_0->Draw();


                drawIt = false;
              }

              // h->Draw();
          }
        }

    }

    c_MIPSummary->cd();
    // c_MIPSummary->Draw();
    c_MIPSummary->SaveAs(TString::Format("plots/MIPs/MIPSummary_layer%i.pdf", i));
  }

}

