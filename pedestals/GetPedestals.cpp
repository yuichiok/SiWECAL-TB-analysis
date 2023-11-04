
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include "../macros/Styles.hh"

template <typename T>
TCanvas* PlotPedestals(T hist, TString option)
{
  TString canvas_name = "c_" + TString(hist->GetName());
  TString pad_name = "pad_" + TString(hist->GetName());
  TCanvas* canvas = new TCanvas(canvas_name,canvas_name, 800, 800);
  TPad *pad = new TPad(pad_name,pad_name,0,0,1,1);
  pad->SetGrid(0,0);
  StylePad(pad, 0.1, 0.15, 0.17, 0.15);
  hist->GetZaxis()->SetTitle(hist->GetTitle());
  hist->GetZaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->SetTitle(";Layer #times 20 + Chip;SCA #times 100 + Channel");
  hist->Draw(option);
  return canvas;
}

void GetPedestals(const char* filename) {
  // Open the ROOT file
  TFile* file = TFile::Open(filename);
  if (!file) {
    std::cerr << "Error: failed to open file " << filename << std::endl;
    return;
  }

  // Extract the histograms
  TH2F* ped_all = (TH2F*) file->Get("ped_all");
  TH2F* width_all = dynamic_cast<TH2F*>(file->Get("width_all"));
  if (!ped_all || !width_all) {
    std::cerr << "Error: failed to extract histograms from file " << filename << std::endl;
    file->Close();
    return;
  }

  TCanvas *c_ped_all   = PlotPedestals(ped_all, "colz");
  TCanvas *c_width_all = PlotPedestals(width_all, "colz");

  // c_ped_all->Draw();
  // c_width_all->Draw();

  c_ped_all->SaveAs("plots/ped_all.pdf");
  c_width_all->SaveAs("plots/ped_width_all.pdf");

}
