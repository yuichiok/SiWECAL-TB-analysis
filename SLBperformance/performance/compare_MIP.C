#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include <iostream>
#include <fstream>

// #include "../../style/Style.C"
// #include "../../style/Labels.C"

using namespace std;

void compare_MIP()
{ // TString filename, int slabadd){

  TString filename = "plots/SLB_3_Run_ILC_slab14_HV100_12192019_15h07min_Ascii.root";

  TFile *_file_1 = TFile::Open(filename);

  TH1F *mip_1[16][64];
  for (int iasic = 0; iasic < 16; iasic++)
  {
    for (int ichn = 0; ichn < 64; ichn++)
    {

      if (_file_1->Get(TString::Format("histograms_mip/mip_chip%i_chn%i", iasic, ichn)))
        mip_1[iasic][ichn] = (TH1F *)_file_1->Get(TString::Format("histograms_mip/mip_chip%i_chn%i", iasic, ichn));
      else
        mip_1[iasic][ichn] = NULL;
    }
  }

  filename = "plots/SLB_3_Run_ILC_slab14_HV180_01022020_15h10min_Ascii.root";

  TFile *_file_2 = TFile::Open(filename);

  TH1F *mip_2[16][64];
  for (int iasic = 0; iasic < 16; iasic++)
  {
    for (int ichn = 0; ichn < 64; ichn++)
    {

      if (_file_2->Get(TString::Format("histograms_mip/mip_chip%i_chn%i", iasic, ichn)))
        mip_2[iasic][ichn] = (TH1F *)_file_2->Get(TString::Format("histograms_mip/mip_chip%i_chn%i", iasic, ichn));
      else
        mip_2[iasic][ichn] = NULL;
    }
  }

  TCanvas *canvas_asic[16];

  for (int iasic = 5; iasic < 6; iasic++)
  {

    canvas_asic[iasic] = new TCanvas(TString::Format("asic_%i", iasic), TString::Format("asic_%i", iasic), 1600, 1000);
    canvas_asic[iasic]->Divide(8, 8);

    TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    bool write_leg1 = false;
    bool write_leg2 = false;

    for (int ichn = 0; ichn < 64; ichn++)
    {
      if (mip_1[iasic][ichn] == NULL)
        continue;

      canvas_asic[iasic]->cd(1 + ichn);
      mip_1[iasic][ichn]->GetXaxis()->SetTitle("DAC");
      mip_1[iasic][ichn]->GetYaxis()->SetTitle("Nhits");
      mip_1[iasic][ichn]->SetLineColor(1);
      mip_1[iasic][ichn]->SetTitle(TString::Format("ASIC %i", iasic));
      if (write_leg1 == false)
      {
        leg2->AddEntry(mip_1[iasic][ichn], "HV=100V", "l");
        write_leg1 = true;
      }
      mip_1[iasic][ichn]->DrawNormalized("histo");

      if (mip_2[iasic][ichn] != NULL)
      {
        mip_2[iasic][ichn]->SetLineColor(2);
        mip_2[iasic][ichn]->DrawNormalized("histosame");
        if (write_leg2 == false)
        {
          leg2->AddEntry(mip_2[iasic][ichn], "HV=180V", "l");
          write_leg2 = true;
        }
      }
      leg2->Draw();
    }
  }
}
