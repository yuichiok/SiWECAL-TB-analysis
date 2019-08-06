#define Layers_compar_cxx
#include "Layers_compar.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Layers_compar::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Layers_compar.C
//      root> Layers_compar t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
  TFile *newFile= new TFile("run_32015_histograms.root", "NEW");
  Long64_t nentries = fChain->GetEntriesFast();
  Int_t NSLAB=9;
  Int_t NSCA=15;
  Int_t NCHIP=4;
  Int_t var_bcid[NCHIP][NSLAB][NSLAB][NSCA][NSCA];
  TH1F *histo_coincidence_dif_slab[NCHIP][NSLAB][NSLAB];
  TH2F *histo_coincidence_dif_slab_2D[NCHIP][NSLAB][NSLAB];
  TCanvas* canvas1[NCHIP];
  TCanvas* canvas2[NCHIP];
  // TCanvas* canvas3;
  // TCanvas* canvas4;
  // TCanvas* canvas5;
  

  // canvas3=new TCanvas("run_32004_dif3_slab","coincidences_dif_3", 1200,600);
  // canvas3->Divide(2,2);
  // canvas4=new TCanvas("run_32004_dif4_slab","coincidences_dif_4", 1200,600);
  // canvas4->Divide(2,2);
  // canvas5=new TCanvas("run_32004_dif5_slab","coincidences_dif_5", 1200,600);
  // canvas5->Divide(2,2);
  Long64_t nbytes = 0, nb = 0;
  for(Int_t chip=0; chip<NCHIP; chip++) {
    Int_t chip0=chip;
    Int_t chip1=15-chip;
    canvas1[chip]=new TCanvas(TString::Format("run_32001_slab_chip_%i_chip_%i", chip0, chip1),TString::Format("slab_chip%i_chip_%i",chip0, chip1), 1800,1400);
  canvas1[chip]->Divide(9,9);
  
  canvas2[chip]=new TCanvas(TString::Format("run_32001_slab_2D_chip_%i_chip_%i", chip0, chip1),TString::Format("var_bcid_chip_%i_chip_%i", chip0, chip1), 1800,1400);
  
  canvas2[chip]->Divide(9,9);
    for (Int_t slab_init0=0; slab_init0<NSLAB; slab_init0++) {
      for (Int_t slab_init1=0; slab_init1<NSLAB; slab_init1++) {
	Int_t chip0=chip;
	Int_t chip1=chip;
	if (slab_init0<5 && slab_init1>4) {
	  Int_t chip0=15-chip;
	  Int_t chip1=chip;
	  histo_coincidence_dif_slab[chip][slab_init0][slab_init1]=new TH1F(TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), 8001, -4000.5, 4000.5);
	  histo_coincidence_dif_slab_2D[chip][slab_init0][slab_init1]=new TH2F(TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), 501,-5,5005,451,-2005,2505);
	}
	if (slab_init0>4 && slab_init1<5) {
	  Int_t chip0=chip;
	  Int_t chip1=15-chip;
	  histo_coincidence_dif_slab[chip][slab_init0][slab_init1]=new TH1F(TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), 8001, -4000.5, 4000.5);
	  histo_coincidence_dif_slab_2D[chip][slab_init0][slab_init1]=new TH2F(TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), 501,-5,5005,451,-2005,2505);
	}
	if (slab_init0<5 && slab_init1<5) {
	  Int_t chip0=15-chip;
	  Int_t chip1=15-chip;
	  histo_coincidence_dif_slab[chip][slab_init0][slab_init1]=new TH1F(TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), 8001, -4000.5, 4000.5);
	  histo_coincidence_dif_slab_2D[chip][slab_init0][slab_init1]=new TH2F(TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), 501,-5,5005,451,-2005,2505);
	}
	if (slab_init0>4 && slab_init1>4) {
	  Int_t chip0=chip;
	  Int_t chip1=chip;
	  histo_coincidence_dif_slab[chip][slab_init0][slab_init1]=new TH1F(TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), TString::Format("number_coincidences_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), 8001, -4000.5, 4000.5);
	  histo_coincidence_dif_slab_2D[chip][slab_init0][slab_init1]=new TH2F(TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0,slab_init1, chip1), TString::Format("var_bcid_slab_%i_chip_%i_slab_%i_chip_%i", slab_init0, chip0, slab_init1, chip1), 501,-5,5005,451,-2005,2505);
	}
      }
    }

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // if (Cut(ientry) < 0) continue;
      for (Int_t slab0=0; slab0<NSLAB; slab0++) {
	for (Int_t sca0=0; sca0<NSCA; sca0++) {
	  for (Int_t slab1=0; slab1<NSLAB; slab1++) {
	    for (Int_t sca1=0; sca1<NSCA; sca1++) {
	      var_bcid[chip][slab0][slab1][sca0][sca1]=10000;
	      Int_t chip0=chip;
	      Int_t chip1=chip;
	      if (slab0<5) {
		chip0=15-chip;
	      }
	      if (slab1<5) {
		chip1=15-chip;
	      }	
	      if (bcid[slab0][chip0][sca0]>50 && bcid[slab1][chip1][sca1]>50)
		if ((bcid[slab0][chip0][sca0]<890 || bcid[slab0][chip0][sca0]>930 ) &&  (bcid[slab1][chip1][sca1]<890 || bcid[slab1][chip1][sca1]>930 ) )
		  if( bcid[slab0][chip0][sca0]<4900 && bcid[slab1][chip1][sca1]<4900)
		    if ( (bcid[slab0][chip0][sca0]<2265 || bcid[slab0][chip0][sca0]>2285 ) &&  (bcid[slab1][chip1][sca1]<2265 || bcid[slab1][chip1][sca1]>2285 ) )
		      if ((bcid[slab0][chip0][sca0]<3095 || bcid[slab0][chip0][sca0]>3105 ) &&  (bcid[slab1][chip1][sca1]<3095 || bcid[slab1][chip1][sca1]>3105 ) ) {
			var_bcid[chip][slab0][slab1][sca0][sca1]=bcid[slab0][chip0][sca0]-bcid[slab1][chip1][sca1];
			//if (var_bcid[dif0][slab0][sca0][sca1]<=-10 || var_bcid[dif0][slab0][sca0][sca1]>=10) {
			histo_coincidence_dif_slab[chip][slab0][slab1]->Fill(var_bcid[chip][slab0][slab1][sca0][sca1]);
		      }
	      
	      if (bcid[slab0][chip0][sca0]>-1 && bcid[slab1][chip1][sca1]>-1 ) {
		var_bcid[chip][slab0][slab1][sca0][sca1]=bcid[slab0][chip0][sca0]-bcid[slab1][chip1][sca1];
		histo_coincidence_dif_slab_2D[chip][slab0][slab1]->Fill(bcid[slab1][chip1][sca1], var_bcid[chip][slab0][slab1][sca0][sca1]);
		//}
	      }
	    }
	  }
	}
      }
    }
  
    for (Int_t slab0=0; slab0<NSLAB; slab0++) {
      for (Int_t slab1=0; slab1<NSLAB; slab1++) {
      canvas1[chip]->cd(slab0*9+slab1+1);
      canvas1[chip]->SetLogy(1);
      gPad->SetLogy(1);
      histo_coincidence_dif_slab[chip][slab0][slab1]->SetLineColor(1);
      histo_coincidence_dif_slab[chip][slab0][slab1]->SetLineStyle(1);
      histo_coincidence_dif_slab[chip][slab0][slab1]->SetLineWidth(1);
      histo_coincidence_dif_slab[chip][slab0][slab1]->GetXaxis()->SetTitle("delta_bcid");
      histo_coincidence_dif_slab[chip][slab0][slab1]->GetYaxis()->SetTitle("entries");
      histo_coincidence_dif_slab[chip][slab0][slab1]->Draw();

      canvas2[chip]->cd(slab0*9+slab1+1);
      gPad->SetLogz(1);
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->SetMarkerColor(1);
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->SetMarkerStyle(1);
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->GetXaxis()->SetTitle("bcid_slb");
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->GetYaxis()->SetTitle("delta_bcid");
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->GetZaxis()->SetRangeUser(5,300);
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->GetYaxis()->SetRangeUser(-50,300);
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->Draw("colz");


      histo_coincidence_dif_slab[chip][slab0][slab1]->Write();
      histo_coincidence_dif_slab_2D[chip][slab0][slab1]->Write();
      
      }
    }
    canvas1[chip]->Write();
    canvas2[chip]->Write();
  }
    // canvas2->cd(slab0+1);
    // histo_coincidence_dif_slab[1][slab0]->SetLineColor(1);
    // histo_coincidence_dif_slab[1][slab0]->SetLineStyle(1);
    // histo_coincidence_dif_slab[1][slab0]->SetLineWidth(1);
    // histo_coincidence_dif_slab[1][slab0]->GetXaxis()->SetTitle("delta_bcid");
    // histo_coincidence_dif_slab[1][slab0]->GetYaxis()->SetTitle("coincidences dif_slab");
    // histo_coincidence_dif_slab[1][slab0]->Draw();
    // canvas3->cd(slab0+1);
    // histo_coincidence_dif_slab[2][slab0]->SetLineColor(1);
    // histo_coincidence_dif_slab[2][slab0]->SetLineStyle(1);
    // histo_coincidence_dif_slab[2][slab0]->SetLineWidth(1);
    // histo_coincidence_dif_slab[2][slab0]->GetXaxis()->SetTitle("delta_bcid");
    // histo_coincidence_dif_slab[2][slab0]->GetYaxis()->SetTitle("coincidences dif_slab");
    // histo_coincidence_dif_slab[2][slab0]->Draw();
    
    // canvas4->cd(slab0+1);
    // histo_coincidence_dif_slab[3][slab0]->SetLineColor(1);
    // histo_coincidence_dif_slab[3][slab0]->SetLineStyle(1);
    // histo_coincidence_dif_slab[3][slab0]->SetLineWidth(1);
    // histo_coincidence_dif_slab[3][slab0]->GetXaxis()->SetTitle("delta_bcid");
    // histo_coincidence_dif_slab[3][slab0]->GetYaxis()->SetTitle("coincidences dif_slab");
    // histo_coincidence_dif_slab[3][slab0]->Draw();
    
    // canvas5->cd(slab0+1);
    // histo_coincidence_dif_slab[4][slab0]->SetLineColor(1);
    // histo_coincidence_dif_slab[4][slab0]->SetLineStyle(1);
    // histo_coincidence_dif_slab[4][slab0]->SetLineWidth(1);
    // histo_coincidence_dif_slab[4][slab0]->GetXaxis()->SetTitle("delta_bcid");
    // histo_coincidence_dif_slab[4][slab0]->GetYaxis()->SetTitle("coincidences dif_slab");
    // histo_coincidence_dif_slab[4][slab0]->Draw();
      // }
  
  newFile->Close();
}
