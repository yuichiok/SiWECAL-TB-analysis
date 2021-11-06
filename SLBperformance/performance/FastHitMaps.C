#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"

using namespace std;

void FastHitMaps(TString fileName = "None", TString outName = "HitMaps.root", bool PRplots=true, int zmax=6000){//TString filename, int slabadd){

  gStyle->SetOptStat(0);
  
  if(fileName == "None") {
    cout << "You need to give a filename" << endl;
    return;
  }

  TFile* inputFile = new TFile(fileName);
  if(inputFile == nullptr) {
    cout << "File not found" << endl;
    return;
  }

  TFile* outFile = new TFile(outName, "RECREATE");
  
  TCanvas* allLayersBeamSpot = new TCanvas("AllLayersBeamSpot");
  allLayersBeamSpot->Divide(4,4);

  TH2F* accumulatedBeamSpot = new TH2F("AccumulatedBeamSpot", "AccumulatedXY", 32, -90., 90., 32, -90., 90.);
  
  TH1F* entriesPerLayerHist = new TH1F("EntriesPerLayer", "EntriesPerLayer;Layer;Entries", 15, 0., 15.);
  
  for(int iLayer = 0; iLayer < 15; iLayer++) {

    TH2F* XYHits = (TH2F*)inputFile->Get(TString::Format("layer_%i/trig_xy_layer_%i",iLayer,iLayer));
    if(PRplots==true) {
      TCanvas* PRcanvas = new TCanvas("PRcanvas");
      PRcanvas->cd();
      XYHits->GetZaxis()->SetRangeUser(10,zmax);
      XYHits->SetTitle(TString::Format("Layer-%i",iLayer));
      XYHits->Draw("COLZ");
      if(iLayer<10) PRcanvas->Print(TString::Format("%s_hitmap_layer0%i.png",outName.Data(),iLayer));
      else  PRcanvas->Print(TString::Format("%s_hitmap_layer%i.png",outName.Data(),iLayer));
    }
    XYHits->Write();
    allLayersBeamSpot->cd(iLayer+1);
    XYHits->Draw("COLZ");


    entriesPerLayerHist->SetBinContent(iLayer, XYHits->GetEntries());

    for(int iX = 0; iX < 32; iX++) {

      for(int iY = 0; iY < 32; iY++) {

	accumulatedBeamSpot->SetBinContent(iX,iY, accumulatedBeamSpot->GetBinContent(iX,iY) + XYHits->GetBinContent(iX,iY));
	
      }
      
    }
    
  }

  

  accumulatedBeamSpot->Write();
  delete accumulatedBeamSpot;
  
  entriesPerLayerHist->Write();
  delete entriesPerLayerHist;
  
  allLayersBeamSpot->Write();
  delete allLayersBeamSpot;

  outFile->Close();
  
}





