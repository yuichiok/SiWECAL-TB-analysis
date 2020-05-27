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
#include "TF1.h"
#include "../scurves/savehistos.C"
#include "../scurves/fithistos.C"

void scurves(TString date, int iteration) {

  savehistos(date);
  std::vector<int> nslboards;
  nslboards.push_back(3);
  nslboards.push_back(12);
  //  nslboards.push_back(13);
  
  fithistos(date,nslboards,iteration);
  
}
