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

#include "../conf_struct.h"

#include "../scurves/savehistos.C"
#include "../scurves/fithistos.C"

void scurves(TString date, int iteration) {

  savehistos(date);
  std::vector<int> nslboards;
  // stack
  // for(int i=0; i<15; i++)  nslboards.push_back(i);
  
  // test bench
  nslboards.push_back(0);
  
  //first  found slboard
  // nslboards.push_back(1);
  // nslboards.push_back(2);
  // nslboards.push_back(3);
  // nslboards.push_back(4);
  // nslboards.push_back(5);
  // nslboards.push_back(6);
  // nslboards.push_back(7);
  // nslboards.push_back(10);
  // nslboards.push_back(11);
  // nslboards.push_back(13);
  // nslboards.push_back(14);

  //    nslboards.push_back(3);
  //nslboards.push_back(12);
  //  nslboards.push_back(13);

  // nslboards.push_back(3);



  // debug (TURN ON LATER)
  fithistos(date,nslboards,iteration);
  
}
