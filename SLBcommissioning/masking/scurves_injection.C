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

void scurves_injection(TString date="15102021", int iteration=9) {

  //savehistos(date+"_PROTO15_ADC85_diag1",date);
  // savehistos(date+"_PROTO15_ADC85_diag2",date);
  //savehistos(date+"_PROTO15_ADC85_diag3",date);
  // savehistos(date+"_PROTO15_ADC85_diag4",date);

  // savehistos(date+"_PROTO15_ADC120_diag1",date);
  //savehistos(date+"_PROTO15_ADC120_diag2",date);
  //savehistos(date+"_PROTO15_ADC120_diag3",date);
  //savehistos(date+"_PROTO15_ADC120_diag4",date);


  std::vector<int> nslboards;
  for(int i=0; i<15; i++)  nslboards.push_back(i);

  fithistos_injection(date,nslboards,iteration);
  
}
