//#include "MeanRMS_FEV10.cc"
#include "RAW2ROOT_TDC_test.cc"

using namespace std;

int ConvertDataTDC(TString filename){
    RAW2ROOT_TDC_test ss;
    ss.ReadFile(filename, true);//, 10, 4);
    /*    if (gener_rms) {
      TString filename2=filename;
      filename2+=".root";
      MeanRMS_FEV10 MRMS=MeanRMS_FEV10(filename2);
      }*/
    gSystem->Exit(0);
    return 0;
}
