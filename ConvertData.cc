//#include "MeanRMS_FEV10.cc"
#include "RAW2ROOT.cc"

using namespace std;

int ConvertData(TString filename){
    RAW2ROOT ss;
    ss.ReadFile(filename, true);//, 10, 4);
    /*    if (gener_rms) {
      TString filename2=filename;
      filename2+=".root";
      MeanRMS_FEV10 MRMS=MeanRMS_FEV10(filename2);
      }*/
    return 0;
}
