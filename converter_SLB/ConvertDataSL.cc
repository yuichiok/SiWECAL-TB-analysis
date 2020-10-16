//#include "MeanRMS_FEV10.cc"
#include "SLBdecoded2ROOT.cc"

using namespace std;

int ConvertDataSL(TString filename, int slboard, bool zerosupression=false){
    SLBdecoded2ROOT ss;
    ss.ReadFile(filename, true, "default", zerosupression);
    gSystem->Exit(0);
    return 0;
}
