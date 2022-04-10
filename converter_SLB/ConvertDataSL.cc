#include "SLBdecoded2ROOT.cc"

using namespace std;

int ConvertDataSL(TString filename, bool zerosupression=false, TString outputname="default"){
    SLBdecoded2ROOT ss;
    ss.ReadFile(filename, true, outputname, zerosupression);
    gSystem->Exit(0);
    return 0;
}
