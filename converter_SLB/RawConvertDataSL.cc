#include "SLBraw2ROOT.cc"

using namespace std;

int RawConvertDataSL(TString filename, bool zerosupression=false, TString outputname="default", bool getbadbcid_bool=false){
    SLBraw2ROOT ss;
    ss.ReadFile(filename, true, outputname, zerosupression, getbadbcid_bool);
    gSystem->Exit(0);
    return 0;
}
