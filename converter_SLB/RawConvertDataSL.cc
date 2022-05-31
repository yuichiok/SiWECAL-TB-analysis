#include "SLBraw2ROOT.cc"

using namespace std;

int RawConvertDataSL(TString filename, bool zerosupression=false, TString outputname="default", bool getbadbcid_bool=true){
    SLBraw2ROOT ss;
    ss._maxReadOutCycleJump=10;
    bool result=false;
    while(result==false) {
      result=ss.ReadFile(filename, true, outputname, zerosupression, getbadbcid_bool);
      ss._maxReadOutCycleJump*=10;
    }
    gSystem->Exit(0);
    return 0;
}
