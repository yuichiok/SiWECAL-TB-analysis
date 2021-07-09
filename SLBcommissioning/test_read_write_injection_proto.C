//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

void test_read_write_injection_proto(TString filename="Run_Settings.txt", bool debug=false) {

  for(int irow=0; irow<8; irow++) {

    read_configuration_file(filename+".txt",false);
    for(int islab=0; islab<15; islab++) {
      for (int ichip=0; ichip<16; ichip++) {
	disable_trig_otherrows(0,islab,0,ichip,irow);
      }
    }
    write_configuration_file(TString::Format("InjectionTest_%i.txt",irow));
  }
}
