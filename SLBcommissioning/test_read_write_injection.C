//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

void test_read_write_injection(TString filename="Run_Settings.txt", bool debug=true) {

  for(int irow=0; irow<8; irow++) {

    read_configuration_file("06102020/Run_Settings.txt",false);
    for(int islab=0; islab<15; islab++) {
      if(islab==8 || islab==12) {
	for (int ichip=0; ichip<4; ichip++) {
	  mask_full_chip(0,islab,0,ichip);
	  enable_trig_row(0,islab,0,ichip,irow);
	}
      } else {
	for (int ichip=0; ichip<16; ichip++) {
	  mask_full_chip(0,islab,0,ichip);
	  enable_trig_row(0,islab,0,ichip,irow);
	}
      }
    }
    write_configuration_file(TString::Format("06102020/Run_Settings_InjEnabled_row%i.txt",irow));
  }
}
