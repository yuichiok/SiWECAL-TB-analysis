//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

void test_read_write(TString filename="Run_Settings.txt", bool debug=true) {

  read_configuration_file("Run_Settings.txt",false);
  write_configuration_file("Run_Settings_new.txt");
}
