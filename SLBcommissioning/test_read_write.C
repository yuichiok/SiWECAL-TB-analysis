#include "conf_struct.h"

void test_read_write(TString filename="Run_Settings.txt", bool debug=true) {

  read("Run_Settings.txt",false);
  write("Run_Settings_new.txt");
}
