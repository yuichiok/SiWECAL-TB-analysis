//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "conf_struct.h"

detector_t detector_new;

void inverse_slab();

void test_read_write_disablePA(TString filename="Run_Settings.txt", bool debug=true) {

  // read_configuration_file("Run_Settings.txt",false);
  read_configuration_file(filename,false);
  disable_PA_mask();
  // inverse_slab();
  write_configuration_file("Run_Settings_comm_it12.txt");
}

void inverse_slab(){

  for(int idaughter=0; idaughter < detector.n_core_daughters; idaughter++) {

    for(int islab=0; islab<detector.core_daughter_n_slabs[idaughter]; islab++) {

      int nslab = detector.core_daughter_n_slabs[idaughter] - 1;
      // std::cout << nslab - islab << std::endl;

      for (int iasu = 0; iasu < detector.slab[idaughter][islab].nb_asus; ++iasu) {

        for (int ichip = 0; ichip < detector.slab[idaughter][islab].asu[iasu].n_chips; ++ichip) {

          for(int ichn=0; ichn<detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].n_channels; ichn++) {

            if( detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[ichn]==1 ) {
              detector_new.slab[idaughter][nslab-islab].asu[iasu].skiroc[ichip].mask[ichn] = detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[ichn];
              detector_new.slab[idaughter][nslab-islab].asu[iasu].skiroc[ichip].preamplifier_mask[ichn] = detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].preamplifier_mask[ichn];
            }

          } // channel

        }   // chip

      }     // asu

    }       // slab

  }         // daughter

  // cout << "test out -> " << detector_new.slab[0][2].asu[0].skiroc[8].mask[3] << endl;


  for(int idaughter=0; idaughter < detector.n_core_daughters; idaughter++) {

    for(int islab=0; islab<detector.core_daughter_n_slabs[idaughter]; islab++) {

      for (int iasu = 0; iasu < detector.slab[idaughter][islab].nb_asus; ++iasu) {

        for (int ichip = 0; ichip < detector.slab[idaughter][islab].asu[iasu].n_chips; ++ichip) {

          for(int ichn=0; ichn<detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].n_channels; ichn++) {

            detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[ichn] = detector_new.slab[idaughter][islab].asu[iasu].skiroc[ichip].mask[ichn];
            detector.slab[idaughter][islab].asu[iasu].skiroc[ichip].preamplifier_mask[ichn] = detector_new.slab[idaughter][islab].asu[iasu].skiroc[ichip].preamplifier_mask[ichn];

          } // channel

        }   // chip

      }     // asu

    }       // slab

  }         // daughter


}