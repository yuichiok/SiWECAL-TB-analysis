//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "SimpleNoiseChecks.h"

int SimpleNoiseChecks(TString filename_in="test", TString round="first", int iteration=0, bool thresholds=false){
  //, TString conf_init="../Run_Settings_it0.txt" , TString conf_out="../Run_Settings_it2.txt"){

 
 
  init();
  // --------------FIRST-------------------------------------
  // Two short runs of 2-3min with acq of 1ms and 5 Hz, DAC=350
  // Chech that the noisy channels are noisy in both runs !
  
  if(round=="first") {
    cout<<" First Round of masking runs: analyze file: "<<filename_in<<endl;
    double_check(filename_in,iteration,1, 5, 5);
  }
  
  //-----------------------------------------------------------------------------
  
  // --------------SECOND round of iterations-------------------------------------
  // Trigger =325/300, several iterations (after preliminary masking). Not search for coincidences between masking
  // 1ms and 5 Hz
  if(round=="second") {
    cout<<" Second Round of masking runs: analyze file: "<<filename_in<<"  iteration="<<iteration<<endl;
    if(iteration==0)  {
      cout<<" EXIT --> Tell me at what iteration of round "<<round<<" are you"<<endl;
      return -1;
    }
    if(iteration<3) single_check(filename_in,iteration,1, 2, 10, 10);
    else single_check(filename_in,iteration,3, 1, 5, 5);
  }//
  
 
  
  //-----------------------------------------------------------------------------

  // --------------THIRD round of iterations-------------------------------------
  // Primero elegimos los thresholds
  // 3 long iterations for cosmics (150ms)
  if(round=="cosmic") {
    cout<<" Cosmic Round of masking runs: analyze file: "<<filename_in<<"  iteration="<<iteration<<endl;
    cosmics_check(filename_in,iteration,1, 100, thresholds);
  }


  return 0;
}



