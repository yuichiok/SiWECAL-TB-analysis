//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "SimpleNoiseChecks.h"

int SimpleNoiseChecks(TString filename_in="05182020", TString round="first", int iteration=0){
  //, TString conf_init="../Run_Settings_it0.txt" , TString conf_out="../Run_Settings_it2.txt"){

 
 
  init(filename_in,true);
  // --------------FIRST-------------------------------------
  // Two short runs of 1-2min with acq of 1ms and 10 Hz, DAC=350
  // We check that the noisy channels are noisy in both runs !
  
  if(round=="first1") {
    cout<<" First Round of masking runs: analyze file: "<<filename_in<<endl;
    double_check(filename_in,iteration,1, 1, 0, 0, 0);
  }

  if(round=="first2") {
    cout<<" First Round of masking runs: analyze file: "<<filename_in<<endl;
    double_check(filename_in,iteration,2,1, 0., 0.5,1);
  }

  if(round=="first3") {
    cout<<" First Round of masking runs: analyze file: "<<filename_in<<endl;
    double_check(filename_in,iteration,3,1, 0., 1,2);
  } 
  
  //-----------------------------------------------------------------------------
  
  // --------------SECOND round of iterations-------------------------------------
  // Trigger =350/325/300/2750, several iterations (after preliminary masking). Not search for coincidences between masking
  // 1ms and 10 Hz
  // 1-2 min (or a fixed number of cycles of at least 600 cycles)
  if(round=="second") {
    cout<<" Second Round of masking runs: analyze file: "<<filename_in<<"  iteration="<<iteration<<endl;
    if(iteration<3 || iteration >7)  {
      cout<<" EXIT --> Tell me at what iteration of round "<<round<<" are you (should be between 3, 7, both included)"<<endl;
      return -1;
    }
    if(iteration==3) single_check(filename_in,iteration,1, 1,1,1); // DAC=350
    if(iteration==4) single_check(filename_in,iteration,1, 1,1,1); // DAC=325
    if(iteration==5) single_check(filename_in,iteration,1, 1,1,1); // DAC=300
    if(iteration==6) single_check(filename_in,iteration,1, 1, 2, 5);//DAC=275
    if(iteration==7) single_check(filename_in,iteration,1, 1, 1, 2.5);//DAC=275
  }//
  
 
  
  //-----------------------------------------------------------------------------

  // --------------THIRD round of iterations-------------------------------------
  // Primero elegimos los thresholds
  // 3 long iterations for cosmics (150ms)
  if(round=="cosmic") {
    cout<<" Cosmic Round of masking runs: analyze file: "<<filename_in<<"  iteration="<<iteration<<endl;
    if(iteration<13) cosmics_check(filename_in,iteration, 100,3); // first iterations (9,10,11,12) are for fine threshold adjusting... I am assuming here that all slabs are equipped with SK2a
    else cosmics_check(filename_in,iteration, 100,5);
    
  }


  return 0;
}



