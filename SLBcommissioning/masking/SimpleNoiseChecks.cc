//# Copyright 2020  Adrián Irles (IJCLab, CNRS/IN2P3)

#include "SimpleNoiseChecks.h"

int SimpleNoiseChecks(TString filename_in="05182020", TString round="first", int iteration=0){
  //, TString conf_init="../Run_Settings_it0.txt" , TString conf_out="../Run_Settings_it2.txt"){

 
 
  init(filename_in,true);
  
  //each channel number of riggers, udnerflows, retriggers, etc is compared with the "expected number of cosmics" expected, via a threshold definded in triple check
  
  int voting=3;

  if(round=="masking") {
    if(iteration<3)   {
      cout<<" ------------------------------------------------------------------------"<<endl;
      cout<<" ITERATION = "<<iteration<<endl;
      cout<<" FIRST Round of MASKING runs (before threshold DAC optimization). File to analyze: "<<filename_in<<endl;
      cout<<" We take runs with very high trhesholds and set noise level definitions very low (up to zero). "<<endl;
      cout<<" In addition, we also look at channels that have pedestals with ADC ~ 0 (or very large) "<<endl;
    } 
    // Sets of short runs of 1min with acq of 1ms and 10 Hz, DAC=350,
    // We check that the noisy channels by a triple voting (three consecutive runs of the type:
    //          Run_ILC_06312020_masking_it0_0_Ascii.dat
    //          Run_ILC_06312020_masking_it0_1_Ascii.dat
    //          Run_ILC_06312020_masking_it0_2_Ascii.dat!
    filename_in=filename_in+"_"+round+"_it"+iteration;
    if(iteration==0) triple_check(filename_in,iteration,voting,acqwindow, 0., 0., 0.5, 1.);// underflow trig, underflow trig, retrig, trig 
    if(iteration==1) triple_check(filename_in,iteration,voting,acqwindow, 0., 0., 0.5, 1.);//
    if(iteration==2) triple_check(filename_in,iteration,voting,acqwindow, 0., 0., 1.,  2.);//
    
    
    // Trigger =350/325/300/275, several iterations (after preliminary masking). 
    // 1ms and 10 Hz
    // 1 min 
    // We check that the noisy channels by a triple voting (three consecutive runs of the type:
    //          Run_ILC_06312020_masking_it3_0_Ascii.dat
    //          Run_ILC_06312020_masking_it3_1_Ascii.dat
    //          Run_ILC_06312020_masking_it3_2_Ascii.dat!
    if(iteration>2)   {
      cout<<" ------------------------------------------------------------------------"<<endl;
      cout<<" ITERATION = "<<iteration<<endl;
      cout<<" SECOND Round of MASKING runs: analyze file: "<<filename_in<<endl;
      cout<<" We take runs with very high trhesholds and iteratively decrease it at the same time that we relax our noise level definitions. "<<endl;
      cout<<" We only look at triggered cells this time."<<endl;
    }
    if(iteration==3) triple_check(filename_in,iteration,voting,acqwindow, 1.,9999.,1.,1.); // DAC=350
    if(iteration==4) triple_check(filename_in,iteration,voting,acqwindow, 1.,9999.,1.,1.); // DAC=325
    if(iteration==5) triple_check(filename_in,iteration,voting,acqwindow, 1.,9999.,1.,1.); // DAC=300
    if(iteration==6) triple_check(filename_in,iteration,voting,acqwindow, 1.,9999.,2.,5.); // DAC=275
    if(iteration==7) triple_check(filename_in,iteration,voting,acqwindow, 1.,9999.,1.,2.5); // DAC=275
    if(iteration>7 && iteration<11) triple_check(filename_in,iteration,voting,acqwindow, 1.,9999.,1.,2.5); // DAC=optimzied
    if(iteration>10) triple_check(filename_in,iteration,voting,acqwindow, 5.,9999.,5.,5); // DAC=optimzied                                                                                                 



  }//
  
 
  
  //-----------------------------------------------------------------------------
  // ATTENTION: HAVE YOU PERFORMED THE ITERATION 8 ? It consists on a quick run with optimized threholds
  // ex: Run_ILC_06312020_masking_it8_Ascii.dat
  
  // -------------- OPTIMIZATION OF THE MASKING WITH COSMICs
  // few, long iterations for cosmics (100ms)
  if(round=="cosmic") {
    cout<<" Cosmic Round of masking runs: analyze file: "<<filename_in<<"  iteration="<<iteration<<endl;
    if(iteration<13) cosmics_check(filename_in,iteration, acqwindow_cosm,5,true); // first iterations (9-12) are for threshold adjusting... If too many channels in one chip are going to be masked (>25% of the available ones), the global chip trehshold in increased by 5, instead of masking the channels
    else cosmics_check(filename_in,iteration, acqwindow_cosm,5,true);
    //run few,times before launching a long cosmics...
    
  }


  return 0;
}



