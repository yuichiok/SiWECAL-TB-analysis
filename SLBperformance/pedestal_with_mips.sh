#!/bin/bash

run_counter=0
initial=${PWD}
cd ../converter_SLB/convertedfiles/
for run in MIPScan_1.2pF_0.8pF_*merged.root
do

    cd ${initial}
    echo " ------- " $run
    run2=${run::-12}
    echo " ------- " $run2
    root -l -q NoiseMips.cc\(\"../converter_SLB/convertedfiles/${run2}_merged.root\",\"from${run2}\",1\) &
    root -l -q NoiseMips.cc\(\"../converter_SLB/convertedfiles/${run2}_merged.root\",\"from${run2}\",0\) 
    cd -
done

cd $initial
cd results_noise
for sca in {0..14}
do
    hadd -f NoiseCovariance_fromMIPScan_LowEnergyElectrons_highgain_sca${sca}.root NoiseCovariance_from*_highgain_sca${sca}.root &
    hadd -f NoiseCovariance_fromMIPScan_LowEnergyElectrons_lowgain_sca${sca}.root NoiseCovariance_from*_lowgain_sca${sca}.root
done

cd $initial

cd noise_covariance_matrix
source analysis.sh

#/mnt/HardDrive/cernbox/SiWECAL/TB2022/SiWECAL-TB-analysis/SLBperformance/noise_covariance_matrix
#cd ../converter_SLB/convertedfiles/
#hadd MIPScan_1.2pF_0.8pF.root MIPScan_1.2pF_0.8pF*.root

cd $initial


#cd TBchecks
#source analysis.sh
