#!/bin/bash


for run in "Pedestal_ch0_injection_run_050575" "Pedestal_ch20_injection_run_050576" "Pedestal_ch40_injection_run_050577" "Pedestal_ch60_injection_run_050578"
do
    initial=${PWD}
    
    data_folder="/mnt/HardDrive/beamData/Pedestals/"${run}"/"
    cd ${data_folder}
    FILE_start=${run}"_raw.bin"
    if test -f "${FILE_start}"; then
	FILE_new=${FILE_start}"_0000"
	mv ${FILE_start} ${FILE_new}
    fi
    
    cd -

    output="../converter_SLB/convertedfiles/"${run}
    if [ -d "$output" ]; then
	cd $output
	ls -1tr | tail -n -1 | xargs -d '\n' rm -f --
	cd -
    else
	mkdir $output
    fi
    
    cd ../converter_SLB
    
    #root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",2\) &
    #root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",1\) &
    #root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",0\)
    
    cd -
    
done

cd ../converter_SLB/convertedfiles/
#hadd -f Pedestal_run_050575_injection_merged.root Pedestal*/*root

cd -

cd noise_covariance_matrix

for run in "Pedestal_run_050575_injection_merged"
do
    cd ../
    #root -l -q Noise.cc\(\"../converter_SLB/convertedfiles/${run}.root\",\"${run}\",0\) &
    #root -l -q Noise.cc\(\"../converter_SLB/convertedfiles/${run}.root\",\"${run}\",1\)
    #sleep 1m
    cd -
    cd analysis
    #root -l -q NoiseStudy.C\(\"$run\",4\) &
    #root -l -q NoiseStudy.C\(\"$run\",2\) 
    #root -l -q NoiseStudy.C\(\"$run\",3\) &
    #root -l -q NoiseStudy.C\(\"$run\",1\) 
    root -l -q SummaryPlots.C\(\"$run\"\)
    cd -
done

cd -
