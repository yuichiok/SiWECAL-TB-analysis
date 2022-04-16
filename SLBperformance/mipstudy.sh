#!/bin/bash

run_counter=0
initial=${PWD}
cd /mnt/HardDrive/beamData/mipscan/
for run in *
do
    
    run_counter=$((run_counter+1))
    
    data_folder="/mnt/HardDrive/beamData/mipscan/"${run}"/"
    cd ${data_folder}
    FILE_start=${run}"_raw.bin"
    if test -f "${FILE_start}"; then
	FILE_new=${FILE_start}"_0000"
	mv ${FILE_start} ${FILE_new}
    fi
    
    cd $initial

    output="../converter_SLB/convertedfiles/"MIPScan_1.2pF_0.8pF_${run_counter}
    if [ -d "$output" ]; then
	cd $output
	ls -1tr | tail -n -1 | xargs -d '\n' rm -f --
	cd -
    else
	mkdir $output
    fi
    
    cd ../converter_SLB
    
    root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/MIPScan_1.2pF_0.8pF_${run_counter}\",2\) &
    root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/MIPScan_1.2pF_0.8pF_${run_counter}\",1\) &
    root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/MIPScan_1.2pF_0.8pF_${run_counter}\",0\)

    sleep 2m
    
    if [ ${run_counter} -gt 9 ]
    then
	hadd convertedfiles/MIPScan_1.2pF_0.8pF_${run_counter}_merged.root convertedfiles/MIPScan_1.2pF_0.8pF_${run_counter}/*root
    else
	hadd convertedfiles/MIPScan_1.2pF_0.8pF_0${run_counter}_merged.root convertedfiles/MIPScan_1.2pF_0.8pF_${run_counter}/*root
    fi
    cd ../SLBperformance

    if [ ${run_counter} -gt 9 ]
    then
	root -l -q Proto.cc\(\"../converter_SLB/convertedfiles/MIPScan_1.2pF_0.8pF_${run_counter}_merged\",\"MIPScan_1.2pF_0.8pF_${run_counter}\",0\) &
    else
	root -l -q Proto.cc\(\"../converter_SLB/convertedfiles/MIPScan_1.2pF_0.8pF_0${run_counter}_merged\",\"MIPScan_1.2pF_0.8pF_0${run_counter}\",0\) &
    fi
    
    cd /mnt/HardDrive/beamData/mipscan/
    
done

cd $initial

cd ../converter_SLB/convertedfiles/
hadd MIPScan_1.2pF_0.8pF.root MIPScan_1.2pF_0.8pF*.root

cd $initial


cd TBchecks
source analysis.sh
