#!/bin/bash


for run in  "Run_ILC_11052022_datarate_aq100" "Run_ILC_11052022_datarate_aq10"
do
    initial=${PWD}
    
#    data_folder="/home/airqui/cernbox/SiWECAL/TB2022/commissioning/"${run}"/"
    data_folder="/eos/project/s/siw-ecal/TB2022-06/commissioning/"${run}"/"
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
    
    root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",2\) &
    root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",1\) &
    root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",0\)
    
    cd -


    #sleep 20s
    hadd ../converter_SLB/convertedfiles/${run}.root  ../converter_SLB/convertedfiles/${run}/*root
    root -l -q Proto.cc\(\"../converter_SLB/convertedfiles/${run}\",\"${run}\",1\)  &
done
