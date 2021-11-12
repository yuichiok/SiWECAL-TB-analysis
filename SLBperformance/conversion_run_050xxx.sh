#!/bin/bash

run=$1
data_folder="/mnt/win1/Run_Data/"${run}"/"
output="../converter_SLB/convertedfiles/"${run}"/"

initial=${PWD}

cd ${data_folder}
FILE_start=${run}".dat"
if test -f "${FILE_start}"; then
    FILE_new=${FILE_start}"_0000"
    mv ${FILE_start} ${FILE_new}
fi
cd -

if [ -d "$output" ]; then
    cd $output
    ls -1tr | tail -n -1 | xargs -d '\n' rm -f --
    cd -
else
    mkdir $output
fi

cd ../converter_SLB

root -l -q ConvertDirectorySL_TB.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",2\) &
root -l -q ConvertDirectorySL_TB.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",1\) &
root -l -q ConvertDirectorySL_TB.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",0\)

cd -

cd ${initial}/../
# source compress_TB2021_endedrun.sh $run &
cd ${initial}

