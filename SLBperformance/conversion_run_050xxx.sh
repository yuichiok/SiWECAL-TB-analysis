#!/bin/bash

run=$1
data_folder="/mnt/win1/Run_Data/run_"${run}"/"
output="../converter_SLB/convertedfiles/run_"${run}"/"

cd ${data_folder}
FILE_start="run_"${run}".dat"
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
root -l -q ConvertDirectorySL_TB.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/run_${run}\"\)
cd -


