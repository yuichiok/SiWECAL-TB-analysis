#!/bin/bash

#run
run="run_050010"
run_file="converted.dat"
output=${PWD}"/../converter_SLB/convertedfiles/run_050010_07172021_13h52min_Ascii/"

initial_folder=$PWD


for i in {324..1774}
do
    j="000"$i
    if [ $i -gt 9 ]; then
        j="00"$i
    fi
    
    if [ $i -gt 99 ]; then
        j="0"$i
    fi

    if [ $i -gt 999 ]; then
        j=$i
    fi
    #analysis
    cd $initial_folder
    root -l -q Proto.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\",\"mip\",\"../pedestals/Pedestal_${run}.txt\"\)

done

cd results_proto
hadd MIPs_15_layers_${run}_0.root MIPs_15_layers_${run}_0*.root
hadd MIPs_15_layers_${run}_1.root MIPs_15_layers_${run}_1*.root
#hadd MIPs_15_layers_${run}_2.root MIPs_15_layers_${run}_2*.root

hadd MIPs_15_layers_${run}.root MIPs_15_layers_${run}_0.root MIPs_15_layers_${run}_1.root
#MIPs_15_layers_${run}_2.root 
source analysis.sh ${run}

cd -
source analysis_run_050xxx_monitoring.sh
