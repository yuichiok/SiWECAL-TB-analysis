#!/bin/bash

#run
run="run_050001"
run_file="run_050001_25062021_19h25min_Ascii.dat"
#ascii folder
data_folder="/lustre/ific.uv.es/prj/ific/flc/SiWECAL/TB2021/run_050XXX/run_050001_25062021_19h25min_Ascii"
#rootfiles
output="/mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis/converter_SLB/convertedfiles/"${run}"/"
#mkdir $output

initial_folder=$PWD

#source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-01/init_ilcsoft.sh      

for i in {1..2339}
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
    #conversion
    #source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-01/init_ilcsoft.sh
    #analysis
    cd $initial_folder
    root -l -q Proto.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\",\"mip\",\"../pedestals/Pedestal_${run}.txt\"\)

done

cd results_proto
hadd MIPs_15_layers_${run}_0.root MIPs_15_layers_${run}_0*.root
hadd MIPs_15_layers_${run}_1.root MIPs_15_layers_${run}_1*.root
hadd MIPs_15_layers_${run}_2.root MIPs_15_layers_${run}_2*.root

hadd MIPs_15_layers_${run}.root MIPs_15_layers_${run}_0.root MIPs_15_layers_${run}_1.root MIPs_15_layers_${run}_2.root 
source analysis.sh ${run}
