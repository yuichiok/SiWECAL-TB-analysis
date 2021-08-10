#!/bin/bash

#Run_ILC_test_cosmic_02192020_18h30min_Ascii"
data_folder="/mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/commissioning/injection_07282021/"
run="07282021_dac1.2V_small"

for chn in "_chn0to7" "_chn8to15" "_chn16to23" "_chn24to31" "_chn32t39" "_chn40to47" "_chn56to63"
do
    cd ../../converter_SLB
    for irun in {20..160..20}
    do
	#irun=80
	folder=${data_folder}${run}${chn}"_hold"${irun}"_Ascii/"
	output="/mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis-dev/converter_SLB/convertedfiles/"${run}${chn}"_Ascii/"
	mkdir ${output}
	root -l -q ConvertDirectorySL.cc\(\"${folder}\",false,\"${output}\"\)
	hadd ${output}/holdscan_hold${irun}.root ${output}/conv*.root
	rm ${output}/conv*.root
    done
    cd -
done

mkdir "/mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis-dev/converter_SLB/convertedfiles/"${run}

for irun in {20..160..20}
do
    hadd /mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis-dev/converter_SLB/convertedfiles/${run}/holdscan_hold${irun}.root /mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis-dev/converter_SLB/convertedfiles/${run}${chn}_Ascii/holdscan_hold${irun}.root
done

root -l HoldscanGraphs.cc
