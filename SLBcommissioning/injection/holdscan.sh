#!/bin/bash
#Expect files in this format:
# /path/run_XXXXX_SLB_2.root (or SLB_3)
# run example:
#> source analysis.sh XXXXX convert
# where:
#          XXXXX is the run number
#          convert=0 if no conversiom =1 if yes

#source analysis.sh 21010 1

#Run_ILC_test_cosmic_02192020_18h30min_Ascii"
data_folder="/mnt/win2/Run_data/"
run="20200226_dac1.15V_chn0to3"

cd ../../converter_SLB


#for irun in {120..160..20}
#do
irun=80
    folder=${data_folder}${run}"_hold"${irun}"_Ascii/"
    output="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis_TB2020/converter_SLB/convertedfiles/"${run}"_Ascii/"
    mkdir ${output}
    root -l -q ConvertDirectorySL.cc\(\"${folder}\",false,\"${output}\"\)
    hadd ${output}/holdscan_hold${irun}.root ${output}/conv*.root
    rm ${output}/conv*.root
#done

cd -
#root -l Holdscan4Channels.cc+

#hadd ${output}/${run}.root ${output}/*.root 

#root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",1\) &
