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
run="calib_02272020_pa1.2fb_trig230_chn0.8.16.24_calib_small"

cd ../../converter_SLB


for irun in 0.0 0.1 0.2 0.3 0.4 0.5 0.65 0.7 0.8 0.9 1.0 1.15 1.3 1.5 1.7
do
    #irun=120
    folder=${data_folder}${run}"_"${irun}"V_Ascii/"
    output="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis_TB2020/converter_SLB/convertedfiles/"${run}"_Ascii/"
    mkdir ${output}
    root -l -q ConvertDirectorySL.cc\(\"${folder}\",false,\"${output}\"\)
    hadd ${output}/injection_small_${irun}V.root ${output}/conv*.root
    rm ${output}/conv*.root
done

cd -
#root -l Linearity4Channels.cc+
