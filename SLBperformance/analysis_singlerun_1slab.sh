#!/bin/bash
#Expect files in this format:
# /path/run_XXXXX_SLB_2.root (or SLB_3)
# run example:
#> source analysis.sh XXXXX convert
# where:
#          XXXXX is the run number
#          convert=0 if no conversiom =1 if yes

#source analysis.sh 21010 1

run="Run_ILC_test_cosmic_02192020_18h50min_Ascii"
#Run_ILC_test_cosmic_02192020_18h30min_Ascii"
data_folder="/mnt/win2/Run_Data/"${run}"/"

mkdir "/home/calice/TB2020/commissioning/SiWECAL-TB-analysis_TB2020/converter_SLB/convertedfiles/"${run}
output="/home/calice/TB2020/commissioning/SiWECAL-TB-analysis_TB2020/converter_SLB/convertedfiles/"${run}"/"

cd ../converter_SLB
root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\) 
hadd ${output}/${run}.root ${output}/*.root 

cd -
root -l Monitoring.cc\(\"${output}/${run}\",\"${run}\",1\) 

root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",1\) &
root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",2\) 
root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",4\) &
root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",5\)


sleep 5

cd performance
source analysis.sh $run 1 &
source analysis.sh $run 2
source analysis.sh $run 4 &
source analysis.sh $run 5

