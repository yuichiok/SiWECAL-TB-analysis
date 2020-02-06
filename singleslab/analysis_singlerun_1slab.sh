#!/bin/bash
#Expect files in this format:
# /path/run_XXXXX_SLB_2.root (or SLB_3)
# run example:
#> source analysis.sh XXXXX convert
# where:
#          XXXXX is the run number
#          convert=0 if no conversiom =1 if yes

#source analysis.sh 21010 1

run="Run_ILC_slab21_02042020_19h29min_Ascii"
#Run_ILC_slab14_HV100_12192019_14h58min_Ascii"
data_folder="/home/calice/TB2020/commissioning/data_slboard/"${run}"/"
runname=$run

cd ../converter_SLB
root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",1,1\) 
hadd ${data_folder}/../${runname}_SLB_1.root ${data_folder}/${run}*SLB_1.root 


#data_folder1="/media/irles/783C4D5D3C4D1790/data_slboard/"${run1}"/"
#data_folder1="/home/irles/WorkAreaECAL/2019/data_slboard/"${run1}"/"

##root -l -q ConvertDirectorySL.cc\(\"${data_folder1}\",0,3\) &
##root -l -q ConvertDirectorySL.cc\(\"${data_folder1}\",1,3\) &
##root -l -q ConvertDirectorySL.cc\(\"${data_folder1}\",2,3\)

#cd -
##sleep 180
cd -


root -l CommissioningAnalysis.cc\(\"${data_folder}/../${runname}\",\"${runname}\",\"_SLB_1\"\) 

cd performance
source analysis.sh $runname 1

#source analysis_singlerun_1slab_bis.sh
