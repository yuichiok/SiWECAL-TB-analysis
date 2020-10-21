#!/bin/bash
#Expect files in this format:
# /path/run_XXXXX_SLB_2.root (or SLB_3)
# run example:
#> source analysis.sh XXXXX convert
# where:
#          XXXXX is the run number
#          convert=0 if no conversiom =1 if yes

#source analysis.sh 21010 1

for run in "Run_ILC_10162020_17h38min_cosmic_Ascii" "Run_ILC_10192020_17h38min_cosmic_Ascii" "Run_ILC_10202020_11h44min_cosmic_Ascii"
do
#run="Run_ILC_10162020_17h38min_cosmic_Ascii"
    
    data_folder="/mnt/win2/Run_Data/"${run}"/"
    
    mkdir "../converter_SLB/convertedfiles/"${run}
    output="../converter_SLB/convertedfiles/"${run}"/"
    
    cd ../converter_SLB
    root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
    hadd ${output}/${run}.root ${output}/*.root
    cd -
    
    root -l SlowControlMonitoring.cc\(\"${output}/${run}\",\"${run}\"\) 
    root -l DummyDisplay.cc\(\"${output}/${run}\",\"${run}\"\)
    
    #root -l Monitoring.cc\(\"${output}/${run}\",\"${run}\",1\) 
    
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",0\) &
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",1\) &
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",2\) &
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",3\) 
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",4\) &
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",5\) 
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",6\) &
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",7\) &
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",8\) 
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",9\) &
    root -l SingleSlabAnalysis.cc\(\"${output}/${run}\",\"${run}\",10\) 
    
    sleep 10
    
    cd performance
    source analysis.sh $run 1 &
    source analysis.sh $run 2 &
    source analysis.sh $run 3
    source analysis.sh $run 4 &
    source analysis.sh $run 5 &
    source analysis.sh $run 6 &
    source analysis.sh $run 7
    source analysis.sh $run 10 &
    source analysis.sh $run 11
    source analysis.sh $run 13 &
    source analysis.sh $run 14

    cd -
done
