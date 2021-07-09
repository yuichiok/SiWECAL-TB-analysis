#!/bin/bash

run="run_050001"
run_file="run_050001_25062021_19h25min_Ascii.dat"
data_folder="/lustre/ific.uv.es/prj/ific/flc/SiWECAL/TB2021/run_050XXX/run_050001_25062021_19h25min_Ascii"
output="/lhome/ific/a/airqui/SiWECAL/SiWECAL-TB-analysis/converter_SLB/convertedfiles/"${run}"/"
mkdir $output

initial_folder=$PWD

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
    #    cat > ${initial_folder}/${run}_${j}_conv.sh <<EOF
    #source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-01/init_ilcsoft.sh
    cd ${output}/
    tar xzvf ${data_folder}/${run_file}_$j.tar.gz
    cd $initial_folder/../converter_SLB
    root -l -q ConvertDataSL.cc\(\"${output}/${run_file}_$j\",false\)
    rm ${output}/${run_file}_$j
    #EOF

    #    cat > ${initial_folder}/${run}_${j}_conv.sub <<EOF
    ## Unix submit description file
    ## ${run}_${j}.sub -- simple scriptjob
    #executable              = ${run}_${j}_conv.sh
    #log                     = log/${run}_${j}_conv.log
    #output                  = log/outfile_${run}_${j}_conv.txt
    #error                   = log/errros_${run}_${j}_conv.txt
    #should_transfer_files   = Yes
    #when_to_transfer_output = ON_EXIT
    #queue 1
    #EOF

    #condor_submit ${run}_${j}_conv.sub
done


