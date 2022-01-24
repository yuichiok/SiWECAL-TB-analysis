

run=$1
gain=$2
sca=$3
output="/lustre/ific.uv.es/prj/ific/flc/SiWECAL/TB2021/rootfiles/3GeVMIPscan/"

initial_folder=$PWD

cd $initial_folder
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-01/init_ilcsoft.sh
root -l -q /lhome/ific/a/airqui/SiWECAL/SiWECAL-TB-analysis-Analysis2021/SLBperformance/Noise.cc\(\"${output}/${run}_merged\",\"${run}\",${gain},${sca}\) 
#root -l -q Noise.cc\(\"${output}/${run}_merged\",\"${run}\",0\)

