run=$1
output="/eos/project/s/siw-ecal/TB2021-11/beamData/rootfiles/3GeVMIPscan/"
initial_folder=$PWD/../

cd $initial_folder
for sca in {0..14}
do
    root -l -q Noise.cc\(\"${output}/${run}_merged.root\",\"${run}\",0,${sca}\)
    root -l -q Noise.cc\(\"${output}/${run}_merged.root\",\"${run}\",1,${sca}\)
done
cd -

	
