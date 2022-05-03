run=$1
output="/eos/project/s/siw-ecal/TB2021-11/beamData/rootfiles/3GeVMIPscan/"

for sca in {0..14}
do
    cd ../
    root -l -q Noise.cc\(\"${output}/${run}_merged.root\",\"${run}\",0,${sca}\)
    root -l -q Noise.cc\(\"${output}/${run}_merged.root\",\"${run}\",1,${sca}\)
    cd -
done


	
