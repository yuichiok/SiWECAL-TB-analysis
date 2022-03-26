for sca in {0..14}
do
    run="03102022_pedestal_13slabs"
    cd ../
    root -l -q Noise.cc\(\"../converter_SLB/convertedfiles/${run}.root\",\"${run}\",0,${sca}\) &
    root -l -q Noise.cc\(\"../converter_SLB/convertedfiles/${run}.root\",\"${run}\",1,${sca}\)
    cd -
done


	
