energy=$1

if [ $energy == 40 ]
then
## 40 GeV Pions +
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744111_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744112_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744113_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744119_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744124_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744125_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_40GeV/ 744339_    
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_40GeV/ 744340_

    cd ../built_files/
    hadd ECAL_40GeV_pions_plus_build.root run744111*build.root run744112*build.root run744113*build.root run744119*build.root run744124*build.root run744124*build.root run744125*build.root 744339_*build.root 744340_*build.root
    mv *root ../Common/ECAL/offset_twiki/PiPlus_40GeV/.
    cd -
fi

if [ $energy == 50 ]
then
## 50 GeV Pions +
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744331_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744328_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744327_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744326_
    
    cd ../built_files/
    hadd  ECAL_50GeV_pions_plus_build.root  744326_*build.root 744327_*build.root 744328_*build.root 744331_*build.root
    mv *root ../Common/ECAL/offset_twiki/PiPlus_50GeV/.
    cd -
fi

if [ $energy == 60 ]
then
## 60 GeV Pions +                                                                                                                                       
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run74744132_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744134_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744141_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744143_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_60GeV/ 744323_

    cd ../built_files/
    hadd  ECAL_60GeV_pions_plus_build.root  run74744132*build.root run744134*build.root run744141*build.root run744143*build.root 744323_*build.root
    mv *root ../Common/ECAL/offset_twiki/PiPlus_60GeV/.
    cd -
fi


## 70 GeV Pions +
if [ $energy == 70 ]
then
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744145_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744152_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744164_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_70GeV/ 744318_

    cd ../built_files/
    hadd  ECAL_70GeV_pions_plus_build.root  run744145*build.root run744152*build.root run744164*build.root 744318_*build.root
    mv *root ../Common/ECAL/offset_twiki/PiPlus_70GeV/.
    cd -
fi

if [ $energy == 80 ]
then
## 80 GeV Pions +                                                                                                                                       
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_80GeV/ 744317_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_80GeV/ 744316_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_80GeV/ 744315_

    cd ../built_files/
    hadd  ECAL_80GeV_pions_plus_build.root  744317_*build.root 744316_*build.root 744315_*build.root 
    mv *root ../Common/ECAL/offset_twiki/PiPlus_80GeV/.
    cd -
fi


if [ $energy == "Muon" ]
then
## Muons 200GeV
 
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744211_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744214_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744218_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744221_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744226_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744230_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744233_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744237_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744242_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744246_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744249_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744254_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744258_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744263_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744269_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744273_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744278_
    source build_script.sh false true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/muon/ 744283_

    cd ../built_files/
    hadd  ECAL_200GeV_muons_build.root  7442*build.root 
    mv *root ../Common/ECAL/offset_twiki/Muon_200GeV/.
    cd -
fi

