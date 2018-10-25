energy=$1

if [ $energy == 40 ]
then
## 40 GeV Pions +
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744111_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744112_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744113_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744119_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744124_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744125_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_40GeV/ 744339_    
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_40GeV/ 744340_
fi

if [ $energy == 50 ]
then
## 50 GeV Pions +
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744331_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744328_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744327_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_50GeV/ 744326_
fi

if [ $energy == 60 ]
then
## 60 GeV Pions +                                                                                                                                       
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run74744132_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744134_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744141_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744143_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_60GeV/ 744323_
fi


## 70 GeV Pions +
if [ $energy == 70 ]
then
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744145_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744152_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/default/ run744164_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_70GeV/ 744318_
fi

if [ $energy == 80 ]
then
## 80 GeV Pions +                                                                                                                                       
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_80GeV/ 744317_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_80GeV/ 744316_
    source build_script_muon.sh true true /eos/project/s/siw-ecal/TB2018-09/CERN2018_pc2/Pions_pos_80GeV/ 744315_

fi

