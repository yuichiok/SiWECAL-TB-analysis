run=$1
#root -l -q analysis.cc\(\"${run}\",\"high\",true,false,1\) &
root -l -q analysis.cc\(\"${run}\",\"high\",false,true,1\)  &
#root -l -q analysis.cc\(\"${run}\",\"low\",true,false,1\) &
root -l -q analysis.cc\(\"${run}\",\"low\",false,true,1\)

