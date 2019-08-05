#root -l pedanalysis0.cc\(21002,2\) pedanalysis0.cc\(21003,2\) pedanalysis0.cc\(21004,2\) pedanalysis0.cc\(21005,2\) pedanalysis0.cc\(21007,2\) pedanalysis0.cc\(21010,2\) pedanalysis0.cc\(21011,2\) pedanalysis0.cc\(21012,2\) pedanalysis0.cc\(21013,2\) 

for i in {21001..21021}
do
    root -l -q pedanalysis0.cc\($i,2\) &
    root -l -q pedanalysis0.cc\($i,3\) 
done

for i in {21027..21048}
do
    root -l -q pedanalysis0.cc\($i,2\) &
    root -l -q pedanalysis0.cc\($i,3\)
done
