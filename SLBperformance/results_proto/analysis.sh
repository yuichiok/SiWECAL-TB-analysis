run=$1
for slab in {0..14}
do
    root -l -q analysis.cc\(\"$run\",$slab\) 
done

