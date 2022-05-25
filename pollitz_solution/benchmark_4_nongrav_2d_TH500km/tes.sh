#! /bin/bash

for i in $(ls strainx*.inTHRUST)
do
    cp $i ${i}-dis
    #awk '{if (NR == 1007) print 1; else print $1,$2,$3,$4,$5}' $i > ${i}-vel
    sed  "1007s/0/1/g" $i > $i-vel
done
