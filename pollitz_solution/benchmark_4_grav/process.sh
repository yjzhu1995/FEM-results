#! /bin/bash
awk '{if (NR>6 && NR<=1006) print $1,$2}' strainx-0.5.inTHRUST > coor
for i in 0.5 5 50 500
do
python ./read.py strainA-${i}.outTHRUST 1000 temp
paste coor temp > pollitz${i}_grav.dis
done
