#!/bin/bash
gfortran *.f90
mkdir ../viscoGreen/GREEN
fl=('0' '10' '20' '29' '35' '40' '50' '60' '70' '80' '90' '100')
for item in ${fl[*]}
do
	cp '../viscoGreen/GREEN_'$item/* ../viscoGreen/GREEN/
	./a.out
	mv displacements.dat ./disp_$item.dat
done
rm -rf ../viscoGreen/GREEN
