#!/bin/bash
gfortran *.f90 *.f -fopenmp -llapack -lblas -o omp.out
fl=(0 10 20 29 35 40 50 60 70 80 90 100)
for item in ${fl[*]}
do 
	mkdir GREEN
	./omp.out << EOF
$item
EOF
	mv ./GREEN ./GREEN_$item
done
