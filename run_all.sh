#/usr/bin/bash

p_min=2
p_max=8
p_step=1

n_min=1000
n_max=10000
n_fac=2

ruler_dist_min=3
ruler_dist_max=6
ruler_dist_step=1

for dist in {10,20,30,50,75,100,200,300,400,500,600,700,800,900,1000}; do

	mpirun -np $1 build/Code ruling_set $2 $dist
done

