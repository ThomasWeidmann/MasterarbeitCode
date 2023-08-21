#/usr/bin/bash



for dist in {50,100,200,300,400,500,600,700,800,900,1000,1200,1500,2000,3000}; do

	mpirun -np $1 build/Code ruling_set $2 $dist
done

