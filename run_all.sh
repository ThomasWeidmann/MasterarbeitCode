#/usr/bin/bash


for n in {1000000,10000000,100000000,1000000000}; do
	for p in {2,3,4,6,8,12,16,24,32}; do
		mpirun -np $p build/Code sequential $n
		mpirun -np $p build/Code pointer_doubling $n
		for dist in {50,100,200,300,400,500,600,700,800,900,1000,1200,1500,2000,3000}; do
			echo "mpirun -np 1 build/Code ruling_set " $n " " $dist
			mpirun -np $p build/Code ruling_set $n $dist
		done
	done
done 