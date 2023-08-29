#/usr/bin/bash



for n in {1000000,10000000}; do
	for p in {2,4,6,8,10,12,16,20,24,28,32,36,40,44,48,52,56,60,65,70,75,80,85,90,95,100,110,128}; do
		for dist in {8,10,15,20,30,40,50,60,70,80,90,100,150}; do
			echo "mpirun -np " $p " build/Code ruling_set " $n " " $dist
			mpirun -np $p build/Code ruling_set $n $dist
		done
	done
done 