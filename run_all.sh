#/usr/bin/bash



for n in {1000000,10000000}; do
	for p in {2,4,6,8,10,12,16,20,24,28,32,36,40,44,48,52,56,60,65,70,75,80,85,90,95,100,105,110,115,120,128}; do
		for dist in {10,20,30,40,50,60,70,80,100,150,200,300,400,500,600,700,800,1000,1200,1500,2000}; do
			echo "mpirun -np " $p " build/Code ruling_set2 " $n " " $dist
			mpirun -np $p build/Code ruling_set2 $n $dist
		done
	done
done 