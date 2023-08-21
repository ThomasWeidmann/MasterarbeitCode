#/usr/bin/bash
#i=10
#while [[ $i -le 1000 ]] ; do
for i in {10,20,50,100,500,600,1000} ; do
   mpirun -np 4 code 10000000 $i
  (( i += i ))
done