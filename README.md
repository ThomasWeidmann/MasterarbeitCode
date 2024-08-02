# Master's thesis

This is the code of my [**master's thesis**](https://ae.iti.kit.edu/theses.php) (my thesis is not uploaded yet as of 02.08.2024) where i coded a distributed algorithm for list ranking and tree rooting. The algorithms were evaluated on different instances up to 16384 cores on the SuperMUC-NG supercomputer at the Leibniz Supercomputing Center.

The problem of list ranking determines for every vertex in a list the distance to the end of the list. The input can be generalized to a forest. This problem is denoted by tree rooting. Here each vertex aims to determine its distance to the tree’s root and identify the root itself.

## Usage
The algorithms use the Message Passing Interface (MPI) to communicate. The different algorithms and test instances are further explained in my [**master's thesis**](https://ae.iti.kit.edu/theses.php).