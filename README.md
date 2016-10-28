# MPI_sendrecv
This is a simple FORTRAN program for sending and receiving an array in MPI. 

The program would be useful for scientific computations.


--> Command for compiling (intel compilers): <mpiifort array_sendrecv_test.f90> (ignore <>)
--> If the above command did not work, try with <mpiifort array_sendrecv_test.f90 -heap-arrays>

--> To run the program <mpiexec -np <num_of_processors> ./a.out>
