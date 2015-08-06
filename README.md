# expFit
Fits single term, double term and triple term exponentials to data.

Example usage:

mpirun -n 4 expFit.x data3.in 2

Above runs through mpi, 4 processes, executable is expFit.x, data table is in data3.in and it's trying to fit a two term exponential.

*6th August 2015*

The program no longer needs mpi.  Just compile with make and run "mpirun -n 1 expFit.x" or "./expFit.x"
