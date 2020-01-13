# HAIL-CAESAR

This is the HPC development version of HAIL-CAESAR arrived at by porting the code at https://github.com/dvalters/hail-caesar to use LibGeoDecomp for MPI-based multi-node parallelism.

# Requirements

To build the HPC version of HAIL-CAESAR you will need the following installed:

- g++ (more recent than version 6)
- an MPI library (e.g. MPICH, OpenMPI, ...)
- Boost
- LibGeoDecomp (currently the hail-caesar branch, available from https://github.com/aproeme/libgeodecomp)
- ruby
- doxygen


# Building

To build, modify the Makefile to provide your MPI compiler commands (e.g. mpic++), the paths to your installations of the MPI library, Boost, and LibGeoDecomp
Then, run:

make typemaps

make

You will be able to run the resulting executable ./bin/HAIL-CAESAR.mpi as follows (e.g. on ARCHER):

aprun -n 24 HAIL-CAESAR.mpi path_to_params_file params_filename

The params file is used as in the original version of HAIL-CAESAR as described as http://hail-caesar.readthedocs.io/en/latest/

For simple synthetic test cases, see /test/synthetic