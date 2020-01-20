# HAIL-CAESAR

This is the HPC development version of HAIL-CAESAR arrived at by porting the code at https://github.com/dvalters/hail-caesar to use LibGeoDecomp for MPI-based multi-node parallelism.

# Requirements

To build the HPC version of HAIL-CAESAR you will need the following installed:

- g++ (more recent than version 6)
- an MPI library (e.g. MPICH, OpenMPI, ...)
- Boost
- LibGeoDecomp (currently the hail-caesar branch, available from https://github.com/aproeme/libgeodecomp)


# Building

- Run ```cp ./include/make/make.inc.template ./make.inc```
- Modify ```make.inc``` to specify the locations of Boost, LibGeoDecomp, and the MPI library, and the MPI compiler wrapper command to compile C++ code (e.g. mpic++)
- Run ```make```

You will be able to run the resulting executable ./bin/HAIL-CAESAR.mpi as follows (e.g. on ARCHER):

aprun -n 24 HAIL-CAESAR.mpi path_to_params_file params_filename

The params file is used as in the original version of HAIL-CAESAR as described as http://hail-caesar.readthedocs.io/en/latest/

For simple synthetic test cases, see /test/synthetic