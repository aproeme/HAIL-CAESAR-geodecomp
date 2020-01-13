#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>

#include <libgeodecomp/communication/mpilayer.h>

#include "catchmentmodel/LSDCatchmentModel.hpp"
#include "catchmentmodel/LSDUtils.hpp"

#ifndef GIT_REVISION
#define GIT_REVISION "N/A"
#endif
#define CHM_VERS 1.0

using namespace LSDUtils;

#ifdef COMPILE_FOR_PARALLEL
#include <mpi.h>
#include <typemaps.h>
#include "typemaps.h"
#include <libgeodecomp/communication/typemaps.h>

#endif


int main(int argc, char *argv[])
{
#ifdef COMPILE_FOR_PARALLEL
  MPI_Init(&argc, &argv);
  std::string pname(argv[1]);
  std::string pfname(argv[2]);
  // The path name and the parameter file name, respectively.

  LibGeoDecomp::Typemaps::initializeMaps(); // initialize LibGeoDecomp native typemaps (this commits MPI types)
  Typemaps::initializeMaps(); // initialize custom typemaps for HAIL-CAESAR    
  LibGeoDecomp::MPILayer().barrier();
  
  if(LibGeoDecomp::MPILayer().rank() == 0)
    {
      std::cout << "##################################" << std::endl;
      std::cout << "#  CATCHMENT HYDROGEOMORPHOLOGY  #" << std::endl;
      std::cout << "#        MODEL version ?.?       #" << std::endl;
      std::cout << "#          (HAIL-CAESAR)         #" << std::endl;
      std::cout << "##################################" << std::endl;
      std::cout << " Version: "<< CHM_VERS << std::endl;
      std::cout << " at git commit number: " GIT_REVISION << std::endl;
      std::cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
      
      std::cout << "The pathname is: " << pname
		<< " and the parameter file is: " << pfname << std::endl;
    }
  
  LibGeoDecomp::MPILayer().barrier();

  runSimulation(pname, pfname);
  MPI_Finalize();
#endif
  
  return 0; 
}


