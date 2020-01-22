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

#include <mpi.h>
#include <typemaps.h>
#include "typemaps.h"
#include <libgeodecomp/communication/typemaps.h>



int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  
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
      
      if (argc < 2)
	{
	  std::cout << "\n###################################################" << std::endl;
	  std::cout << "No parameter file supplied" << std::endl;
	  std::cout << "You must supply a path and parameter file!" << std::endl;
	  std::cout << "see https://dvalters.github.io/HAIL-CAESAR/" << std::endl;
	  std::cout << "for assistance." << std::endl;
	  std::cout << "###################################################" << std::endl;
	  
	  exit(0);
	}
      
      if (argc > 2)
	{
	  std::cout << "Too many input arguments supplied (should be 3...)" << std::endl;
	  exit(0);
	}

    }

  LibGeoDecomp::MPILayer().barrier();
  
  std::string pfname(argv[1]);
  
  if(LibGeoDecomp::MPILayer().rank() == 0)
    {
      std::cout << "Parameter file is: " << pfname << std::endl;  
    }
    
  runSimulation(pfname);
  MPI_Finalize();
  
  return 0; 
}


