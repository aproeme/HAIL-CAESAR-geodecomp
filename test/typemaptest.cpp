#include <mpi.h>

#include <typemaps.h>

#include "catchmentmodel/LSDCatchmentModel.hpp"

int main(int argc, char *argv[])
{
  std::cout << "running typemaps test\n";
  
  MPI::Init(argc, argv);
  Typemaps::initializeMaps();
  
  Cell *sendcell = new Cell(Cell::CellType::INTERNAL, 0.0, 0.0, 0.0, 0.0, false);
  Cell *recvcell = new Cell(Cell::CellType::INTERNAL, 0.0, 0.0, 0.0, 0.0, false);
  
  
  int tag = 13513;
  MPI::Request requests[2];
    
  std::cout << "sending...\n";
  requests[0] = MPI::COMM_WORLD.Isend(&sendcell, 1, MPI_CELL, 0, tag);
  std::cout << "receiving...\n";
  requests[1] = MPI::COMM_WORLD.Irecv(&recvcell, 1, MPI_CELL, 0, tag);
  MPI::Request::Waitall(2, requests);
  std::cout << "done.\n";
  
  MPI::Finalize();
  return 0;
}
