#include <cmath>
#include <sys/stat.h> 

#include <boost/assign/std/vector.hpp>

#include "typemaps.h"
#include <libgeodecomp/io/ppmwriter.h>
#include <libgeodecomp/misc/apitraits.h>
#include <libgeodecomp/io/tracingwriter.h>
#include <libgeodecomp/io/simpleinitializer.h>
#include <libgeodecomp/parallelization/serialsimulator.h>
#include <libgeodecomp/io/asciiwriter.h>
#include <libgeodecomp/io/mpiiowriter.h>
#include <libgeodecomp/communication/typemaps.h>
#include <libgeodecomp/io/collectingwriter.h>
#include <libgeodecomp/io/mpiioinitializer.h>
#include <libgeodecomp/io/parallelmpiiowriter.h>
#include <libgeodecomp/parallelization/stripingsimulator.h>
#include <libgeodecomp/parallelization/hiparsimulator.h>
#include <libgeodecomp/loadbalancer/noopbalancer.h>
#include <libgeodecomp/io/bovwriter.h>
#include <libgeodecomp/geometry/partitions/recursivebisectionpartition.h>

#include "catchmentmodel/cell.hpp"
#include "catchmentmodel/LSDCatchmentModel.hpp"
#include "catchmentmodel/LSDUtils.hpp"

  

// Initialising default values for static Cell variables
double LSDCatchmentModel::DX = 1.0;
double LSDCatchmentModel::DY = LSDCatchmentModel::DX;
double LSDCatchmentModel::no_data_value = -9999;
double LSDCatchmentModel::water_depth_erosion_threshold = 1.0;
double LSDCatchmentModel::edgeslope = 0.001;
double LSDCatchmentModel::hflow_threshold = 0.00001;
double LSDCatchmentModel::mannings = 0.04;
double LSDCatchmentModel::froude_limit = 0.8;
double LSDCatchmentModel::time_factor = 1;
double LSDCatchmentModel::courant_number = 0.7;
double LSDCatchmentModel::maxdepth = 10;


using namespace LSDUtils;



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// CREATE FUNCTIONS
// These define what happens when an LSDCatchmentModel object is created
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDCatchmentModel::create()
{
  std::cout << "You are trying to create an LSDCatchmentModel object with no"
            << "supplied files or parameters." << std::endl
            << "Exiting..." << std::endl;
  exit(EXIT_FAILURE);
}

void LSDCatchmentModel::create(std::string pfname)
{
  LSDCatchmentModel::initialise_variables(pfname);
  if(LibGeoDecomp::MPILayer().rank() == 0)
    {
      std::cout << "The user-defined parameters have been"
		<< " ingested from the param file." << std::endl;
    }
}
















// Cell default constructor
Cell::Cell(Cell::CellType celltype_in = INTERNAL,	\
	   double elevation_in = 0.0,			\
	   double water_depth_in = 0.0,			\
	   double qx_in = 0.0,				\
	   double qy_in = 0.0) : celltype(celltype_in), elevation(elevation_in), water_depth(water_depth_in), qx(qx_in), qy(qy_in)
{}




template<typename Cell>
class LibGeoDecomp::APITraits::SelectMPIDataType<Cell>  
{
public:
  static inline MPI_Datatype value()
  {
    return MPI_CELL;  // MPI_CELL is initialized by custom-generated Typemaps::initializeMaps() for HAIL-CAESAR Cell class
  }
};



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Overall update routine - this is what LibGeoDecomp calls each time step
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename COORD_MAP>
void Cell::update(const COORD_MAP& neighborhood, unsigned nanoStep)
{
  set_global_timefactor();
  
  // Hydrological and flow routing processes
  // Add water to the catchment from rainfall input file
  catchment_waterinputs(neighborhood);
  // Distribute the water with the LISFLOOD Cellular Automaton algorithm
  flow_route_x(neighborhood);
  flow_route_y(neighborhood);
  // Calculate the new water depths in the catchment
  depth_update(neighborhood);
  // Water outputs from edges/catchment outlet 
  water_flux_out(neighborhood);
}












void Cell::set_global_timefactor()
{
  if (LSDCatchmentModel::maxdepth <= 0.1)
    {
      LSDCatchmentModel::maxdepth = 0.1;
    }
  if (LSDCatchmentModel::time_factor < (LSDCatchmentModel::courant_number * (LSDCatchmentModel::DX / std::sqrt(Cell::gravity * (LSDCatchmentModel::maxdepth)))))
    {
     LSDCatchmentModel::time_factor = (LSDCatchmentModel::courant_number * (LSDCatchmentModel::DX / std::sqrt(Cell::gravity * (LSDCatchmentModel::maxdepth))));
    }
  /*    if (input_output_difference > in_out_difference_allowed &&LSDCatchmentModel::time_factor > (LSDCatchmentModel::courant_number * (LSDCatchmentModel::DX / std::sqrt(Cell::gravity * (LSDCatchmentModel::maxdepth)))))
	{
  	time_factor =LSDCatchmentModel::courant_number * (LSDCatchmentModel::DX / std::sqrt(Cell::gravity * (LSDCatchmentModel::maxdepth)));
  	}*/

}



template<typename COORD_MAP>
void Cell::catchment_waterinputs(const COORD_MAP& neighborhood) // refactor - incomplete (include runoffGrid for complex, i.e. spatially variable rainfall)
{
  
  //waterinput = 0;
  double local_time_factor = set_local_timefactor();
  /*
  if (spatially_complex_rainfall == true)
    {
      catchment_water_input_and_hydrology(local_time_factor, runoff);
    }
    else */
  //  {
  catchment_water_input_and_hydrology(neighborhood, local_time_factor);
  //    }
}

template<typename COORD_MAP>
void Cell::catchment_water_input_and_hydrology(const COORD_MAP& neighborhood, double local_time_factor)
{
  /*
  for (unsigned i = 1; i<=rfnum; i++)
    {
      waterinput += j_mean[i] * nActualGridCells[i] * LSDCatchmentModel::DX * LSDCatchmentModel::DX;
    }
  
    for (int z=1; z <= totalinputpoints; z++)
    {
      int i = catchment_input_x_coord[z];
      int j = catchment_input_y_coord[z];
      double water_add_amt = (j_mean[rfarea[i][j]] * nActualGridCells[rfarea[i][j]]) /
        (catchment_input_counter[rfarea[i][j]]) * local_time_factor;    //
      if (water_add_amt > ERODEFACTOR)
	{
	  water_add_amt = ERODEFACTOR;
	}
      
  */
	
  water_depth = water_depth;
  
      /*
	}
      */
  
  // if the input type flag is 1 then the discharge is input from the hydrograph
  /*  if (cycle >= time_1)
    {
      do
	{
	  time_1++;
	  topmodel_runoff(time_1);  // calc_J is based on the rainfall rate supplied to the cell

	  if (time_factor > max_time_step) // && new_j_mean[1] > (0.2 / (jmax * imax * LSDCatchmentModel::DX * LSDCatchmentModel::DX)))
	    {
	      // Find the current maximum runoff amount
	      double j_mean_max_temp = 0;
	      for (unsigned n = 1; n <= rfnum; n++)
		{
		  if (new_j_mean[n] > j_mean_max_temp)
		    {
		      j_mean_max_temp = new_j_mean[n];
		    }
		}

	      // check after the variable rainfall area has been added
	      // stops code going too fast when there is actual flow in the channels greater than 0.2cu
	      if (j_mean_max_temp > (0.2 / (imax * jmax * LSDCatchmentModel::DX * LSDCatchmentModel::DX)))
		{
		  cycle = time_1 + (max_time_step / 60);
		 LSDCatchmentModel::time_factor = max_time_step;
		}
	    }
	} while (time_1 < cycle);
    }
  

  calchydrograph(time_1 - cycle);

  double jmeanmax =0;
  for (unsigned n=1; n <= rfnum; n++)
    {
      if (j_mean[n] > jmeanmax)
	{
	  jmeanmax = j_mean[n];
	}
    }
  */


  /*
  // DV - This is for reading the discharge direct from an input file
  if (jmeaninputfile_opt == true)
    {
      j_mean[1] = ((hourly_rain_data[(static_cast<int>(cycle / rain_data_time_step))][0] //check in original
		    / std::pow(LSDCatchmentModel::DX, 2)) / nActualGridCells[1]);
    }

  if (jmeanmax >= baseflow) // > baseflow
    {
      baseflow = baseflow * 3;    // Magic number 3!? - DAV
      zero_and_calc_drainage_area();         // Could this come from one of the LSDobject files? - DAV
      get_catchment_input_points();
    }

  if (baseflow > (jmeanmax * 3) && baseflow > 0.000001) // DV reverted to match CL 1.8f (formerly 10^-7)
    // Could make the baseflow threshold a parameter in param file?
    {
      baseflow = jmeanmax * 1.25;   // Where do these magic numbers come from? DAV
      zero_and_calc_drainage_area();
      get_catchment_input_points();
    }
  */

}



void LSDCatchmentModel::zero_values()
{
  for(unsigned i=0; i <= imax+1; i++)
    {
      for(unsigned j=0; j <= jmax+1; j++)
	{
	  Vel[i][j] = 0;
	  area[i][j] = 0;
	  elev[i][j] = -9999;
	  bedrock[i][j] = -9999;
	  init_elevs[i][j] = elev[i][j];
	  water_depth[i][j] = 0;
	  index[i][j] = -9999;
	  inputpointsarray[i][j] = false;

	  qx[i][j] = 0;
	  qy[i][j] = 0;

	  qxs[i][j] = 0;
	  qys[i][j] = 0;

	  for(int n=0; n<=8; n++)
	    {
	      vel_dir[i][j][n]=0;
	    }

	  Vsusptot[i][j] = 0;
	  rfarea[i][j] = 1;
	}
    }

  // These arrays are different dimensions
  // So they need to be zeroed differently
  for (unsigned i=0; i<=imax; i++)
    {
      for(unsigned j=0; j<=jmax; j++)
	{
	  if (vegetation_on)
	    {
	      veg[i][j][0] = 0;// elevation
	      veg[i][j][1] = 0; // densitj
	      veg[i][j][2] = 0; // jw density
	      veg[i][j][3] = 0; // height
	    }
	  edge[i][j] = 0;
	  edge2[i][j] = 0;
	}
    }

  for(unsigned i=1; i<((jmax*imax)/LIMIT); i++)
    {
      if (!hydro_only)
	{
	  for(unsigned j=0; j<=G_MAX; j++)
	    {
	      grain[i][j] = 0;
	    }
	  for(int z=0; z <= 9; z++)
	    {
	      for(unsigned j=0; j <= G_MAX-2; j++)
		{
		  strata[i][z][j] =0;
		}
	    }
	}

      catchment_input_x_coord[i] = 0;
      catchment_input_y_coord[i] = 0;
    }
}



// Initialise the arrays (as done in initialise() )
// Not sure the point of having them declared in header file if you
// can't resize them...surely this is duplicating array creation?? DV
// TO DO DAV - address above comment

// Initialise the relevant arrays
void LSDCatchmentModel::initialise_arrays()
{
  std::cout << "Cartesian imax (no. of rows): " << imax << \
               " Cartesian jmax (no. of cols): " << jmax << std::endl;

  // Need to change this so it does not waste memory assigning arrays
  // when running in hydro mode etc.
  elev = TNT::Array2D<double> (imax+2,jmax+2, -9999);
  water_depth = TNT::Array2D<double> (imax+2,jmax+2, 0.0);

  // Cast to int and then double, what?
  //old_j_mean_store = new double[(int)((maxcycle*60)/input_time_step)+10];
  old_j_mean_store = std::vector<double>
    (static_cast<int>((maxcycle*60)/input_time_step)+10);

  qx = TNT::Array2D<double> (imax + 2, jmax + 2, 0.0);
  qy = TNT::Array2D<double> (imax + 2, jmax + 2, 0.0);

  qxs = TNT::Array2D<double> (imax + 2, jmax + 2, 0.0);
  qys = TNT::Array2D<double> (imax + 2, jmax + 2, 0.0);

  Vel = TNT::Array2D<double> (imax + 2, jmax + 2, 0.0);

  area = TNT::Array2D<double> (imax+2, jmax + 2, 0.0);
  index = TNT::Array2D<int> (imax +2, jmax + 2, 0);
  elev_diff = TNT::Array2D<double> (imax + 2, jmax + 2);

  bedrock = TNT::Array2D<double> (imax+2, jmax+2, -9999);
  tempcreep = TNT::Array2D<double> (imax+2,jmax+2);
  init_elevs = TNT::Array2D<double> (imax+2,jmax+2, -9999);

  vel_dir = TNT::Array3D<double> (imax+2, jmax+2, 9, 0.0);

  Vsusptot = TNT::Array2D<double> (imax+2,jmax+2, 0.0);

  // Will come back to this later - DAV
  //if (vegetation_on)
  //{
    veg = TNT::Array3D<double> (imax+1, jmax+1, 4, 0.0);
  //}


  //cross_scan = TNT::Array2D<int> (imax+2,jmax+2, 0);
  down_scan = TNT::Array2D<int> (jmax+2, imax+2, 0);

  // line to stop max time step being greater than rain time step
  if (rain_data_time_step < 1) rain_data_time_step = 1;
  if (max_time_step / 60 > rain_data_time_step)
  {
    max_time_step = static_cast<int>(rain_data_time_step) * 60;
  }
  // see StackOverflow for how to set size of nested vector
  // http://stackoverflow.com/questions/2665936/
  //     is-there-a-way-to-specify-the-dimensions-of-a-nested-stl-vector-c
  // DEBUG

  #ifdef DEBUG
  long int xdim = static_cast<int>(maxcycle * (60 / rain_data_time_step)) + 100;
  std::cout << "Size of rain vector: " << xdim << std::endl;

  int short_xdim = static_cast<int>(maxcycle * (60 / rain_data_time_step)) \
                    + 100;
  std::cout << "Int size of rain vector: " << short_xdim << std::endl;
  #endif

  hourly_rain_data = std::vector< std::vector<float> >
    ( (static_cast<int>(maxcycle * (60 / rain_data_time_step)) + 100),
                        vector<float>(rfnum+1) );

  hourly_m_value = std::vector<double>
    (static_cast<int>(maxcycle * (60 / rain_data_time_step)) + 100);

  inputpointsarray = TNT::Array2D <bool> (imax + 2, jmax + 2);

  edge = TNT::Array2D<double> (imax+1,jmax+1, 0.0);
  edge2 = TNT::Array2D<double> (imax+1,jmax+1, 0.0);

  Tau = TNT::Array2D<double> (imax+2,jmax+2, 0.0);

  catchment_input_x_coord = std::vector<int> (jmax * imax, 0);
  catchment_input_y_coord = std::vector<int> (jmax * imax, 0);

  area_depth = TNT::Array2D<double> (imax + 2, jmax + 2, 0.0);


  // Grain arrays - still to port
  

  // Distributed Hydrological Model Arrays
  // j and jo need a very small initial value to get them going
  j = std::vector<double> (rfnum + 1, 0.000000001);
  jo = std::vector<double> (rfnum + 1, 0.000000001);
  j_mean = std::vector<double> (rfnum + 1);
  // std::vectors will be default initalised to 0, unless specified
  old_j_mean = std::vector<double> (rfnum + 1);
  new_j_mean = std::vector<double> (rfnum + 1);
  rfarea = TNT::Array2D<int> (imax + 2, jmax + 2, 0);
  nActualGridCells = std::vector<int> (rfnum + 1);
  catchment_input_counter = std::vector<int> (rfnum + 1);

  // Only need these ones for erosion-enabled simulation runs
  if (!hydro_only)
  {
    sr = TNT::Array3D<double> (imax + 2, jmax + 2, 10, 0.0);
    sl = TNT::Array3D<double> (imax + 2, jmax + 2, 10, 0.0);
    su = TNT::Array3D<double> (imax + 2, jmax + 2, 10, 0.0);
    sd = TNT::Array3D<double> (imax + 2, jmax + 2, 10, 0.0);
    ss = TNT::Array2D<double> (imax + 2, jmax + 2, 0.0);

    strata = TNT::Array3D<double> ( ((imax+2)*(jmax+2))/LIMIT , 10, G_MAX+1, 0.0);
    grain = TNT::Array2D<double> ( ((2+imax)*(jmax+2))/LIMIT, G_MAX+1 , 0.0);
    temp_grain = std::vector<double> (G_MAX+1, 0.0);
  }
  // Initialise suspended fraction vector
  // Tells program whether sediment is suspended fraction or not
  isSuspended = std::vector<bool>(G_MAX+1, false);
  //fallVelocity = std::vector<double>(G_MAX+1, 0.0);
  //set_fall_velocities();
  // Not entirely sure this is necessary? - TODO DAV
  zero_values();

}








// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// WATER FLUXES OUT OF CATCHMENT
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Calculate the water coming out and zero any water depths at the edges
// This will actually set it to the minimum water depth
// This must be done so that water can still move sediment to the edge of the catchment
// and hence remove it from the catchment. (otherwise you would get sediment build
// up around the edges.
//
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
template<typename COORD_MAP>
void Cell::water_flux_out(const COORD_MAP& neighborhood)
{
  // must still figure out best way to accumulate total across cells during update executed by LibGeodecomp
  switch(celltype){
  case Cell::INTERNAL:
    break;
  case Cell::EDGE_WEST:
  case Cell::EDGE_NORTH:
  case Cell::EDGE_EAST:
  case Cell::EDGE_SOUTH:
  case Cell::CORNER_NW:
  case Cell::CORNER_NE:
  case Cell::CORNER_SE:
  case Cell::CORNER_SW:
    if (here_old.water_depth > LSDCatchmentModel::water_depth_erosion_threshold) water_depth = LSDCatchmentModel::water_depth_erosion_threshold;
    break;
  default:
    std::cout << "\n\n WARNING: no water_flux_out rule specified for cell type " << celltype << "\n\n";
    break;
  }
}






// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// THE WATER ROUTING ALGORITHM: LISFLOOD-FP
//
//           X direction
//
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
template<typename COORD_MAP>
void Cell::flow_route_x(const COORD_MAP& neighborhood)
{
  double hflow;
  double tempslope;
  double west_elevation;
  double west_water_depth;
  double local_time_factor = set_local_timefactor();

  switch (celltype){
  case Cell::INTERNAL:
  case Cell::EDGE_NORTH:
  case Cell::EDGE_SOUTH:
    west_elevation = west_old.elevation;
    west_water_depth = west_old.water_depth;
    tempslope = ((west_elevation + west_water_depth) - (here_old.elevation + here_old.water_depth)) / LSDCatchmentModel::DX;
    break;
  case Cell::EDGE_WEST:
  case Cell::CORNER_NW:
  case Cell::CORNER_SW:
    west_elevation = LSDCatchmentModel::no_data_value;
    west_water_depth = 0.0;
    tempslope = LSDCatchmentModel::edgeslope;
    break;
  case Cell::EDGE_EAST:
  case Cell::CORNER_NE:
  case Cell::CORNER_SE:
    west_elevation = west_old.elevation;
    west_water_depth = west_old.water_depth;
    tempslope = LSDCatchmentModel::edgeslope;
    break;
  default:
    std::cout << "\n\n WARNING: no x-direction flow route rule specified for cell type " << celltype << "\n\n";
    break;
  }


  if (here_old.water_depth > 0 || west_water_depth > 0)
    {
      hflow = std::max(here_old.elevation + here_old.water_depth, west_elevation + west_water_depth) - std::max(here_old.elevation, west_elevation);
      if (hflow > LSDCatchmentModel::hflow_threshold)
	{
	  update_qx(neighborhood, hflow, tempslope, local_time_factor);
	  froude_check(qx, hflow);
	  discharge_check(neighborhood, qx, west_water_depth, LSDCatchmentModel::DX);
	}
      else
	{
	  qx = 0.0;
	  // qxs = 0.0;
	}
    }




  // calc velocity now
  //if (qx > 0)
  //  {
  //	vel_dir[7] = qx / hflow;
  //  }
  // refactor: old code tries to update vel_dir belonging to neighbour sites - can't do this because neighborhood is passed as a CONST. Should make flow_route only modify its own vel_dir components, based on neighbouring qx qy (i.e. flip around the current logic of modifying neigbouring vel_dir based on local qx qy)
  // But need all cells to record their hflow (or precompute and store -qx/hflow)
  /*	if (qx < 0)
  	{
  	west_old.vel_dir[3] = (0 - qx) / hflow;
  	}

  	}*/

}



// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// THE WATER ROUTING ALGORITHM: LISFLOOD-FP
//
//           Y direction
//
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
template<typename COORD_MAP>
void Cell::flow_route_y(const COORD_MAP& neighborhood)
{
  double hflow;
  double tempslope;
  double north_elevation;
  double north_water_depth;
  double local_time_factor = set_local_timefactor();

  switch (celltype){
  case Cell::INTERNAL:
  case Cell::EDGE_WEST:
  case Cell::EDGE_EAST:
    north_elevation = north_old.elevation;
    north_water_depth = north_old.water_depth;
    tempslope = ((north_elevation + north_water_depth) - (here_old.elevation + here_old.water_depth)) / LSDCatchmentModel::DY;
    break;
  case Cell::EDGE_NORTH:
  case Cell::CORNER_NW:
  case Cell::CORNER_NE:
    north_elevation = LSDCatchmentModel::no_data_value;
    north_water_depth = 0.0;
    tempslope = LSDCatchmentModel::edgeslope;
    break;
  case Cell::EDGE_SOUTH:
  case Cell::CORNER_SW:
  case Cell::CORNER_SE:
    north_elevation = north_old.elevation;
    north_water_depth = north_old.water_depth;
    tempslope = LSDCatchmentModel::edgeslope;
    break;
  default:
    std::cout << "\n\n WARNING: no y-direction flow route rule specified for cell type " << static_cast<int>(celltype) << "\n\n";
    break;
  }


  if (here_old.water_depth > 0 || north_water_depth > 0)
    {
      hflow = std::max(here_old.elevation + here_old.water_depth, north_elevation + north_water_depth) - std::max(here_old.elevation, north_elevation);
      if (hflow > LSDCatchmentModel::hflow_threshold)
	{
	  update_qy(neighborhood, hflow, tempslope, local_time_factor);
	  froude_check(qy, hflow);
	  discharge_check(neighborhood, qy, north_water_depth, LSDCatchmentModel::DY);
	}
      else
	{
	  qy = 0.0;
	  // qys = 0.0;
	}
    }




  // calc velocity now
  //if (qy > 0)
  //  {
  //vel_dir[1] = qy / hflow;
  // }
  //refactor: old code tries to update vel_dir belonging to neighbour sites - can't do this because neighborhood is passed as a CONST. Should make flow_route only modify its own vel_dir components, based on neighbouring qx qy (i.e. flip around the current logic of modifying neigbouring vel_dir based on local qx qy)
  // But need all cells to record their hflow (or precompute and store -qx/hflow)

  /*		if (qx < 0)
  		{
  		north_old.vel_dir[5] = (0 - qy) / hflow;
  		}

  		}*/
}







// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// DEPTH UPDATE
//
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
template<typename COORD_MAP>
void Cell::depth_update(const COORD_MAP& neighborhood)
{
  double east_qx;
  double south_qy;

  double local_time_factor = set_local_timefactor();

  switch (celltype){
  case Cell::INTERNAL:
  case Cell::EDGE_NORTH:
  case Cell::EDGE_WEST:
  case Cell::CORNER_NW:
    east_qx = east_old.qx;
    south_qy = south_old.qy;
    update_water_depth(neighborhood, east_qx, south_qy, local_time_factor);
    break;
  case Cell::EDGE_EAST:
  case Cell::CORNER_NE:
    east_qx = 0.0;
    south_qy = south_old.qy;
    update_water_depth(neighborhood, east_qx, south_qy, local_time_factor);
    break;
  case Cell::EDGE_SOUTH:
  case Cell::CORNER_SW:
    east_qx = east_old.qx;
    south_qy = 0.0;
    update_water_depth(neighborhood, east_qx, south_qy, local_time_factor);
    break;
  case Cell::CORNER_SE:
    east_qx = 0.0;
    south_qy = 0.0;
    update_water_depth(neighborhood, east_qx, south_qy, local_time_factor);
    break;
  case Cell::NODATA:
    water_depth = 0.0;
    break;
  default:
    std::cout << "\n\n WARNING: no depth update rule specified for cell type " << Cell::CellType(int(celltype))<< "\n\n";
    break;
  }
}

  
template<typename COORD_MAP>
void Cell::update_water_depth(const COORD_MAP& neighborhood, double east_qx, double south_qy, double local_time_factor)
{
  water_depth = 0.001 + here_old.water_depth + local_time_factor * ( (east_qx - here_old.qx)/LSDCatchmentModel::DX + (south_qy - here_old.qy)/LSDCatchmentModel::DY );
}



template<typename COORD_MAP>
void Cell::update_qx(const COORD_MAP& neighborhood, double hflow, double tempslope, double local_time_factor)
{
  update_q(here_old.qx, qx, hflow, tempslope, local_time_factor);
}


template<typename COORD_MAP>
void Cell::update_qy(const COORD_MAP& neighborhood, double hflow, double tempslope, double local_time_factor)
{
  update_q(here_old.qy, qy, hflow, tempslope, local_time_factor);
}


void Cell::update_q(const double &q_old, double &q_new, double hflow, double tempslope, double local_time_factor)
{
  q_new = ((q_old - (Cell::gravity * hflow * local_time_factor * tempslope)) \
	   / (1.0 + Cell::gravity * hflow * local_time_factor * (LSDCatchmentModel::mannings * LSDCatchmentModel::mannings) \
	      * std::abs(q_old) / std::pow(hflow, (10.0 / 3.0))));
}
    

  // FROUDE NUMBER CHECK
  // need to have these lines to stop too much water moving from
  // one cell to another - resulting in negative discharges
  // which causes a large instability to develop
  // - only in steep catchments really
void Cell::froude_check(double &q, double hflow)
{
  if ((std::abs(q / hflow) / std::sqrt(Cell::gravity * hflow)) > LSDCatchmentModel::froude_limit) // correctly reads newly calculated value of q, not here_old.q
    {
      q = std::copysign(hflow * (std::sqrt(Cell::gravity*hflow) * LSDCatchmentModel::froude_limit), q);
    }
}


// DISCHARGE MAGNITUDE/TIMESTEP CHECKS
// If the discharge is too high for this timestep, scale back...
template<typename COORD_MAP>
void Cell::discharge_check(const COORD_MAP& neighborhood, double &q, double neighbour_water_depth, double Delta)
{
  double local_time_factor = set_local_timefactor();
  double criterion_magnitude = std::abs(q * local_time_factor / Delta);

  if (q > 0 && criterion_magnitude > (here_old.water_depth / 4.0))
    {
      q = ((here_old.water_depth * Delta) / 5.0) / local_time_factor;
    }
  else if (q < 0 && criterion_magnitude > (neighbour_water_depth / 4.0))
    {
      q = -((neighbour_water_depth * Delta) / 5.0) / local_time_factor;
    }
}
  



double Cell::set_local_timefactor()
{
  double local_time_factor = LSDCatchmentModel::time_factor;
  if (local_time_factor > (LSDCatchmentModel::courant_number * (LSDCatchmentModel::DX  / std::sqrt(Cell::gravity * (LSDCatchmentModel::maxdepth)))))
    {
      local_time_factor = LSDCatchmentModel::courant_number * (LSDCatchmentModel::DX / std::sqrt(Cell::gravity * (LSDCatchmentModel::maxdepth)));
    }
  return local_time_factor;
}



  // /*  void set_inputoutput_diff()
  // {
  //   input_output_difference = std::abs(waterinput - waterOut);
  //   }*/








class CellInitializer : public LibGeoDecomp::SimpleInitializer<Cell>
{
public:
  using LibGeoDecomp::SimpleInitializer<Cell>::gridDimensions; 
  
  CellInitializer(LSDCatchmentModel *catchment_in, int nsteps) : LibGeoDecomp::SimpleInitializer<Cell>(LibGeoDecomp::Coord<2>(catchment_in->elev.dim2(), catchment_in->elev.dim1()), nsteps)
  {
    catchment = catchment_in;
  }
  
  void grid(LibGeoDecomp::GridBase<Cell, 2> *subgrid)
  {
    Cell::CellType celltype;
    
    for (int x=0; x<=catchment->jmax-1; x++)
      {
	for (int y=0; y<=catchment->imax-1; y++)
	  {
	    if (y == 0)
	      {
		if (x == 0) celltype = Cell::CORNER_NW;
		else if (x == catchment->jmax-1) celltype = Cell::CORNER_NE;
		else celltype = Cell::EDGE_NORTH;
	      }
	    else if (y == catchment->imax-1)
	      {
		if (x == 0) celltype = Cell::CORNER_SW;
		else if (x == catchment->jmax-1) celltype = Cell::CORNER_SE;
		else celltype = Cell::EDGE_SOUTH;
	      }
	    else  
	      {
		if (x == 0) celltype = Cell::EDGE_WEST;
		else if (x == catchment->jmax-1) celltype = Cell::EDGE_EAST;
		else celltype = Cell::INTERNAL;
	      }
	    
	    LibGeoDecomp::Coord<2> coordinate(x, y);
	    subgrid->set(coordinate, Cell(celltype, catchment->elev[y][x], catchment->water_depth[y][x], catchment->qx[y][x], catchment->qx[y][x]));
	  }
      }
  }
private:
  LSDCatchmentModel *catchment;
};






void LSDCatchmentModel::initialise_model_domain_extents()
{
  std::string FILENAME = read_path + "/" + read_fname + "."	\
    + dem_read_extension;
  
  if (!does_file_exist(FILENAME))
    {
      std::cout << "No terrain DEM found by name of: " << FILENAME
		<< std::endl
		<< "You must supply a correct path and filename "
		<< "in the input parameter file" << std::endl;
      exit(EXIT_FAILURE);
    }
  
  try
    {
      std::cout << "\n\nLoading DEM header info, the filename is "
		<< FILENAME << std::endl;
      
      // open the data file
      std::ifstream data_in(FILENAME.c_str());
      
      //Read in raster data
      std::string str;            // a temporary string for discarding text
      
      // read the georeferencing data and metadata
      data_in >> str >> jmax;
      std::cout << "NCols: " << jmax << " str: " << std::endl;
      data_in >> str >> imax;
      std::cout << "NRows: " << imax << " str: " << std::endl;
      data_in >> str >> xll
	      >> str >> yll
	      >> str >> LSDCatchmentModel::DX // cell size or grid resolution
	      >> str >> LSDCatchmentModel::no_data_value;

      LSDCatchmentModel::DY = LSDCatchmentModel::DX;
    }
  catch(...)
    {
      std::cout << "Something is wrong with your initial elevation raster file."
		<< std::endl
		<< "Common causes are: " << std::endl
		<< "1) Data type is not correct"
		<< std::endl << "2) Non standard raster format" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::cout << "The model domain has been set from "
            << "reading the elevation DEM header." << std::endl;
}



void LSDCatchmentModel::check_DEM_edge_condition()
{
  std::cout << "Checking edge cells for suitable catchment outlet point..." << std::endl;
  //check for no_data_values on RH edge of DEM
  double temp = LSDCatchmentModel::no_data_value;
  std::cout << temp << " " << imax << " " << jmax << std::endl;
  
  for (unsigned n = 0; n < jmax; n++)
    {
      // Check bottom edge (row major!)
      if (elev[imax-1][n] > LSDCatchmentModel::no_data_value) temp = elev[imax-1][n];
      // check top edge
      if (elev[1][n] > LSDCatchmentModel::no_data_value) temp = elev[1][n];
    }
  for (unsigned n = 0; n < imax; n++)
    {
      // check LH edge
      if (elev[n][1] > LSDCatchmentModel::no_data_value) temp = elev[n][1];
      // check RH edge
      if (elev[n][jmax-1] > LSDCatchmentModel::no_data_value) temp = elev[n][jmax-1];
    }
  
  if (temp < -10)
    {
      std::cout << "DEM EDGE CONDITION ERROR: LSDCatchmentModel may not function  \
properly, as the edges of the DEM are all nodata (-9999)		\
values. This will prevent any water or sediment from leaving		\
the edge of the model domain (DEM)" << std::endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      std::cout << "Suitable outlet point on edge found. " << std::endl;
    }
}



void LSDCatchmentModel::load_data()
{
  std::string DEM_FILENAME = read_path + "/" + read_fname + "."	+ dem_read_extension;
  
  LSDRaster elevR;
  
  if (!does_file_exist(DEM_FILENAME))
    {
      std::cout << "No terrain DEM found by name of: " << DEM_FILENAME
		<< std::endl
		<< "You must supply a correct path and filename "
		<< "in the input parameter file" << std::endl;
      exit(EXIT_FAILURE);
    }
  
  try
    {
      elevR.read_ascii_raster(DEM_FILENAME);
      elev = elevR.get_RasterData_dbl();
      // Check that there is an outlet for the catchment water
      check_DEM_edge_condition();
    }
  catch(...)
    {
      std::cout << "Something is wrong with your initial elevation raster file."
		<< std::endl
		<< "Common causes are: " << std::endl
		<< "1) Data type is not correct" << std::endl
		<< "2) Non standard raster format" << std::endl;
      exit(EXIT_FAILURE);
    }


  
  // PORT SPATIALLY VARIABLE RAINFALL BY READING HYDROINDEX FILE AND RAINDATA FILE etc.



  

  
}

















//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function gets all the data from a parameter file
//
// Update: It also intialises the other params that are set internally
// Some functions have been taken out of mainloop()
//
// DAV - this is a bit of a clunky method - perhaps replace it with the
//   paramter ingestion method used in CHILD one day?
//
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDCatchmentModel::initialise_variables(std::string pfname)
{
  if(LibGeoDecomp::MPILayer().rank() == 0)
    {
      std::cout << "Initialising the model parameters..." << std::endl;
    }
  // Concatenate the path and parameter file name to get the full file path
  string full_name = pfname;  
  
  std::ifstream infile;
  // Open the parameter file
  infile.open(full_name.c_str());
  string parameter, value, lower, lower_val;
  string bc;
  
  if(LibGeoDecomp::MPILayer().rank() == 0)
    {
      std::cout << "Parameter filename is: " << full_name << std::endl;
    }

  // now ingest parameters
  while (infile.good())
  {
    parse_line(infile, parameter, value);
    lower = parameter;
    if (parameter == "NULL")
      continue;
    for (unsigned int i=0; i<parameter.length(); ++i)
    {
      lower[i] = std::tolower(parameter[i]);  // converts to lowercase
    }

    // get rid of control characters
    value = RemoveControlCharactersFromEndOfString(value);
    
    if (lower == "dem_read_extension")
    {
      dem_read_extension = value;
      dem_read_extension = RemoveControlCharactersFromEndOfString(
                              dem_read_extension);
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "dem_read_extension: " << dem_read_extension << std::endl;
	} 
    }
    else if (lower == "dem_write_extension")
    {
      dem_write_extension = value;
      dem_write_extension = RemoveControlCharactersFromEndOfString(
                              dem_write_extension);
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "dem_write_extension: " << dem_write_extension << std::endl;
	}
    }
    else if (lower == "write_path")
    {
      write_path = value;
      write_path = RemoveControlCharactersFromEndOfString(write_path);
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "output write path: " << write_path << std::endl;
	}
    }
    else if (lower == "write_fname")
    {
      write_fname = value;
      write_fname = RemoveControlCharactersFromEndOfString(write_fname);
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "write_fname: " << write_fname << std::endl;
	}
    }
    else if (lower == "read_path")
    {
      read_path = value;
      read_path = RemoveControlCharactersFromEndOfString(read_path);
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "read_path: " << read_path << std::endl;
	}
    }
    else if (lower == "read_fname")
    {
      read_fname = value;
      read_fname = RemoveControlCharactersFromEndOfString(read_fname);
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "read_fname: " << read_fname << std::endl;
	}
    }
    
    //=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // parameters for landscape numerical modelling
    // (LSDCatchmentModel: DAV 2015-01-14)
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Supplementary Input Files
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "hydroindex_file")
    {
      hydroindex_fname = value;
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "hydroindex_file: " << hydroindex_fname << std::endl;
	}
    }
    else if (lower == "rainfall_data_file")
    {
      rainfall_data_file = value;
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "rainfall_data_file: " << rainfall_data_file << std::endl;
	}
    }
    
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Numerical
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "no_of_iterations")
    {
      no_of_iterations = atoi(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "no of iterations: " << no_of_iterations << std::endl;
	}
    }
    
    
    
    
    // Other numerical parameters, output and save options, sedimentation,
    // erosion, grain size options, still to be ported from original code



    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Hydrology and Flow
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    else if (lower == "hydro_model_only")
    {
      hydro_only = (value == "yes") ? true : false;
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "run hydro model only (NO EROSION): "
		    << hydro_only << std::endl;
	}
    }

    else if (lower == "rainfall_data_on")
    {
      rainfall_data_on = (value == "yes") ? true : false;
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "rainfall_data_on: " << rainfall_data_on << std::endl;
	}
    }

    else if (lower == "topmodel_m_value")
    {
      M = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "topmodel m value: " << M << std::endl;
	}
    }

    else if (lower == "num_unique_rain_cells")
    {
      rfnum = atoi(value.c_str());
      {
	std::cout << "Number of unique rain cells: " << rfnum << std::endl;
      }
    }

    else if (lower == "rain_data_time_step")
    {
      rain_data_time_step = atof(value.c_str());
      {
	std::cout << "rainfall data time step: "
		  << rain_data_time_step << std::endl;
      }
    }

    else if (lower == "spatial_var_rain")
    {
      spatially_var_rainfall = (value == "yes") ? true : false;
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Spatially variable rainfall: "
		    << spatially_var_rainfall << std::endl;
	}
    
    }

    else if (lower == "in_out_difference")
    {
      in_out_difference_allowed = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "in-output difference allowed (cumecs): "
		    << in_out_difference_allowed << std::endl;
	}
    }

    else if (lower == "min_q_for_depth_calc")
    {
      MIN_Q = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "minimum discharge for depth calculation: "
		    << MIN_Q << std::endl;
	}
    }

    else if (lower == "max_q_for_depth_calc")
    {
      MIN_Q_MAXVAL = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "max discharge for depth calc: "
		    << MIN_Q_MAXVAL << std::endl;
	}
    }

    else if (lower == "hflow_threshold")
    {
      LSDCatchmentModel::hflow_threshold = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Horizontal flow threshold: "
		    << LSDCatchmentModel::hflow_threshold << std::endl;
	}
    }

    else if (lower == "water_depth_erosion_threshold")
    {
      LSDCatchmentModel::water_depth_erosion_threshold = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Water depth for erosion threshold: "
		    << LSDCatchmentModel::water_depth_erosion_threshold << std::endl;
	}
    }

    else if (lower == "slope_on_edge_cell")
    {
      LSDCatchmentModel::edgeslope = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Slope on model domain edge: "
		    << LSDCatchmentModel::edgeslope << std::endl;
	}
    }

    else if (lower == "evaporation_rate")
    {
      k_evap = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Evaporation rate: " << k_evap << std::endl;
	}
    }
    
    else if (lower == "courant_number")
    {
      LSDCatchmentModel::courant_number = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Courant number: " << LSDCatchmentModel::courant_number << std::endl;
	}
    }
    
    else if (lower == "froude_num_limit")
    {
      LSDCatchmentModel::froude_limit = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Froude number limit: " << LSDCatchmentModel::froude_limit << std::endl;
	}
    }

    else if (lower == "mannings_n")
    {
      LSDCatchmentModel::mannings = atof(value.c_str());
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "mannings: " << LSDCatchmentModel::mannings << std::endl;
	}
    }
    else if (lower == "spatially_complex_rainfall_on")
    {
      spatially_complex_rainfall = (value == "yes") ? true : false;
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Spatially complex rainfall option: "
		    << spatially_complex_rainfall << std::endl;
	}
    }
    
    
    // Vegetation, hillslope & debugging options still to be ported
    // from original code (master branch)


    //=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // LibGeoDecomp Options
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=

    else if (lower == "simulator")
      {
	LSDCatchmentModel::simulator = value;
      }
    
    
    // Visualisation
    else if (lower == "write_elevation_ppm")
      {
	LSDCatchmentModel::write_elevation_ppm = (value == "yes") ? true : false;
      }
    else if (lower == "write_water_depth_ppm")
      {
	LSDCatchmentModel::write_water_depth_ppm = (value == "yes") ? true : false;
      }
    else if (lower == "water_depth_ppm_interval")
      {
	LSDCatchmentModel::water_depth_ppm_interval = atoi(value.c_str());
      }
    else if (lower == "write_water_depth_bov")
      {
	LSDCatchmentModel::write_water_depth_bov = (value == "yes") ? true : false;
      }
    else if (lower == "water_depth_bov_interval")
      {
	LSDCatchmentModel::water_depth_bov_interval = atoi(value.c_str());
      }
    else if (lower == "write_water_depth_visit")
      {
	LSDCatchmentModel::write_water_depth_visit = (value == "yes") ? true : false;
      }
    else if (lower == "water_depth_visit_interval")
      {
	LSDCatchmentModel::water_depth_visit_interval = atoi(value.c_str());
      }


    
  }

  if(LibGeoDecomp::MPILayer().rank() == 0)
    {
      std::cout << "No other parameters found, parameter ingestion complete."
		<< std::endl;
      std::cout << "Initialising hard coded-constants." << std::endl;
    }
  // Initialise the other parameters (those not set by param file)
  tx = output_file_save_interval;
  
  if (spatially_var_rainfall == false)
    {
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  std::cout << "Making sure no of rain cells is set to 1, "
		    << "for uniform rainfall input.."
		    << std::endl;
	}
      rfnum = 1;
    }
}










void runSimulation(std::string pfname)
{
  
  std::string initial_snapshot_name = "initial_snapshot";
  std::string initial_snapshot_file = initial_snapshot_name + ".mpiio";
  LSDCatchmentModel *catchment = new LSDCatchmentModel(pfname); // so we have model params available on each rank (MPI_Bcast if overhead too large)
  
  // Run one step serially to load input data, initialise grid, and write MPI-IO snapshot
  if (LibGeoDecomp::MPILayer().rank() == 0)
    {
      catchment->initialise_model_domain_extents(); 
      catchment->initialise_arrays(); 
      catchment->load_data();
      CellInitializer *initialiser = new CellInitializer(catchment, 0);
      LibGeoDecomp::SerialSimulator<Cell> init(initialiser);
      
      std::cout << "\n Generating initial MPI-IO snapshot...";
      init.addWriter(new LibGeoDecomp::MPIIOWriter<Cell>("", 1, 1, MPI_COMM_SELF, MPI_CELL));
      init.run();
      std::cout << "done \n";
    }
  
  LibGeoDecomp::MPILayer().barrier();  // Make sure rank 0 has finished writing initial MPI IO snapshot
  LibGeoDecomp::MPIIOInitializer<Cell> *mpiio_initialiser = new LibGeoDecomp::MPIIOInitializer<Cell>("00000.mpiio", MPI_CELL);
  LibGeoDecomp::MPILayer().barrier();
  if( LibGeoDecomp::MPILayer().rank() == 0){ std::cout << "\n MPIIO snapshot loaded, initialising parallel simulation... \n\n"; }


  LibGeoDecomp::DistributedSimulator<Cell> *sim = 0;
  if(catchment->simulator == "striping")
    {
      sim = new LibGeoDecomp::StripingSimulator<Cell>(mpiio_initialiser, LibGeoDecomp::MPILayer().rank()? 0 : new LibGeoDecomp::NoOpBalancer(), 1);
    }
  else if(catchment->simulator == "hipar")
    {
      sim = new LibGeoDecomp::HiParSimulator<Cell, LibGeoDecomp::RecursiveBisectionPartition<2> >(mpiio_initialiser, LibGeoDecomp::MPILayer().rank() ? 0 : new LibGeoDecomp::NoOpBalancer(), 1, 1);
    }
  
  sim->setMaxSteps(catchment->no_of_iterations);
  LibGeoDecomp::MPILayer().barrier();

  
  // Set up visualisation outputs  
  LibGeoDecomp::PPMWriter<Cell> *elevationPPMWriter = 0;
  LibGeoDecomp::PPMWriter<Cell> *water_depthPPMWriter = 0;
  if(catchment->write_elevation_ppm)
    {
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  system("mkdir -p elevation/ppm");
	  elevationPPMWriter = new LibGeoDecomp::PPMWriter<Cell>(&Cell::elevation, 0.0, 255.0, "elevation/ppm/elevation", \
								 catchment->elevation_ppm_interval, LibGeoDecomp::Coord<2>(catchment->pixels_per_cell, catchment->pixels_per_cell));
	}
    }
  if(catchment->write_water_depth_ppm)
    {
      if(LibGeoDecomp::MPILayer().rank() == 0)
	{
	  system("mkdir -p water_depth/ppm");
	  water_depthPPMWriter = new LibGeoDecomp::PPMWriter<Cell>(&Cell::water_depth, 0.0, 1.0, "water_depth/ppm/water_depth", \
								   catchment->water_depth_ppm_interval, LibGeoDecomp::Coord<2>(catchment->pixels_per_cell, catchment->pixels_per_cell));
	}
      LibGeoDecomp::CollectingWriter<Cell> *water_depthPPMCollectingWriter = new LibGeoDecomp::CollectingWriter<Cell>(water_depthPPMWriter);
      sim->addWriter(water_depthPPMCollectingWriter);
    }
  if(catchment->write_water_depth_bov)
    {
      system("mkdir -p water_depth/bov");
      sim->addWriter(new LibGeoDecomp::BOVWriter<Cell>(LibGeoDecomp::Selector<Cell>(&Cell::water_depth, "water_depth"), "water_depth/bov/water_depth", \
						      catchment->water_depth_bov_interval));
    }
  
  
  if (LibGeoDecomp::MPILayer().rank() == 0){ sim->addWriter(new LibGeoDecomp::TracingWriter<Cell>(1, catchment->no_of_iterations)); }
  
  sim->run();
      
  LibGeoDecomp::MPILayer().barrier(); 
  
}

