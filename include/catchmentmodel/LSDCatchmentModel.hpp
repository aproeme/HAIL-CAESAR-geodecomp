// LSDCatchmentModel.hpp
//
// Header file for the LSDCatchmentModel

#include <vector>
#include <cmath>
#include <string>
#include <array>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>

// Include for OpenMP
#include <omp.h>

// Topotools includes for stats tools and raster handling
#include "topotools/LSDRaster.hpp"
#include "topotools/LSDStatsTools.hpp"

#include "LSDGrainMatrix.hpp"
#include "LSDRainfallRunoff.hpp"

#include "TNT/tnt.h"   // Template Numerical Toolkit library

#include <libgeodecomp/misc/apitraits.h>
#include <libgeodecomp/geometry/coord.h>

#include <typemaps.h>

#ifndef LSDCatchmentModel_geodecomp_H
#define LSDCatchmentModel_geodecomp_H


void runSimulation(std::string pfname);



class LSDCatchmentModel: public LSDRaster
{
  friend class Cell;
  friend class CellInitializer;
  friend class Typemaps;
  
public:

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // Constructors
  //
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  /// @brief The default constructor.
  /// @details this asks for a pathname and a filename of the parameter file
  /// It then opens the paramter file and ingests the information
  /// @author DAV
  /// @date 2015-01-16
  LSDCatchmentModel()
  {
    create();
  }

  /// @brief this constructor just reads the param file given by the path and
  /// filename. You must give the parameter file extension!
  /// @param pfname the filename of the parameter file !!INCLUDING EXTENSION!!
  /// @author DAV
  /// @date 2015-01-16
  LSDCatchmentModel(string pfname)
  {
    create(pfname);
  }

  
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // MODEL DOMAIN METHODS
  // Set up of model domain, initialisation
  // of arrays.
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  /// @brief Initialises the size of the arrays holding the various
  /// model fields such as elevation, water depth etc.
  void initialise_model_domain_extents();

  /// @brief Checks that there is a real terrain point on at least one side of the DEM
  /// and also counts the number of actual grid cells in the catchment.
  /// @details This only currently checks for an edge that is not NODATA on at least one side
  /// It does not check that the DEM has its lowest point on this edge. This
  /// should probably be added.
  void check_DEM_edge_condition();

  /// @brief reads data values from the parameter file into the relevant maps
  /// @return
  void initialise_variables(std::string pfname);

  /// @brief initialises array sizes based on DEM dimensions
  /// @details also sets 'hard-coded' parameters to start the model
  void initialise_arrays();
  
  int get_imax() const { return imax; }
  int get_jmax() const { return jmax; }

  // MODEL OPERATION FLAGS

  /// @brief Is this a hydrology only simulation?
  /// I.e. no erosion methods.
  bool is_hydro_only() const { return hydro_only; }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // INPUT/OUTPUT
  // Methods for loading and manipulating files
  // (should probably go in separate class/file)
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  /// @brief Loads the required data files based on the parameters set in the parameter file
  /// @author dav
  void load_data();

  /// Prints the initial values set by user from param file
  /// as well as those default initial values in the code.
  void print_parameters();

  /// @brief Loads the rainfall data file which is in a special format (headerless text file)
  /// @author DAV
  /// @details Rainfall data file is not too big, so think it's okay to use the vector<vector>
  /// method instead of TNT array. Easier to dynamically resize as the rainfall data contains
  /// no header and can vary in size. Saves the user having to count the rows and cols. Reads
  /// in the rainfall data file specified in the parameter list file as floats.
  /// @return Returns a vector of vector<float>. (A 2D-like vector).
  std::vector< std::vector<float> > read_rainfalldata(std::string FILENAME);

    /// @brief Reads in the grain data file, note, that this is not a raster and in
  /// a special format like the rainfall file.
  /// @author DAV
  /// @details The grain data file consists of multiple columns of data in a text file
  /// that store the grain size fractions for the surface and the subsurface strata. Also
  /// contains an index number and the x y location of each grid cell containing grain data.
  void ingest_graindata_from_file(std::string FILENAME);

  /// @brief Prints the contents of the rainfall data for checking
  /// @author DAV
  /// @details Uses iterators <iterator> header to iterate through
  /// the vector of vectors that is raingrid.
  /// @details Actually, this is a generic function for printing a 2D vector
  void print_rainfall_data();

  /// @brief Calls the various save functions depending on the data types to be saved (Raster Output)
  /// @author DAV
  /// @details dependent on the LSDRaster class calling the overloaded write_raster func. If you are looking
  /// for the function that writes the hydrograph/sediment time series, see the write_output() function.
  void save_raster_data(double tempcycle);

  /// @brief Writes the timeseries file for current timestep.
  /// @detail Writes discharge and sediment flux according to
  /// the same format as found in the CAESAR-Lisflood catchmetn
  void write_output_timeseries(runoffGrid& runoff);

  void save_raster_output();

  /// @brief Writes the time series of catchment output data.
  void output_data();

  /// @brief Writes the time series of catchment output data.
  /// @details Overloaded to take a reference to a runoff object
  /// to allow calculation from the OOP method.-
  /// @author DAV
  void output_data(double temptotal, runoffGrid& runoff);


  // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // MODEL TIMING CONTROL
  // Routines for incrementing counters,
  // setting the timestep etc.
  // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  void set_time_counters();

  void set_loop_cycle();

  void set_global_timefactor();

  double set_local_timefactor();

  void increment_counters();

  void print_cycle();

  double get_cycle() const { return cycle; }
  int get_maxcycle() const { return maxcycle; }

  /// @brief Zeros certain arrays which have to be reset every timestep
  /// or every certain number of timesteps.
  void zero_values();

  
  // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // CATCHMENT MORPHOMETRICS
  // Methods that return information about
  // catchment size, drainage area, wetted area
  // etc.
  // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  /// @brief Counts the number of cells within the catchment boundary.
  /// @details For
  /// spatially variable rainfall, this counts the number of cells
  /// within each hydroindex region that are within the boundary. It
  /// modifies nActualGridCells.
  void count_catchment_gridcells();

  /// @brief Calculate which cells within the catchment are
  /// underwater.
  void check_wetted_area(int scan_area_interval_iter);

  void initialise_drainage_area();

  /// @brief Calculates the area and calls area4().
  /// @details This is basically a wrapper function now. It sets area = 0,
  ///  and area_depth = 1 where there is actual elevation data. Then calls
  /// get_area4() which does the actual work.
  /// @author Translated by DAV
  void zero_and_calc_drainage_area();

  /// @brief Calculates the drainage area, after being called by get_area()
  /// @author Translated by DAV
  void drainage_area_D8();


  // Skipped erosion & slope process components - still to be ported



  
  // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // HYDROLOGY COMPONENTS
  // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  /// @brief Initialises the rainfall runoff grid if using
  /// spatially complex rainfall runoff object.
  void initialise_rainfall_runoff(runoffGrid& runoff);

  /// @brief Updates the water depths (and susp sedi concentrations)
  void depth_update();

  /// @brief Performs the water routing method
  /// @details Uses the Bates et al. (2010) simplification of the St. Venant
  /// shallow water equations. Does not assume steady state flow across the
  /// landscape. (i.e. this is the LISFLOOD algorithm)
  void flow_route();
  
  /// @brief Wrapper that determines which water input routine to call,
  /// either the default one, or the object-based one with spatially complex
  /// runoff paterns.
  /// @todo This logic needs simplifying, why bother creating a runoff object
  /// if it is never used for the simple runoff case (which is most uses.)
  //  void catchment_waterinputs(runoffGrid& runoff);
  
  /// @brief Calculates the amount of runoff on a grid-cell from the rainfall
  /// timeseries input
  //  void catchment_water_input_and_hydrology( double local_time_factor);
  
  /// @brief Overloaded function is for when using the fully distriuted/complex
  /// rainfall patterns option in the model. Takes a reference to the runoffGrid
  /// object.
  //void catchment_water_input_and_hydrology( double local_time_factor, runoffGrid& runoff);

  /// @brief Gets the number of catchment cells that have water input to them
  /// @detail Calculates which cells contain a discharge greater than MIN_Q
  /// and lower than MIN_Q_MAXVAL multiplies by a parameter related to the
  /// contributng draingage area and the BASEFLOW parameter.
  void get_catchment_input_points();

  /// @brief Same as method above but uses the runoff object-based approach
  /// @detail This version is still experimental as of 2017 -DAV TODO
  void get_catchment_input_points(runoffGrid& runoff);

  /// Calculates the amount of water entering grid cells from the rainfall timeseries
  /// and hydroindex if spatially variable rainfall is used.
  /// @details Based on the semi-dsitributed TOPMODEL rainfall runoff model
  void topmodel_runoff( double cycle);

  /// Calculates amount of water entering grid cells when using the fully-distributed
  /// rainfall runoff model (i.e. where every single cell can have diffferent rainfall
  /// and saturation levels.
  /// @details Based on TOPMODEL, modified to fully 2D distributed version
  void topmodel_runoff(double cycle, runoffGrid& runoff);

  /// @brief Calculates the hydrograph values (TOPMODEL) for printing to
  /// the output timeseries file.
  void calchydrograph( double time);

  /// @brief Overloaded function for calculating hydrograph when using the fully distributed
  /// model. Takes an extra reference to the runoff object.
  void calchydrograph(double time, runoffGrid& runoff);

  /// @brief Evaporation routine.
  void evaporate(double time);


  /// @brief Calculates water exiting from the catchment boundaries (on all
  /// four sides of the domain, regardles of where 'true' catchment outlet
  /// point is.)
  void water_flux_out();

  /// @brief Calculates the difference between water entering the catchment
  /// and water leaving the catchment.
  /// If this value is below a user-set threshold, the timestep can be increased
  /// during periods of low water flow. (e.g. inter-storm periods.)
  void set_inputoutput_diff();


  
  // skipped vegetation growth, see original code

  
  TNT::Array2D<double> elev;
  TNT::Array2D<double> water_depth;
  
  // simulation options
  int no_of_iterations = 100;
  std::string simulator = "hipar";
    
  /// set by ncols and nrows
  unsigned int jmax, imax;
  double xll, yll;
  static double no_data_value;
  
  // visualisation options
  bool elevation_ppm = false;
  bool elevation_bov = false;
  bool elevation_visit = false;
  bool water_depth_ppm = false;
  bool water_depth_bov = false;
  bool water_depth_visit = false;
  int elevation_ppm_interval = 1;
  int elevation_bov_interval = 1; 
  int elevation_visit_interval = 1;
  int water_depth_ppm_interval = 1;
  int water_depth_bov_interval = 1;
  int water_depth_visit_interval = 1;
  int pixels_per_cell = 10;

  string dem_read_extension;
  string dem_write_extension;
  string write_path;
  string read_path;
  string write_fname = "catchment.dat";
  string read_fname;

  private:
  
  //constants
  const double root = 7.07;
  const double gravity = 9.81;
  const float g = 9.81F;
  
  static double water_depth_erosion_threshold;
  double time_1 = 1;
  double waterinput = 0;
  double waterOut = 0;
  double input_output_difference = 0;
  double in_out_difference_allowed = 0;
  static double mannings;

  /// no. of rainfall cells
  unsigned rfnum = 1;

  int maxcycle = 1000;

  double ERODEFACTOR=0.05;
  static double DX;
  static double DY;

  static double time_factor;
  std::vector<double> j, jo, j_mean, old_j_mean, new_j_mean;

  /// TOPMODEL 'm'
  double M = 0.005;
  double baseflow = 0.00000005; //end of hyd model variables usually 0.0000005 changed 2/11/05
  // Reverted to match CL 1.8f 17/08/16 - DV

  double cycle = 0; 
  double output_file_save_interval = 60;
  int max_time_step = 0;
  
  static double maxdepth;
  static double courant_number;
  static double edgeslope;
  static double froude_limit;
  static double hflow_threshold;

  double tx = 60;

  TNT::Array2D<int> rfarea;


  std::vector<int> catchment_input_x_coord;
  std::vector<int> catchment_input_y_coord;

  TNT::Array1D<double> vel_dir;

  std::vector<double> hourly_m_value;
  std::vector< std::vector<float> > hourly_rain_data;
  std::vector<double> old_j_mean_store;


  std::vector<int> nActualGridCells;
  double Jw_newvol = 0.0;
  double Jw_oldvol = 0.0;
  double Jw_lastvol = 0.0;
  double Jw_stepvol = 0.0;
  double Jw_hourvol = 0.0;
  double Jw_hour = 0.0;
  double Jw_overvol = 0.0;
  double k_evap = 0.0;

  double rain_data_time_step = 60; // time step for rain data - default is 60.

  // lisflood caesar adaptation globals
  std::vector<int> catchment_input_counter;
  int totalinputpoints = 0;

  /// Option Bools
  bool jmeaninputfile_opt = false;
  // This is for reading discharge direct from input file - DAV 2/2016
  bool recirculate_opt = false;
  bool reach_mode_opt = false;
  bool allow_in_out_diff = true;

  bool rainfall_data_on = false;
  bool hydro_only = false;
  bool spatially_var_rainfall = false;
  bool spatially_complex_rainfall = false;

  // Bools for writing out files
  bool write_elev_file = false;
  bool write_params_file = false;
  bool write_flowvel_file = false;
  bool write_waterd_file = false;
  bool write_elevdiff_file = false;

  /// input file names
  std::string rainfall_data_file;

  /// output file names
  std::string elev_fname;
  std::string params_fname;
  std::string waterdepth_fname;
  std::string flowvel_fname;
  std::string hydroindex_fname;
  std::string timeseries_fname;
  std::string elevdiff_fname;
  std::string runoffgrid_fname;
  
  std::vector< std::vector<float> > raingrid;	 // this is for the rainfall data file

  // Mainly just the definitions of the create() functions go here:
  // The implementations are in the .cpp file.

  void create();
  void create(std::string pfname);
};

#endif


