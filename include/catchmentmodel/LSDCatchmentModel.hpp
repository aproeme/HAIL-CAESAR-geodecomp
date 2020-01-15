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

#ifdef COMPILE_FOR_PARALLEL
#include <typemaps.h>
#endif

#ifndef LSDCatchmentModel_geodecomp_H
#define LSDCatchmentModel_geodecomp_H


void runSimulation(std::string pname, std::string pfname);



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
  /// @param pname the pathname to the parameter file
  /// @param pfname the filename of the parameter file !!INCLUDING EXTENSION!!
  /// @author DAV
  /// @date 2015-01-16
  LSDCatchmentModel(string pname, string pfname)
  {
    create(pname, pfname);
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
  void initialise_variables(std::string pname, std::string pfname);

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

  /// @brief Calculates which cells have water content and marks these
  /// cells with a flag.
  /// @detail By creating a mask of cells that contain water, later calculations
  /// can be sped up by only iterating over cells that contain water (and there
  /// -fore will have discharges, water to route, erosion occuring etc.)
  /// Interesting to see how useful this is in catchments where majority of
  /// cells contain water, i.e. during flood events, etc. and whether there
  /// is worth to this method which has overhead but ultimately reduces
  /// compute time, especially in situations where there is large areas of
  /// catchment with minimal water content.
  void scan_area();

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
  
  int no_of_iterations = 3;

  // visualisation options
  bool write_elevation_ppm = false;
  bool write_elevation_bov = false;
  bool write_elevation_visit = false;
  bool write_water_depth_ppm = false;
  bool write_water_depth_bov = false;
  bool write_water_depth_visit = false;
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
  
  bool uniquefilecheck = false;

  //constants
  const double root = 7.07;
  const double gravity = 9.81;
  const float g = 9.81F;
  const float kappa = 0.4F;
  const int ACTIVE_FACTOR=1;
  const int TRUE=1;
  const int FALSE=0;
  const unsigned int G_MAX=10;
  const std::array<int, 9> deltaX = {{0,  0,  1,  1,  1,  0, -1, -1, -1}};
  const std::array<int, 9> deltaY = {{0, -1, -1,  0,  1,  1,  1,  0, -1}};

  static double water_depth_erosion_threshold;
  int input_time_step = 60;
  int number_of_points = 0;
  double globalsediq = 0;
  double time_1 = 1;
  double save_time = 0;
  double creep_time = 1;
  double creep_time2 = 1;
  double creep_time3 = 1;
  double grass_grow_interval = 1;

  double soil_erosion_time = 1;
  double soil_development_time = 1;

  double bedrock_erosion_threshold = 0;
  double bedrock_erosion_rate = 0;
  double p_b = 1.5;  // detach capacity exponent from CHILD
  double bedrock_erodibility_coeff_ke = 0.002;

  ///int tot_number_of_tracer_points=0;
  int input_type_flag=0; // 0 is water input from points, 1 is input from hydrograph or rainfall file.
  double failureangle=45;
  double saveinterval=1000;
  int counter=0;

  double waterinput = 0;
  double waterOut = 0;
  double input_output_difference = 0;
  double in_out_difference_allowed = 0;
  static double mannings;

  /// no. of rainfall cells
  unsigned rfnum = 1;

  /// set by ncols and nrows
  unsigned int jmax, imax;
  double xll, yll;
  static double no_data_value;

  int maxcycle = 1000;

  double ERODEFACTOR=0.05;
  static double DX;
  static double DY;
  /// memory limit
  int LIMIT=1;
  double MIN_Q=0.01;
  double MIN_Q_MAXVAL=1000.0;
  double CREEP_RATE=0.0025;
  double SOIL_RATE = 0.0025;
  double active=0.2;
  int grain_array_tot =1 ;

  /// Number of passes for edge smoothing filter
  double edge_smoothing_passes = 100.0;
  /// Number of cells to shift lat erosion downstream
  double downstream_shift= 5.0;
  /// Max difference allowed in cross channel smoothing of edge values
  double lateral_cross_channel_smoothing = 0.0001;
  double lateral_constant=0.0000002;

  static double time_factor;
  std::vector<double> j, jo, j_mean, old_j_mean, new_j_mean;

  /// TOPMODEL 'm'
  double M = 0.005;
  double baseflow = 0.00000005; //end of hyd model variables usually 0.0000005 changed 2/11/05
  // Reverted to match CL 1.8f 17/08/16 - DV

  double cycle =0;  // can't initalise static vars in header file!
  double rain_factor = 1;
  double sediQ = 0;

  // speed in which vegetation reaches full maturity in year
  double grow_grass_time = 0;
  double duneupdatetime = 0;

  double output_file_save_interval = 60;
  double min_time_step = 0;
  double vegTauCrit = 100;

  int max_time_step = 0;
  int dune_mult = 5;
  double dune_time = 1;
  double max_vel = 5;
  double sand_out = 0;
  static double maxdepth;
  static double courant_number;
  double erode_call = 0;
  double erode_mult = 1;
  double lateralcounter = 1;
  static double edgeslope;
  double chann_lateral_erosion = 20.0; // formerly 'bed_proportion =0.01'
  double veg_lat_restriction = 0.1;

  static double froude_limit;
  double recirculate_proportion = 1;

  double Csuspmax = 0.05; // max concentration  of SS allowed in a cell (proportion)
  static double hflow_threshold;

  // KAtharine
  int variable_m_value_flag = 0;

  // TO DO: DAV - these could be read from an input file.
  // Swale grainsizes
  double d1=0.000065;
  double d2=0.001;
  double d3=0.002;
  double d4=0.004;
  double d5=0.008;
  double d6=0.016;
  double d7=0.032;
  double d8=0.064;
  double d9=0.128;

  // std::array<double, 11> dprop = {{0.0, 0.144, 0.022, 0.019, 0.029, 0.068, 0.146, 0.22, 0.231, 0.121, 0.0}}; // Default
  double dprop[11] = {0.0, 0.05, 0.05, 0.15, 0.225, 0.25, 0.1, 0.075, 0.05, 0.05, 0.0}; // Swale

  double previous;
  int hours = 0;
  double new_cycle = 0;
  double old_cycle = 0;
  double tx = 60;
  double Tx = 0;
  double tlastcalc = 0;
  double Qs_step = 0;
  double Qs_hour = 0;
  double Qs_over = 0;
  double Qs_last = 0;
  double Qw_newvol = 0;
  double Qw_oldvol = 0;
  double Qw_lastvol = 0;
  double Qw_stepvol = 0;
  double Qw_hourvol = 0;
  double Qw_hour = 0;
  double Qw_overvol = 0;
  double temptotal = 0;
  double old_sediq = 0;

  double temptot = 0;

  std::vector<double> sum_grain, sum_grain2;
  std::vector<double> old_sum_grain,old_sum_grain2;
  std::vector<double> Qg_step, Qg_step2, Qg_hour, Qg_hour2;
  std::vector<double> Qg_over, Qg_over2, Qg_last,Qg_last2;

  TNT::Array2D<double> elev;
  TNT::Array2D<double> bedrock;
  TNT::Array2D<double> init_elevs;
  TNT::Array2D<double> water_depth;
  TNT::Array2D<double> area;
  TNT::Array2D<double> tempcreep;
  TNT::Array2D<double> Tau;
  TNT::Array2D<double> Vel;
  TNT::Array2D<double> qx;
  TNT::Array2D<double> qy;
  TNT::Array2D<double> qxs;
  TNT::Array2D<double> qys;
  // TODO - these are for the dune model which is as of yet unimplemented
  TNT::Array2D<double> area_depth;
  TNT::Array2D<double> sand;
  TNT::Array2D<double> elev2;
  TNT::Array2D<double> sand2;
  TNT::Array2D<double> grain;
  TNT::Array2D<double> elev_diff;

  TNT::Array2D<int> index;
  TNT::Array2D<int> down_scan;
  TNT::Array2D<int> rfarea;

  TNT::Array2D<bool> inputpointsarray;

  std::vector<int> catchment_input_x_coord;
  std::vector<int> catchment_input_y_coord;

  TNT::Array3D<double> vel_dir;
  TNT::Array3D<double> strata;

  std::vector<double> hourly_m_value;
  std::vector<double> temp_grain;
  std::vector< std::vector<float> > hourly_rain_data;
  TNT::Array3D<double> veg;
  TNT::Array2D<double> edge, edge2; //TJC 27/1/05 array for edges
  std::vector<double> old_j_mean_store;
  TNT::Array3D<double> sr, sl, su, sd;
  TNT::Array2D<double> ss;

  // MJ global vars
  std::vector<double> fallVelocity;
  std::vector<bool> isSuspended;
  TNT::Array2D<double> Vsusptot;


  std::vector<int> nActualGridCells;
  double Jw_newvol = 0.0;
  double Jw_oldvol = 0.0;
  double Jw_lastvol = 0.0;
  double Jw_stepvol = 0.0;
  double Jw_hourvol = 0.0;
  double Jw_hour = 0.0;
  double Jw_overvol = 0.0;
  double k_evap = 0.0;

  // sedi tpt flags
  bool einstein = false;
  bool wilcock = false;
  int div_inputs = 1;
  double rain_data_time_step = 60; // time step for rain data - default is 60.

  // lisflood caesar adaptation globals
  std::vector<int> catchment_input_counter;
  int totalinputpoints = 0;

  /// Soil generation variables
  double P1, b1, k1, c1, c2, k2, c3, c4;

  /// Option Bools
  bool soildevoption = false;
  bool suspended_opt = false;
  bool jmeaninputfile_opt = false;
  // This is for reading discharge direct from input file - DAV 2/2016
  bool recirculate_opt = false;
  bool reach_mode_opt = false;
  bool dunes_opt = false;
  bool bedrock_lower_opt = false;
  bool physical_weather_opt = false;
  bool chem_weath_opt = false;
  bool allow_in_out_diff = true;

  bool soil_j_mean_depends_on = false;
  bool rainfall_data_on = false;
  bool hydro_only = false;
  bool vegetation_on = false;
  bool bedrock_layer_on = false;
  bool lateral_erosion_on = false;
  bool spatially_var_rainfall = false;
  bool graindata_from_file = false;

  bool spatially_complex_rainfall = false;

  int erode_timestep_type = 0;  // 0 for default based on erosion amount, 1 for basedon hydro timestep
  int hydro_timestep_type = 0;  // 0 for default

  // Bools for writing out files
  bool write_elev_file = false;
  bool write_grainsz_file = false;
  bool write_params_file = false;
  bool write_flowvel_file = false;
  bool write_waterd_file = false;
  bool write_elevdiff_file = false;

  /// input file names
  std::string rainfall_data_file;
  std::string grain_data_file;
  std::string bedrock_data_file;

  /// output file names
  std::string elev_fname;
  std::string grainsize_fname;
  std::string params_fname;
  std::string waterdepth_fname;
  std::string flowvel_fname;
  std::string hydroindex_fname;
  std::string timeseries_fname;
  std::string elevdiff_fname;
  std::string raingrid_fname;
  std::string runoffgrid_fname;
  
  bool DEBUG_print_cycle_on = false;
  bool DEBUG_write_raingrid = false;
  bool DEBUG_write_runoffgrid = false;

  int timeseries_interval;
  float run_time_start;

  int tempcycle = 0;
  
  
  std::vector< std::vector<float> > raingrid;	 // this is for the rainfall data file

  // Mainly just the definitions of the create() functions go here:
  // The implementations are in the .cpp file.

  void create();
  void create(std::string pname, std::string pfname);
};

#endif


