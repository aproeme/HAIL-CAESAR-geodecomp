#CXX := g++-9 # This is the main compiler
CXX := mpic++

#GEODECOMP_DIR := /Users/aproeme/libgeodecomp/9529db
GEODECOMP_DIR := /Users/aproeme/libgeodecomp/hail-caesar
GEODECOMP_SRC_DIR := /Users/aproeme/libgeodecomp/libgeodecomp
#GEODECOMP_SRC_DIR := /Users/aproeme/libgeodecomp/libgeodecomp-0.4.0

BOOST_DIR := /usr/local/Cellar/boost/1.71.0
MPI_DIR := /usr/local/Cellar/mpich/3.3.2

# CC := clang --analyze # and comment out the linker last line for sanity
GITREV = -D'GIT_REVISION="$(shell git log --pretty=format:'%h' -n 1)"'
SRCDIR := src
BUILDDIR := build

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -fopenmp -std=c++11 $(GITREV) -Wfatal-errors #-O3 or -O2  REMEMBER TO TURN THIS BACK ON BEFORE RUNNING PRODUCTION
INC := -I . -I ./include -I./include/libgeodecomp -I $(GEODECOMP_DIR)/include -I $(BOOST_DIR)/include  
LDFLAGS := -fopenmp -L $(GEODECOMP_DIR)/lib -L $(BOOST_DIR)/lib 
LIBS := -lgeodecomp -lboost_date_time 

TYPEMAP_TEST_OBJECTS := src/catchmentmodel/LSDCatchmentModel.o src/libgeodecomp/typemaps.o test/typemaptest.o

ifeq ($(SERIAL),TRUE)
  build_mode := serial
  TARGET := bin/HAIL-CAESAR.serial
  CFLAGS += -DCOMPILE_FOR_SERIAL
else
  build_mode := parallel
  TARGET := bin/HAIL-CAESAR.mpi
  CFLAGS += -DCOMPILE_FOR_PARALLEL
  INC += -I $(MPI_DIR)/include 
  LDFLAGS += -L $(MPI_DIR)/lib
  LIBS += -lmpi
endif




$(TARGET): $(OBJECTS) 
	@echo -e " \n Linking... \n"
	@echo " $(CXX) $(LDFLAGS) $(INC) $(LIBS) $^ -o $(TARGET)"; $(CXX) $(LDFLAGS) $(INC) $(LIBS) $^ -o $(TARGET) 

typemaptest: $(TYPEMAP_TEST_OBJECTS) $(OBJECTS)
	@echo -e " \n Linking... \n"
	@echo " $(CXX) $(LDFLAGS) $(INC) $(LIBS) $(TYPEMAP_SRC_OBJECTS) -o bin/typemaptest"; $(CXX) $(LDFLAGS) $(INC) $(LIBS) $(TYPEMAP_TEST_OBJECTS) -o bin/typemaptest

typemaptest.o : test/typemaptest.cpp include/catchmentmodel/LSDCatchmentModel.hpp include/libgeodecomp/typemaps.h
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o test/typemaptest.o test/typemaptest.cpp"; $(CXX) $(CFLAGS) $(INC) -c -o test/typemaptest.o test/typemaptest.cpp

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) 
	@mkdir -p bin
	@mkdir -p $(BUILDDIR)/topotools
	@mkdir -p $(BUILDDIR)/catchmentmodel
	@mkdir -p $(BUILDDIR)/libgeodecomp
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

typemaps: # only necessary to generate new MPI typemaps run if the model has been changed, not needed during normal build process
	@echo " Generating xml using doxygen..."; echo "doxygen doxygen.conf"; doxygen doxygen.conf
	@echo " Generating MPI typemaps...";
	@mkdir typemaps
	@$(GEODECOMP_SRC_DIR)/tools/typemapgenerator/generate.rb -S typemaps-doxygen-docs/xml typemaps
	@sed -i '' 's/CellType/int/g' typemaps/typemaps.cpp  # correct mistake in LibGeoDecomp typemap generation script for enums, see https://github.com/gentryx/libgeodecomp/issues/72
	@cp typemaps/*.h include/libgeodecomp/
	@cp typemaps/*.cpp src/libgeodecomp/
	@echo " Typemaps saved in src/libgeodecomp and include/libgeodecomp"
	@rm -rf typemaps
	@rm -rf typemaps-doxygen-docs

print_build_mode:
	@echo -e "\n Building $(build_mode) version of HAIL-CAESAR \n" 

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -rf $(BUILDDIR) $(TARGET) typemaps typemaps-doxygen-docs"; $(RM) -r $(BUILDDIR) $(TARGET) typemaps typemaps-doxygen-docs

# Tests
tester:
	$(CXX) $(CFLAGS) test/tester.cpp $(INC) $(LIBS) -o bin/tester

# Spikes
ticket:
	$(CXX) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIBS) -o bin/ticket

.PHONY: clean
