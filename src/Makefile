### Makefile for the package
### Author: Yan Y. Liu <yanliu@illinois.edu>
### Date: 06/18/2014
# view what happens when make runs
SHELL       += -x
# compiler
CXX          = mpicxx
# compile options
CXXFLGS     += -g -Wall -std=c++0x 
# openmp options
# intel compiler
OPENMP= -qopenmp
# gcc
#OPENMP=-fopenmp -lgomp
# source file pattern
SRCC         = cc
SRCH         = h

# SPTW
#SPTW_HOME    ?= ../../sptw
SPTW_HOME    ?= $(WORK)/tutorial/sptw
GDALTIFF_DIR ?= $(WORK)/software/gdal-1.11.0/frmts/gtiff/libtiff

# functions and classes
TARGETS      = data util #$(SPTW_HOME)/sptw 
OBJS         = $(TARGETS:=.o)
# executables
MAINS        = mpibasic sendrecv mapalg-diff mapalg-dist mapalg-dist-openmp mapalg-dist-mrsw #mapalg-dist-mrmw

# external libs
EXTDIRS      = $(GDAL_HOME)
INCPATH      = $(EXTDIRS:%=-I%/include)
LIBPATH      = $(EXTDIRS:%=-L%/lib)
LIBS         = -lgdal

# extra
EXTRAINC     = #-I$(GDALTIFF_DIR) -I$(SPTW_HOME) 

# default system libs
LIBS_DEFAULT = -lm

all: $(MAINS) Makefile

# compile and link executables
$(MAINS): % : %.$(SRCC) $(DEFS) $(OBJS)
	@$(CXX) $(CXXFLGS) $(OPENMP) -I. $(INCPATH) $(EXTRAINC) -o $@ $< $(OBJS) $(LIBPATH) $(LIBS) $(LIBS_DEFAULT)

# compile functions
$(OBJS): %.o: %.$(SRCC) %.$(SRCH)
	@$(CXX) $(CXXFLGS) $(OPENMP) -I. $(INCPATH) $(EXTRAINC) -o $@ -c $<

# clean
clean:
	@rm -f $(OBJS) $(MAINS)
