# This makefile is used to compile ASTR code.
# The compiler: gfortran compiler
#
#FCFLAGS=  -Wuse-without-only -g
#FC=mpif90
FC=h5pfc
#FC=ftn

SRCDIR = src user_define_module
OBJDIR = obj
BINDIR = bin
CTRDIR = cma

ifeq ($(USER),lthpc)
    export FFTW_INC_DIR=/usr/include
	export FFTW_LIB_DIR=/usr/lib/x86_64-linux-gnu
endif



FCFLAGS= -O3 -fbounds-check -I$(FFTW_INC_DIR) -lfftw3

# OPTIONS1 = -fcheck=all
OPTIONS2 = -J $(OBJDIR)
OPTIONS3 = -DHDF5
# OPTIONS4 = -DCOMB -I$(CTRDIR)/include/cantera 
# OMP = -fopenacc


EXE=astr

LIBS= -lz -lm -L$(FFTW_LIB_DIR) -lfftw3_mpi -lfftw3# -L$(CTRDIR)/lib -lcantera_fortran -lcantera -lstdc++ -pthread 
#LIBS= -lz -lm 

TARGET = $(BINDIR)/$(EXE)

VPATH = $(SRCDIR):$(OBJDIR)

srs=  strings.F90 fdnn.F90 singleton.F90 commtype.F90 stlaio.F90 constdef.F90 tecio.F90     \
      vtkio.F90 interp.F90 commvar.F90 utility.F90 thermchem.F90 commarray.F90 fludyna.F90  \
      parallel.F90 fftwlink.F90 hdf5io.F90 cmdefne.F90 commfunc.F90 commcal.F90 models.F90  \
      statistic.F90 comsolver.F90 \
	  udf_tool.F90 userdefine_2dhit.F90 \
	  bc.F90 readwrite.F90 geom.F90 ibmethod.F90\
	  gridgeneration.F90 riemann.F90 solver.F90 \
	  udf_pp.F90 udf_pp_spectra.F90 udf_pp_SGS.F90 udf_pp_velgrad.F90 udf_pp_hitgen.F90 \
	  udf_pp_force.F90 \
      pp.F90 initialisation.F90 \
      mainloop.F90 test.F90 astr.F90
      
OBJS=$(srs:.F90=.o)

%.o:%.F90
	@mkdir -p $(OBJDIR) 
	$(FC) $(FCFLAGS) $(INCL) $(OPTIONS1) $(OPTIONS2) $(OPTIONS3) $(OPTIONS4) $(OMP) -c -o $(OBJDIR)/$@  $<

default: $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJDIR)/*.o $(LIBS) $(INCL) $(OMP) 

clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET) $(OBJDIR)/*.mod
