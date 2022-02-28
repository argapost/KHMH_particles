#=======================================================================
# Makefile
#=======================================================================

ifeq ($(CMP),gcc)
FC = mpif90
FFLAGS =
NETCDFloc = /usr/local
NETCDFlib = -I${NETCDFloc}/include -L${NETCDFloc}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -lm
else
FC = mpiifort
FFLAGS = -fpp -O3 -xHost -ipo -fp-model fast=2 -mcmodel=large -qopenmp
NETCDFloc =
NETCDFlib = -lnetcdf -lnetcdff
endif

MODDIR = ./mod
SRCDIR = ./src

### List of all files for the main code
SRC = $(SRCDIR)/KHMH_particles.f90 $(SRCDIR)/load_timestep.f90 $(SRCDIR)/save.f90
OBJ = $(SRC:%.f90=%.o)

###### OPTIONS settins ########
OPT = -I$(SRCDIR) $(NETCDFlib) $(FFLAGS)
LINKOPT = $(NETCDFlib) $(FFLAGS)


# -------------------------------------------------

all: KHMH_particles

KHMH_particles : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(LINKOPT) $(OBJ) $(NETCDFlib)

$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) -c $<
	mv $(@F) ${SRCDIR}
	# mv *.mod ${SRCDIR}

%.o : %.f90
	$(FC) $(FFLAGS) -c $<


.PHONY: clean

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod KHMH_particles
