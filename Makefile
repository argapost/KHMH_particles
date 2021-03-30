#=======================================================================
# Makefile
#=======================================================================

CMP = gcc #intel, gcc
FC = mpif90

######## NETCDF include ########
NETCDFloc = /usr/local
NETCDFlib = -I${NETCDFloc}/include -L${NETCDFloc}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -lm
#NETCDFloc = 
#NETCDFlib = -lnetcdf -lnetcdff

MODDIR = ./mod
SRCDIR = ./src

### List of all files for the main code
SRC = $(SRCDIR)/interpolate_mod.f90 $(SRCDIR)/particleTracer.f90 $(SRCDIR)/initialize.f90 $(SRCDIR)/load_velocity.f90 $(SRCDIR)/save.f90 
OBJ = $(SRC:%.f90=%.o)

###### OPTIONS settins ########
OPT = -I$(SRCDIR) $(NETCDFlib)
LINKOPT = $(NETCDFlib)


# -------------------------------------------------

all: particleTracer

particleTracer : $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJ) $(NETCDFlib)

$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(OPT) -c $< 
	mv $(@F) ${SRCDIR}
	# mv *.mod ${SRCDIR}

%.o : %.f90
	$(FC) -c $<


.PHONY: clean

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod particleTracer
