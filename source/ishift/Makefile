#==================================================
# 
# Dr J. King
#
#==================================================

SHELL=/bin/sh
MAKE=gmake
FC=gfortran

FCFLAGS=-W -O3 -m64 -g -fbounds-check -fopenmp
FCA = -m64 -fopenmp

VPATH =../common  

OBJS := kind_parameters.o common_parameter.o common_vars.o 
OBJS += neighbours.o boundaries.o 
OBJS += inputoutput.o 
OBJS += main.o

default: ishift

ishift: $(OBJS)
	$(FC) $(FCA) -o $@ $(OBJS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod ishift fort*
