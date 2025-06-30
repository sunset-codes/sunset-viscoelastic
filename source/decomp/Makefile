#==================================================
# 
# Dr J. King
#
#==================================================

SHELL=/bin/sh
MAKE=gmake
FC=gfortran

FCFLAGS=-W -O3 -m64 -g -fbounds-check -fopenmp
ifeq ($(shifting),1)
FCFLAGS += -Dshft
endif

FCA = -m64 -fopenmp

VPATH =../common  

OBJS := kind_parameters.o common_parameter.o common_vars.o 
OBJS += neighbours.o boundaries.o 
OBJS += inputoutput.o 
OBJS += main.o

ifeq ($(shifting),1)
default: shifting

shifting: $(OBJS)
	$(FC) $(FCA) -o $@ $(OBJS)
	rm *.o

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod decomp shifting fort*
	
else
default: decomp

decomp: $(OBJS)
	$(FC) $(FCA) -o $@ $(OBJS)
	rm *.o

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod decomp shifting fort*
endif	
	

