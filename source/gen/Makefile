#==================================================
# 
# Dr. Xiaohu Guo
# xiaohu.guo@stfc.ac.uk
#
#==================================================

SHELL=/bin/sh
MAKE=gmake
FC=gfortran

FCFLAGS=-W -O3 -m64 -g -fbounds-check
FCA = -m64

VPATH =../common  

OBJS= kind_parameters.o common_parameter.o global_variables.o datclass.o

default: gen2D

gen2D: $(OBJS)
	$(FC) $(FCA) -o $@ $(OBJS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod gen2D fort*
