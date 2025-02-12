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

VPATH = .  

OBJS := main.o

default: vtkbuild

vtkbuild: $(OBJS)
	$(FC) $(FCA) -o $@ $(OBJS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod vtkbuild
