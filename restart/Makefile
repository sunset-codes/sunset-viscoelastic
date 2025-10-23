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

# Three dimensional?
ifeq ($(dim3),1)
FCFLAGS += -Ddim3
endif

#Start 3D simulation with 2D solution?
ifeq ($(changedim),1)
FCFLAGS += -Dchangedim
endif

VPATH = .  

OBJS := redecomp.o

default: redecomp

redecomp: $(OBJS)
	$(FC) $(FCA) -o $@ $(OBJS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod redecomp
