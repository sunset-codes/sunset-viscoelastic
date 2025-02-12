# Makefile for sunset_code
#
#
# ========================== OPTIONS ==============================================================
# -------------------------------------------------------------------------------------------------
# restart    Start from initial conditions (0) or restart file (1)                     (default: 0)
# mpi        Shared only (0) or distributed-shared (1) acceleration                    (default: 1)          
# dim3       Two (0) or three (1) dimensional simulation                               (default: 0)
# pgrad      Drive the flow with a pressure gradient and P.I.D control                 (default: 0)
# allout     If 3D, output the entire domain (1) or just a slice (0)                   (default: 1)
# ceform     Direct integration (0), Cholesky (1), log-cholesky (2), log-conf (3)      (default: 3)
# newt       Newtonian calculations (1) or not (0)                                     (default: 0)
# fenep      FENE-P (1) or sPTT (0)                                                    (default: 1)
# morder     m (order) value = 4,6,8,10                                                (default: 8)
# mcorr      Correct mass conservation (1) or don't (0)                                (default: 1)
# -------------------------------------------------------------------------------------------------
#
# EXAMPLE USAGE:
# make dim3=0 mpi=1 pgrad=0 spd=0
#

#
# Choose compiler depending on whether mpi
ifneq ($(mpi),0)
FC := mpifort
LD := mpifort
else
FC := gfortran
LD := gfortran
endif

# Set compiler flags based on make options
CFLAGS := -Wall -O3 -g -m64
FFLAGS := -fopenmp -fbounds-check -ffpe-trap=zero -O3 -Wall -g -J./obj -I./obj -m64

# Order of the numerical scheme (even, from 4 to 10, default 8)
ifeq ($(morder),4)
FFLAGS += -Dorder=$(morder)
else
ifeq ($(morder),6)
FFLAGS += -Dorder=$(morder)
else
ifeq ($(morder),10)
FFLAGS += -Dorder=$(morder)
else
FFLAGS += -Dorder=8
endif
endif
endif

# Restart from dump file.
ifeq ($(restart), 1)
FFLAGS += -Drestart
endif

# Non-Newtonian
ifneq ($(newt),1)
# Use something to ensure SPD (Cholesky or log-conf)
ifeq ($(ceform), 0)
FFLAGS += -Ddi
else
ifeq ($(ceform), 1)
FFLAGS += -Dch
else
ifeq ($(ceform), 2)
FFLAGS += -Dchl
else
FFLAGS += -Dlc
endif
endif
endif

# or Newtonian
else
FFLAGS += -Dnewt
endif

# Multiprocessor? (use mpi?)
ifneq ($(mpi),0)
FFLAGS += -Dmp
endif

# FENE-P?
ifneq ($(fenep),0)
FFLAGS += -Dfenep
endif

# Three dimensional?
ifeq ($(dim3),1)
FFLAGS += -Ddim3
endif

# Flow driven by pressure gradient?
ifeq ($(pgrad),1)
FFLAGS += -Dpgrad
endif

# Correction for mass conservation
ifneq ($(mcorr),0)
FFLAGS += -Dmcorr
endif

# Output entire domain?
ifneq ($(allout),0)
FFLAGS += -Dallout
endif

LDFLAGS := -fopenmp -m64 

# Identify directories
SUB_DIRS := common base
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

# identify object files
#parameters come first, as almost everything depends on them.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/common_vars.o
OBJ_FILES += obj/rbfs.o obj/mirror_boundaries.o obj/derivatives.o 
OBJ_FILES += obj/mpi_transfers.o
OBJ_FILES += obj/neighbours.o obj/output.o obj/statistics.o 
OBJ_FILES += obj/turbulence.o obj/svdlib.o obj/mat2lib.o
OBJ_FILES += obj/load_data.o obj/setup_domain.o obj/setup_flow.o
OBJ_FILES += obj/labf.o obj/fd.o
OBJ_FILES += obj/characteristic_boundaries.o obj/rhs.o
OBJ_FILES += obj/step.o
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.F90,obj/%.o,$(wildcard $(sdir)/*.F90)))

vpath %.F90 $(SRC_DIR)

#-------
default: sunset
sunset: $(OBJ_FILES)
	$(LD) -o $@ $^ $(LDFLAGS)

obj/%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -vf ./obj/*.o
	rm -vf ./obj/*.mod
	rm -vf ./sunset
	rm -rfv fort.*	
	rm -vf ./data_out/fields*
	rm -vf ./data_out/nodes*
	rm -vf ./data_out/flame*
	rm -vf ./data_out/time.out
	rm -vf ./data_out/statistics/*.out
	rm -vf ./paraview_files/LAYER*
	rm -vf ./data_out/grilli*

