# Makefile for sunset_code
#
#
# ========================== OPTIONS ==============================================================
# -------------------------------------------------------------------------------------------------
# restart    Start from initial conditions (0) or restart file (1)                     (default: 0)
# mpi        Shared only (0) or distributed-shared (1) acceleration                    (default: 1)          
# dim3       Two (0) or three (1) dimensional simulation                               (default: 0)
# vpid       P.I.D control of the velocity (1) or don't (0)                            (default: 0)
# pgrad      Body forces converted to pressure gradient (1) or not (0)                 (default: 0)
# allout     If 3D, output the entire domain (1) or just a slice (0)                   (default: 1)
# ceform     Direct integration (0), log-Cholesky (1), log-conf (2)                    (default: 1)
# newt       Newtonian calculations (1) or not (0)                                     (default: 0)
# fenep      FENE-P (1) or sPTT (0)                                                    (default: 0)
# forder     filter (order) value = 4,6,8,10                                           (default: 8)
# morder     main (order) value = 4,6,8,10                                             (default: 8)
# mcorr      Correct mass conservation (1) or don't (0)                                (default: 1)
# tarout     Compress files as they're written (1) or don't (0)                        (default: 0)
# tracers    Include some Lagrangian tracer particles (1) or don't (0)                 (default: 0)
# kernel     Multiple Kernels (2) or not (1)                                           (default: 1)
# changedim  Routine to enable starting 3d simulation from 2d solution                 (default: 0)  
# -------------------------------------------------------------------------------------------------
#
# EXAMPLE USAGE:
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
ifeq ($(forder),4)
FFLAGS += -Dforder=$(forder)
else
ifeq ($(forder),6)
FFLAGS += -Dforder=$(forder)
else
ifeq ($(forder),10)
FFLAGS += -Dforder=$(forder)
else
FFLAGS += -Dforder=8
endif
endif
endif
ifeq ($(morder),4)
FFLAGS += -Dmorder=$(morder)
else
ifeq ($(morder),6)
FFLAGS += -Dmorder=$(morder)
else
ifeq ($(morder),10)
FFLAGS += -Dmorder=$(morder)
else
FFLAGS += -Dmorder=8
endif
endif
endif

# Restart from dump file.
ifeq ($(restart), 1)
FFLAGS += -Drestart
endif

# MULTIPLE KERNELS
ifeq ($(kernel), 2)
FFLAGS += -Dkernel=$(kernel)
else
FFLAGS += -Dkernel=1
endif

# Non-Newtonian
ifneq ($(newt),1)
# Use something to ensure SPD (Cholesky or log-conf)
ifeq ($(ceform), 0)
FFLAGS += -Ddi
else
ifeq ($(ceform), 2)
FFLAGS += -Dlc
else
FFLAGS += -Dchl
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
ifeq ($(fenep),1)
FFLAGS += -Dfenep
endif

# Tar output files
ifeq ($(tarout),1) 
FFLAGS += -Dtarout
endif

# Three dimensional?
ifeq ($(dim3),1)
FFLAGS += -Ddim3
endif

# Flow velocity PID controlled?
ifeq ($(vpid),1)
FFLAGS += -Dvpid
endif

# Body forces converted to pressure gradient
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

# Tracer particles
ifeq ($(tracers),1)
FFLAGS += -Dtracers
endif

#Change dimension restart?
ifeq ($(changedim),1)
FFLAGS += -Dchangedim
endif

LDFLAGS := -fopenmp -m64 

# Identify directories
SUB_DIRS := common base
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

# identify object files
#parameters come first, as almost everything depends on them.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/common_vars.o
OBJ_FILES += obj/rbfs.o obj/mirror_boundaries.o obj/derivatives.o 
OBJ_FILES += obj/mpi_transfers.o obj/interpolation.o
OBJ_FILES += obj/neighbours.o obj/output.o obj/statistics.o obj/tracer_particles.o
OBJ_FILES += obj/turbulence.o obj/svdlib.o obj/conf_transforms.o
OBJ_FILES += obj/load_data.o obj/setup_domain.o obj/setup_flow.o
OBJ_FILES += obj/labf.o obj/fd.o
OBJ_FILES += obj/characteristic_boundaries.o 
OBJ_FILES += obj/rhs.o obj/rhs_di.o obj/rhs_lc.o
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
	rm -vf ./paraview_files/disc0.vtu	
	rm -vf ./data_out/grilli*
	rm -vf ./data_out/*.tar.gz
	rm -vf ./data_out/tracer*
	rm -vf ./restart/fields*
	rm -vf ./restart/nodes*
	rm -vf ./restart/*.o

