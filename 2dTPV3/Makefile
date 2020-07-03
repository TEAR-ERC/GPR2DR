# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#-------------------------------------------------------------------------------------------------------
# COMPILER            Intel             GNU                          Used compiler (and linker)
# MODE                Release           Debug,                       Verbosity and Debug level
#                                       Asserts,PeanoProfile,Profile          
# SHAREDMEM           TBB               OMP, None                    Shared-memory parallelisation
# DISTRIBUTEDMEM      MPI               None                         Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                     Type of boundary conditions
# USE_IPO             Off               On, Unsafe                   IPO during compilation (intel only)
# ******************************************************************************************************

# Helper
SHELL = bash
tolower = $(shell echo $(1) | tr '[:upper:]' '[:lower:]')

# Default Paths and name
PEANO_KERNEL_PEANO_PATH=/import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./Peano/peano
PEANO_KERNEL_TARCH_PATH=/import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./Peano/mpibalancing
EXAHYPE_PATH=/import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./ExaHyPE
PROJECT_PATH=/import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./ApplicationExamples/GPRDR/GPRDR_slip4
EXECUTABLE=ExaHyPE-GPRDR

# Default flags from specification file constants
PROJECT_CFLAGS+= -DDim2
ARCHITECTURE=hsw
PROJECT_CFLAGS+= -DALIGNMENT=32

# Additional flags from specification file
PROJECT_CFLAGS+=-DEQNTYPED99 -DODESOLVER
PROJECT_LFLAGS+=-ltecio

# Default sharedmem/distributedmem required by the specification file
ifeq ($(SHAREDMEM),)
    SHAREDMEM=TBB
endif
ifeq ($(DISTRIBUTEDMEM),)
    DISTRIBUTEDMEM=MPI
endif

# If using TBB only, link pthread (comment this if you want to link a specific pthread)
ifeq ($(call tolower,$(SHAREDMEM)),tbb)
    ifeq ($(call tolower,$(DISTRIBUTEDMEM)),none)
        PROJECT_LFLAGS+= -pthread
    endif
endif

# Default options required by the specification file
PROJECT_CFLAGS += -DOPT_KERNELS
ifneq ($(call tolower,$(MODE)),release)
ifneq ($(call tolower,$(MODE)),)
    PROJECT_CFLAGS += -DTEST_OPT_KERNEL
endif
endif



# Main Makefile include
-include /import/deadlock-data/dli/ExaHyPE/ExaHyPE-Engine/./ExaHyPE/Makefile


all: 
	@echo GPRDR
	@echo =================
	@echo An ExaHyPE solver