#============================================================================#
#                                  GFORTRAN                                  #
#============================================================================#

# 1. Compiler Configuration
# use the gfortran version >= 10
# and put the complier you want to use here
FC = gfortran

# 2. Base Options
# Standard scientific flags
FFLAGS = -g -O3 -fno-second-underscore -ftree-vectorize -falign-loops=16 -ffast-math -fassociative-math -fallow-argument-mismatch
# Debug flags (commented out)
# FFLAGS = -pedantic-errors -Wall -Wconversion -Wextra -Wunreachable-code -fbacktrace -g3 -fbounds-check -O0 -fstack-protector-all -fno-automatic -Wuninitialized -ftrapv -fno-automatic -fimplicit-none
# Options and Path
#  FFLAGS = -g -O3 -fno-second-underscore -ftree-vectorize -falign-loops=16 -ffast-math -fassociative-math -fallow-argument-mismatch
#  FFLAGS = -pedantic-errors -Wall -Wconversion -Wextra -Wunreachable-code -fbacktrace -g3 -fbounds-check -O0 -fstack-protector-all -fno-automatic -Wuninitialized -ftrapv -fno-automatic -fimplicit-none

# 3. Platform Specific Checks (macOS Fixes)
# Auto-detect macOS SDK path and set Deployment Target
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # Fix for newer macOS (Sequoia/Sonoma) clang error "invalid version"
    # This export only affects the shell running these commands, not your global system
    export MACOSX_DEPLOYMENT_TARGET=14.0

    # Find the SDK path if not already set
    ifndef SDKROOT
        SDKROOT := $(shell xcrun --show-sdk-path)
    endif
    
    # Append the sysroot flag to existing FFLAGS (using += preserves the flags above)
    FFLAGS += -isysroot $(SDKROOT)
endif


.MAKEOPTS: -k -s

.SUFFIXES: .f90

# 4. Object Files
OBJ = global_variables.o \
      opkda1.o \
      opkda2.o \
   	  opkdmain.o \
      auxilary.o \
      shielding.o \
      init_chem.o \
      init_struc.o \
      chemistry.o \
      ode_solve.o \
	  disk_output.o
		
PROGS = disk_main disk_extract
# point to the correct dir for the bin file
PROG_PATH = ../bin/

# 5. Build Targets
ALL: $(PROGS) install
	echo "Compilation finished"
	
install: $(PROGS)
	mv $(PROGS) $(PROG_PATH)
	
global_variables.o: global_variables.f90
	$(FC) -c $(FFLAGS) global_variables.f90

%.o : %.f90 global_variables.o
	$(FC) -c $(FFLAGS) $<

disk_main: $(OBJ) disk_main.o
	$(FC) $(FFLAGS) disk_main.o $(OBJ) -o disk_main
		
disk_extract: $(OBJ) disk_extract.o
	$(FC) $(FFLAGS) disk_extract.o $(OBJ) -o disk_extract
			
clean:
	rm -f -r $(OBJ) *.dSYM disk_main disk_extract
