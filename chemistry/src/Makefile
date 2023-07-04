#============================================================================#
#                                  GFORTRAN                                  #
#============================================================================#
# use the gfortran version >= 10. and put the complier you want to use here
  FC = gfortran

# Options and Path
  FFLAGS = -g -O3 -fno-second-underscore -ftree-vectorize -falign-loops=16 -ffast-math -fassociative-math -fallow-argument-mismatch
#  FFLAGS = -pedantic-errors -Wall -Wconversion -Wextra -Wunreachable-code -fbacktrace -g3 -fbounds-check -O0 -fstack-protector-all -fno-automatic -Wuninitialized -ftrapv -fno-automatic -fimplicit-none

.MAKEOPTS: -k -s

.SUFFIXES: .f90


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
