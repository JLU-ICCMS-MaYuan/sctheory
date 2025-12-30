#****s* maketest
# NAME
#  makefile
#
# AUTHOR
#  W. Hergert 
#
# TYPE
#  makefile   
# 
# USAGE
#  Makefile for tests of SImplex method according to Nelder-Mead
#
#  Note, the programs correspond to the user's guide to Pikaia 1.0
#
#  Compilers tested:
#
#  gfortran : use the argument std
#  ifort    : no problems
#
#  For other compilers you  have to check by yourself.
#
# OPTIONS
#  -
#
#**********************************************************************
#
# compiler (uncomment what you need)
#
# gfortran
# 
F90C	= /usr/local/bin/gfortran -std='legacy'
LFLAGS  = -L/usr/local/lib    -llapack     -lblas
#
# ifort 
#
#F90C	= /opt/intel/ifort
#FFLAGS  = -g 
#LFLAGS  = -mkl
#
#---------------------------------
#
#
SP_HOME     = $(PWD)
SRC         = $(SP_HOME)/src
EXE         = $(SP_HOME)/exe
#
help:
	clear
	@echo
	@echo  "--------------------------------------------------------------"
	@echo  " Test programs for Nelder-Mead simplex method test           |"
	@echo  "-------------------------------------------------------------|"
	@echo  " make #                                                      |"
	@echo  "                                                             |"                                              
	@echo  "# = orig     - original as047 test program                   |"
	@echo  "# = mod      - modfied test program                          |"
	@echo  "# = clean    - delete intermediate files (compilation)       |"
	@echo  "# = help     - this overview                                 |"
	@echo  "--------------------------------------------------------------"     
#
.f90.o:
	$(F90C) -c $(FFLAGS) $<
	#        
# clean up
#
clean:
	rm *.o *.mod *.bak               

#
# Nelder-Mead original test
#
orig:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/orig/asa047.f90  $(SRC)/orig/asa047_test.f90  \
	-o $(EXE)/nm_$@.x  
#
# Nelder-Mead modified test program
#
mod:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/typedef1.f90 $(SRC)/timestamp.f90 $(SRC)/asa047.f90 \
	$(SRC)/nm_helical.f90 $(SRC)/nm_powell.f90 $(SRC)/nm_quartic.f90 \
	$(SRC)/nm_rosenbrock.f90 $(SRC)/nm_test.f90 \
	-o $(EXE)/nm_$@.x 
	
#***
    
