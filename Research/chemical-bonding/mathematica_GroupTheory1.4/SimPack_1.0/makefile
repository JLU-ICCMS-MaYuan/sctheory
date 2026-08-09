#****s* make
#
# NAME
#  make
#
# AUTHOR
#  W. Hergert
# 
# TYPE
#  makefile   
# 
# USAGE
#  Makefile for FORTRAN programs in SimPack
#
#  programs for : fitting, bandstructure calculation, DOS calculation  
#
#  Note, the programs correspond to different versions of FORTRAN. Thus,
#  one has to be careful with the compiler.
#
#  gfortran : use the argument std
#  ifort    : no problems
#  Use the appropriate absolute or relatice path names.
#
#  For other compilers you have to check by yourself.
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
F90C	= gfortran -std='legacy'
LFLAGS  = -L/usr/local/lib   -llapack     -lblas 
#
# ifort 
#
#F90C	= ifort 
#FFLAGS  = -I"${MKLROOT}/include"
#LFLAGS  =   -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl
#
#---------------------------------
#
# define one of the  Hamiltonians and the path for executables
#
# example
#
HAMIL  = fcc_s.f90
DIREC  = example
NAME   = s
#
# simple s-band
#
#HAMIL   = htbs.f90
#DIREC  = applications/sb/exe
#NAME    = s
#
# fcc-Hamiltonian d-bands only
#
#HAMIL   = fcc_d.f90
#DIREC   = applications/Cu_d/exe
#NAME    = d
#
# fcc-Hamitonian with basis of -s, p-, and d-electrons
#
#HAMIL   = fcc_spd.f90
#DIREC   =applications/Cu_spd/exe
#NAME    = spd
#
# path variables
#
SP_HOME     = $(PWD)
SRC         = $(SP_HOME)/src
#
# help menue
#
help:
	clear
	@echo
	@echo  "----------------------------------------------------------------"
	@echo  " SimPack - a FORTRAN supplement to GTPack                      |"
	@echo  "---------------------------------------------------------------|"
	@echo  " * uses TB-Hamiltonians exported to FORTRAN from GTPack        |"
	@echo  " * fits TB_Hamiltonians against ab initio data                 |"
	@echo  " * calculates bandstructures and densities of states           |"
	@echo  "---------------------------------------------------------------|"
	@echo  " Change the makefile to include :                              |"
	@echo  " * your favourite compiler                                     |"
	@echo  " * the correct Hamiltonian                                     |"
	@echo  "---------------------------------------------------------------|"
	@echo  " make #                                                        |"
	@echo  "# = tb_fit   - Fit (least squares, genetic algorithm, simplex) |"
	@echo  "# = tb_bands - band structures                                 |"
	@echo  "# = tb_dos   - densities of states                             |"
	@echo  "#                                                              |"                                              
	@echo  "# = clean    - delete intermediate files (compilation)         |"
	@echo  "# = help     - this overview                                   |"
	@echo  "----------------------------------------------------------------"
#
.f90.o:
	$(F90C) -c $(FFLAGS) $<
	#        
# clean up
#
clean:
	rm *.o *.mod *.bak               

#
# tb_fit 
#
tb_fit:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/general/typedef1.f90 \
	$(SRC)/math/asa047.f90  $(SRC)/math/dnlq1.f90 $(SRC)/math/pikaia.f \
	$(SRC)/general/tb_data_1.f90 $(SRC)/fit/tb_abinitio.f90 \
	$(SRC)/hamiltonians/$(HAMIL) \
	$(SRC)/fit/tb_ls.f90 $(SRC)/fit/tb_simplex.f90 \
	$(SRC)/fit/tb_ga.f90  $(SRC)/fit/parmfix.f90 \
	$(SRC)/fit/tb_fit.f90 \
 	-o $(DIREC)/$@_$(NAME).x  
#
# energy bands
#
tb_bands:
	$(F90C) $(FFLAGS) \
	$(SRC)/general/typedef1.f90   $(SRC)/general/tb_data.f90 \
	$(SRC)/hamiltonians/$(HAMIL) $(SRC)/kspace/kpoints_3d.f90 \
        $(SRC)/bands/tb_bands.f90  $(LFLAGS) \
	-o $(DIREC)/$@_$(NAME).x  
#
# DOS and partial DOS
#
tb_dos:
	$(F90C) $(FFLAGS) $(LFLAGS) $(SRC)/general/typedef.f90 \
	$(SRC)/kspace/kpoints_3d.f90 $(SRC)/math/dqsf.f90 \
	$(SRC)/dos/dos_3d.f90 $(SRC)/dos/glbi.f90 \
	 -o $(DIREC)/$@_$(NAME).x  	
#
#***

