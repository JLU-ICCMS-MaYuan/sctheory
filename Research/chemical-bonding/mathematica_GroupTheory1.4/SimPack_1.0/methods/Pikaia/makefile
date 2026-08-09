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
#  Makefile for test of Pikaia co
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
F90C	= gfortran -std='legacy'
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
	@echo  " Test programs for Pikaia see user's guide to Pikaia 1.0     |"
	@echo  "-------------------------------------------------------------|"
	@echo  " make                                                        |"
	@echo  " # = test0  - maximizing function  (sec. 5.1 p 35.ff)        |"
	@echo  " # = test1a - 2 variables lin. fit (sec. 5.2 p 41.ff)        |"
	@echo  " # = test1b - nonlinear fit        (sec. 5.3 p 47.ff)        |"
	@echo  " # = test2  - gen. LSQR            (sec. 5.4 p 56.ff)        |"
	@echo  " # = test3  - robust estimators    (sec. 5.5 p 59.ff)        |"
	@echo  "                                                             |"
	@echo  " # = orig1b - original nonlinear fit                         |"                
	@echo  "                                                             |"                                              
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
# fit0 
#
test0:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/xpk0.f90  $(SRC)/twod.f90 $(SRC)/pikaia.f $(SRC)/ga_input.f90 \
	-o $(EXE)/$@.x  
#
# fit1a
#
test1a:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/xpk1a.f90  $(SRC)/fit1a.f90 $(SRC)/pikaia.f  $(SRC)/ga_input.f90 \
	-o $(EXE)/$@.x  
#
# fit1b
#
test1b:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/xpk1b.f90  $(SRC)/fit1b.f90 $(SRC)/pikaia.f  $(SRC)/ga_input.f90 \
	-o $(EXE)/$@.x  
#
# original fit1b (modified input and rescaling)
#
orig1b:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/pikaia.f  $(SRC)/orig/xpk1bm.f \
	$(SRC)/orig/fit1b.f   $(SRC)/ga_input.f90 \
	-o $(EXE)/$@.x  
#
# fit2
#
test2:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/xpk2.f90  $(SRC)/fit2.f90 $(SRC)/pikaia.f $(SRC)/ga_input.f90  \
	-o $(EXE)/$@.x  
#
# fit3
#
test3:
	$(F90C) $(FFLAGS) $(LFLAGS) \
	$(SRC)/xpk3m.f90  $(SRC)/fit3.f90 $(SRC)/pikaia.f  \
	-o $(EXE)/$@.x  			
 			
#***
    
