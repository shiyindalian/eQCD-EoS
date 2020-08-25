# Start of the makefile
# Defining variables
objects=derivs.o dtVdiff1.o eta.o FRG.o fun.o gapEq.o GAULEG.o initial.o intLin.o IRenhan.o lblf.f90 mass.o newt.o nfdx_com.o nfdx.o odeint.o phypoint.o PolyakovEq.o eQCD.o selfEQ.o thresholdfunc.o
# Makefile
eqcd: $(objects)
	gfortran -o eqcd -O2 $(objects) 
%.o: %.f90
	gfortran -c -O2 -g $< 
