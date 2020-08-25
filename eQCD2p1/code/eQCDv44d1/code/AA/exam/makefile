# Start of the makefile
# Defining variables
objects=chebft.o derivs.o dtVdiff1.o dtVdiff2.o eta.o expaPoint.o FRG.o fun.o gapEq.o GAULEG.o initial.o intLin.o IRenhan.o newt.o nfdx_com.o nfdx.o odeint.o phypoint.o PolyakovEq.o eQCD.o selfEQ.o thresholdfunc.o vInf.o
# Makefile
exe: $(objects)
	gfortran -o exe $(objects) 
%.o: %.f90
	gfortran -c -g $< 
# Cleaning everything
clean:
	rm $(objects)
# End of the makefile
