OBJ=QM.f90 chebft.o 
F77=f95

exe: $(OBJ)
	$(F77) -o exe -O2 $(OBJ) $(LIBS)

chebft.o: chebft.f90
	$(F77)  -c -O2 chebft.f90
