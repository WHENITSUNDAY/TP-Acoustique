FC=gfortran
#FFLAGS = -O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero
FFLAGS = -O3 -march=native

EXE=run

all :	$(EXE)

$(EXE) :	mod_param.o mod_solution.o mod_schemas.o ondes.o
	$(FC) $(FFLAGS) -o $(EXE) mod_param.o mod_solution.o mod_schemas.o ondes.o

mod_param.o : mod_param.f90
	$(FC) $(FFLAGS) -c mod_param.f90

mod_solution.o : mod_solution.f90 mod_param.o
	$(FC) $(FFLAGS) -c mod_solution.f90

mod_schemas.o : mod_schemas.f90 mod_param.o mod_solution.o
	$(FC) $(FFLAGS) -c mod_schemas.f90

ondes.o :	ondes.f90 mod_param.o mod_solution.o mod_schemas.o
	$(FC) $(FFLAGS) -c ondes.f90

clean :
	rm -f *.o *.mod $(EXE)