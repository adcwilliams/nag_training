FC=nagfor

fd1d_heat_explicit.exe: libfmw.a fd1d_heat_explicit.o
	$(FC) -o $@ fd1d_heat_explicit.o -L. -lfmw

Types_mod.o : Types_mod.f90
	$(FC) -I. -c $<
CFL_mod.o : CFL_mod.f90 
	$(FC) -I. -c $<
fd1d_heat_explicit.o : fd1d_heat_explicit.f90 
	$(FC) -c $<
IO_mod.o : IO_mod.f90
	$(FC) -c $<	
RHS_mod.o : RHS_mod.f90
	$(FC) -I. -c $<	
Solver_mod.o : Solver_mod.f90
	$(FC) -c $<
Types_mod.o : Types_mod.f90 
	$(FC) -c $<
libfmw.a :  Types_mod.o CFL_mod.o IO_mod.o RHS_mod.o Solver_mod.o#
	ar rcs $@ $^

.PHONY: clean

clean:
	rm -f *.mod *.o *.png *.exe *.a
