
# gfortran -shared -O2 -fPIC fortran_iso_c.f90 -o fortranlib.so 
#
FC=gfortran

FFLAGS=-Wall -shared -O2 -fPIC

SOURCES=fortran_iso_c.f90

lib: $(SOURCES)
	$(FC) $(FFLAGS) $(SOURCES) -o fortranlib.so

clean:
	rm *.mod *.so
