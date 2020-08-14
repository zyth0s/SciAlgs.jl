
# gfortran -shared -O2 -fPIC fortrancode.f90 -o fortranlib.so 
#
FC=gfortran

FFLAGS=-Wall -shared -O2 -fPIC

SOURCES=fortrancode.f90

lib: $(SOURCES)
	$(FC) $(FFLAGS) $(SOURCES) -o fortranlib.so

clean:
	rm *.mod *.so
