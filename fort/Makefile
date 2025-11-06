
FC = ifort
FFLAGS = -O2
LDFLAGS = -qmkl

all: numpoten cLLSM


numpoten: numpoten.f
	$(FC) $(FFLAGS) numpoten.f -o numpoten
	
cLLSM: cLLSM.f90
	$(FC) $(LDFLAGS) $(FFLAGS) cLLSM.f90 -o cLLSM

.PHONY: clean
clean:
	rm -f numpoten cLLSM constants.mod

