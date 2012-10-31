#define the complier
COMPILER = mpif90
# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 -mp

# libraries
LIBS = -L${FFTW_LINK} -lfftw3 
# source list for main program
SOURCES =  KgSemiImp3d.f90 initialdata.f90 savedata.f90 getgrid.f90 \
			storeold.f90 saveresults.f90 enercalc.f90 forwardfft.f90 \
			backwardfft.f90
		   
test: $(SOURCES)
		${COMPILER} -o kg $(FLAGS) $(SOURCES) $(LIBS)   
		
clean:
	rm *.o 
