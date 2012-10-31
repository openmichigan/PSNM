
	# All settings here for use on Benson Muite's laptop

	#define the complier
	COMPILER = gfortran
	# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 -fopenmp
	# libraries
	LIBS =  -L/usr/local/lib -lfftw3 -lm 		
	



# source list for main program
SOURCES =  NavierStokes3DfftIMR.f90 savedata.f90

test: $(SOURCES)
		${COMPILER} -o NSIncompress $(FLAGS) $(SOURCES) $(LIBS)

clean:
	rm *.o
