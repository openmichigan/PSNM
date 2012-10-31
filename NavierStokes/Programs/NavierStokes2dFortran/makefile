ifeq ($(HOST),flux)
	# All settings here for use on FLUX, a mixed cluster at the University of Michigan
	# Center for Advanced Computing (CAC), primarily using INTEL nehalem I& hardware, 
	# Need to load fftw module

	#define the complier
	COMPILER = ifort
	# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 -openmp
	# libraries
	LIBS =  -L${FFTW_LINK} -lfftw3 -lm 		
	
endif

ifeq ($(HOST),muitelaptop)
	# All settings here for use on Benson Muite's laptop

	#define the complier
	COMPILER = gfortran
	# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 -fopenmp
	# libraries
	LIBS =  -L/usr/local/lib -lfftw3 -lm 		
	
endif



# source list for main program
SOURCES =  navierstokes.f90

test: $(SOURCES)
		${COMPILER} -o NSIncompress $(FLAGS) $(SOURCES) $(LIBS)

clean:
	rm *.o
