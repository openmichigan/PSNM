#define the complier
COMPILER = gfortran
# compilation settings, optimization, precision, parallelization
FLAGS = -O0 -fopenmp 


# libraries
LIBS = -L/usr/local/lib -lfftw3 -lm 
# source list for main program
SOURCES =  NLSsplitting.f90 
		   
test: $(SOURCES)
		${COMPILER} -o NLSsplitting $(FLAGS) $(SOURCES) $(LIBS)   
		
clean:
	rm *.o 

clobber:
	rm NLSsplitting
