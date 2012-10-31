#define the complier
COMPILER = mpif90
# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 

# libraries
LIBS = -L${FFTW_LINK} -lfftw3 -lm 
# source list for main program
SOURCES =  heat.f90 
		   
test: $(SOURCES)
		${COMPILER} -o heat $(FLAGS) $(SOURCES) $(LIBS)   
		
clean:
	rm *.o 
