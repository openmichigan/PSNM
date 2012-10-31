#define the complier
COMPILER = mpif90
# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 

# libraries
LIBS =  
# source list for main program
SOURCES =  montecarloserial.f90 
		   
test: $(SOURCES)
		${COMPILER} -o montecarloserial $(FLAGS) $(SOURCES)    
		
clean:
	rm *.o 

clobber:
	rm  montecarloserial
