#define the complier
COMPILER = mpif90
# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 

# libraries
LIBS =  
# source list for main program
SOURCES =  helloworld.f90 
		   
test: $(SOURCES)
		${COMPILER} -o helloworld $(FLAGS) $(SOURCES)    
		
clean:
	rm *.o 

clobber:
	rm helloworld
