#define the complier
COMPILER = mpif90
# compilation settings, optimization, precision, parallelization
	FLAGS = -O0 

# libraries
LIBS =  
# source list for main program
SOURCES =  montecarloparallel.f90 
		   
test: $(SOURCES)
		${COMPILER} -o montecarloparallel $(FLAGS) $(SOURCES)    
		
clean:
	rm *.o 

clobber:
	rm  montecarloparallel
