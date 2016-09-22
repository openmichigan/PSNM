# Set options to allow for conditional compilation, use only one of the options below
# These are processed using the C pre-processor
#
#
OPTIONS= 
CPP=-cpp
	# All settings here for use Muite Laptop

	COMPILER =  mpif90  
	decompdir= ./2decomp_fft
	FFTW_LINK=./fftwinstall
	# compilation settings, optimization, precision, parallelization
	FLAGS =  -O3  
	LIBS = -I${FFTW_LINK}/include -L${FFTW_LINK}/lib -lfftw3 		
	
	DECOMPLIB = -I${decompdir}/include -L${decompdir}/lib -l2decomp_fft 


# libraries
# source list for main program

SOURCES = brusselator.f90 initialdata.f90 savedata.f90 getgrid.f90 \
		 linear.f90 nonlinear1.f90 nonlinear2.f90 

OBJECTS = brusselator.o initialdata.o savedata.o getgrid.o \
		 linear.o nonlinear1.o nonlinear2.o 


Brusselator: $(OBJECTS)
		${COMPILER} -o Brusselator  $(OPTIONS) $(FLAGS) $(CPP) $(OBJECTS)  $(DECOMPLIB) $(LIBS)

$(OBJECTS): $(SOURCES)
		${COMPILER} -c $(OPTIONS) $(CPP) $(DECOMPLIB) $(LIBS) $(SOURCES)

clean:
	rm -f *.o
clobber:	
	rm -f Brusselator
