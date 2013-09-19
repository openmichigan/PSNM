	# All settings here for use on FLUX, a cluster at the University of Michigan
	# Center for Advanced Computing (CAC), using INTEL nehalem hardware, 
	# Need to load fftw module

	COMPILER =  mpif90  
	decompdir=/home/msrinivasan/develop/others/2decomp_fft
	# compilation settings, optimization, precision, parallelization
	FLAGS =  -O0  
	#LIBS = -L${FFTW_LINK} 		
	
	DECOMPLIB = -I${decompdir}/include -L${decompdir}/lib -l2decomp_fft


# libraries
# source list for main program
SOURCES = KgSemiImp3d.f90 initialdata.f90 savedata.f90 getgrid.f90 \
			storeold.f90 saveresults.f90 enercalc.f90 readinputfile.f90

Kg: $(SOURCES)
		${COMPILER} -o Kg $(FLAGS) $(SOURCES)  $(DECOMPLIB)


clean:
	rm -f *.o
clobber:	
	rm -f Kg
