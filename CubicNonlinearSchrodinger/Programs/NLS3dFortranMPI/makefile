
COMPILER =  mpif90 
decompdir=../2decomp_fft
FLAGS = -O0 
	
DECOMPLIB = -I${decompdir}/include -L${decompdir}/lib -l2decomp_fft
LIBS = -L${FFTW_LINK} -lfftw3_threads -lfftw3 -lm 
SOURCES = NLSsplitting.f90  

NLSsplitting: $(SOURCES)
		${COMPILER} -o NLSsplitting $(FLAGS) $(SOURCES) $(LIBS) $(DECOMPLIB) 


clean:
	rm -f *.o
	rm -f *.mod
clobber:	
	rm -f NLSsplitting

