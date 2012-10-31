COMPILER =  mpif90 
decompdir= ../2decomp_fft
FLAGS = -O0 
	
DECOMPLIB = -I${decompdir}/include -L${decompdir}/lib -l2decomp_fft
LIBS = #-L${FFTW_LINK}  -lfftw3 -lm 
SOURCES = NavierStokes3DfftIMR.f90 savedata.f90  

ns3d: $(SOURCES)
		${COMPILER} -o ns3d $(FLAGS) $(SOURCES) $(LIBS) $(DECOMPLIB) 

clean:
	rm -f *.o
	rm -f *.mod
clobber:	
	rm -f ns3d

