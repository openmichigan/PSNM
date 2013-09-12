# Makefile for Beacon at NICS 
# Assumes: module load xeon/fftw/3.2.1
# Assumes: module load PE-intel
# Assumes: export LD_LIBRARY_PATH=/nics/b/home/mvanmoer/apps/Beacon/offload/Mesa/9.1/lib:$LD_LIBRARY_PATH

FC = mpiifort 
FCWARNINGS =
FCDEBUG = 
FCFLAGS = -O0 $(FCWARNINGS) $(FCDEBUG) 

CXX = mpiicc
CXXWARNINGS =
CXXDEBUG =
CXXFLAGS = -c $(CXXWARNINGS) $(CXXDEBUG)

DECOMPDIR = /nics/b/home/mvanmoer/apps/Beacon/offload/2decomp_fft/1.5.847

PVINSTALLDIR = /nics/b/home/mvanmoer/apps/ParaView/4.0.1
PVINCDIR = -I$(PVINSTALLDIR)/include/paraview-4.0
PVLIBS = -lvtkPVCatalystPython26D-pv4.0 -lvtkPVCatalyst-pv4.0 -lvtkCommonCore-pv4.0 -lvtkCommonDataModel-pv4.0 -lvtkCommonDataModel-pv4.0 -lvtkPVPythonCatalyst-pv4.0

INCDIRS = $(PVINCDIR) -I$(DECOMPDIR)/include
LIBDIRS = -L$(PVINSTALLDIR)/lib/paraview-4.0 -L$(DECOMPDIR)/lib -L${HOME}/apps/Beacon/offload/Mesa/9.1/lib

LIBS =  $(FFTW3_LIB) $(PVLIBS) -lOSMesa32 -l2decomp_fft

BIN = Kg
all: $(BIN)

VTKCellBasedDataSet.o: VTKCellBasedDataSet.cxx
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

KGadaptor.o: KGadaptor.f90 VTKCellBasedDataSet.o
	$(FC) -c $(FCFLAGS) $(INCDIRS) -o $@ $^

%.o: %.f90
	$(FC) -c $< -o $@ $(LIBDIRS) $(INCDIRS) $(LIBS)

OBJECTS = VTKCellBasedDataSet.o KGadaptor.o enercalc.o \
	initialdata.o savedata.o storeold.o getgrid.o readinputfile.o saveresults.o KgSemiImp3d.f90

$(BIN): $(OBJECTS)
	$(FC) -o $@ $(FCFLAGS) $(INCDIRS) $(LIBDIRS) -Wl,-rpath $(PVINSTALLDIR)/lib/paraview-4.0 $(OBJECTS) $(LIBS)

clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) $(BIN)
