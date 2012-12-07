# Nautilus 
# Assumes: PE-intel, since ParaView was built with Intel.
# Assumes: module load paraview/3.14.1
CXX = icpc
FC = ifort

DEBUG = -g
WARNINGS =
FCFLAGS = $(WARNINGS) $(DEBUG)
CXXFLAGS = $(DEBUG)

DECOMPDIR = /nics/c/home/bmuite/2decomp_fft
FFTW_PATH = /sw/analysis/fftw/3.3/sles11.1_intel11.1
PVDIR = /sw/analysis/paraview/3.14.1/sles11.1_intel11.1.038
PVSRCDIR = $(PVDIR)/ParaView-3.14.1-Source
PVBUILDDIR = $(PVDIR)/build
PVINCDIRS = -I$(PVSRCDIR)/CoProcessing/Adaptors/FortranAdaptors -I$(PVBUILDDIR)/CoProcessing/Adaptors/FortranAdaptors -I$(PVSRCDIR)/VTK/Filtering -I$(PVSRCDIR)/VTK/Common -I$(PVBUILDDIR)/VTK -I$(PVSRCDIR)/CoProcessing/CoProcessor

INCDIRS = $(PVINCDIRS) -I$(DECOMPDIR)/include

LIBDIRS = -L$(PVDIR)/lib/paraview-3.14 -L$(DECOMPDIR)/lib -L$(FFTW_PATH)/lib 
LIBS = -lvtkCoProcessorPythonD -l2decomp_fft -lmpi -lstdc++ -lfftw3

BIN = ns3d-pvcp-cells 

all: $(BIN)

FortranAdaptorAPI.o: $(PVSRCDIR)/CoProcessing/Adaptors/FortranAdaptors/FortranAdaptorAPI.cxx
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

VTKCellBasedDataSet.o: VTKCellBasedDataSet.cxx
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

NSadaptor.o: NSadaptor.f90 VTKCellBasedDataSet.o
	$(FC) -c $(FCFLAGS) $(INCDIRS) -c -o $@ $^

$(BIN): NavierStokes3DfftIMR.f90 NSadaptor.o FortranAdaptorAPI.o VTKCellBasedDataSet.o
	$(FC) -Wl,-rpath $(PVDIR)/lib/paraview-3.14 $(FCFLAGS) $(INCDIRS) $(LIBDIRS)  -o $@ $^ $(LIBS)


clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) $(BIN)
