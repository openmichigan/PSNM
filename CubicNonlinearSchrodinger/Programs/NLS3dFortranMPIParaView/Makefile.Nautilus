# Nautilus
# Assumes: module load PE-intel
# Assumes: module load fftw/3.3.0
# Assumes: module load paraview/3.14.1

FC = ifort
CXX = icpc

DEBUG = -g
WARNINGS = 
FCFLAGS = $(WARNINGS) $(DEBUG) -O0
CXXFLAGS = $(DEBUG)

DECOMPDIR = /nics/b/home/mvanmoer/projects/muite/2decomp_fft

PVDIR = /sw/analysis/paraview/3.14.1/sles11.1_intel11.1.038
PVSRCDIR = $(PVDIR)/ParaView-3.14.1-Source
PVBUILDDIR = $(PVDIR)/build
PVINCDIRS = -I$(PVSRCDIR)/CoProcessing/Adaptors/FortranAdaptors -I$(PVBUILDDIR)/CoProcessing/Adaptors/FortranAdaptors -I$(PVSRCDIR)/VTK/Filtering -I$(PVSRCDIR)/VTK/Common -I$(PVBUILDDIR)/VTK -I$(PVSRCDIR)/CoProcessing/CoProcessor

INCDIRS = $(PVINCDIRS) -I$(DECOMPDIR)/include 
LIBDIRS = -L$(PVDIR)/lib/paraview-3.14 -L$(DECOMPDIR)/lib

LIBS = -lmpi -l2decomp_fft -lfftw3 -lvtkCoProcessorPythonD

BIN = NLSsplitting
all: $(BIN)

FortranAdaptorAPI.o: $(PVSRCDIR)/CoProcessing/Adaptors/FortranAdaptors/FortranAdaptorAPI.cxx
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

VTKCellBasedDataSet.o: VTKCellBasedDataSet.cxx
	$(CXX) -c $(CXXFLAGS) $(INCDIRS) -o $@ $<

NLSadaptor.o: NLSadaptor.f90 VTKCellBasedDataSet.o
	$(FC) -c $(FCFLAGS) $(INCDIRS) -o $@ $^

OBJECTS = FortranAdaptorAPI.o VTKCellBasedDataSet.o NLSadaptor.o

$(BIN): $(addsuffix .f90,$(BIN)) $(OBJECTS)
	$(FC) $(INCDIRS) $(LIBDIRS) -o $(BIN) $< $(LIBS) $(OBJECTS)

clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) $(BIN)
