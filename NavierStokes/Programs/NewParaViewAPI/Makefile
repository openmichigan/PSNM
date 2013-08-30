# beacon@NICS
# module load xeon/fftw/3.2.1
# This is for "offload" building.

# compilers, project uses C++ and Fortran
fc = ifort
#fc = gfortran
fcwarnings =
fcdebug =
fcflags = -O0 $(fcwarnings) $(fcdebug)

cxx = icpc
#cxx = g++
cxxwarnings =
cxxdebug =
cxxflags = -c $(cxxdebug) $(cxxwarnings)

pvinstalldir = /nics/b/home/mvanmoer/apps/ParaView/4.0.1
pvincdir = -I$(pvinstalldir)/include/paraview-4.0
pvlibs = -lvtkPVCatalystPython26D-pv4.0 -lvtkPVCatalyst-pv4.0 -lvtkCommonCore-pv4.0 -lvtkCommonDataModel-pv4.0 -lvtkCommonDataModel-pv4.0 -lvtkPVPythonCatalyst-pv4.0


incdirs = $(pvincdir)
libdirs = -L$(pvinstalldir)/lib/paraview-4.0

libs = $(FFTW3_LIB) $(pvlibs)

bin = ns
all: $(bin)

PointBasedDataSet.o: PointBasedDataSet.cxx
	$(cxx) $(cxxflags) $(incdirs) -o $@ $< $(libs)

NSadaptor.o: NSadaptor.f90 PointBasedDataSet.o
	$(fc) -c $(fcflags) -stand f03 $(incdirs) $(libdirs) -o $@ $^ $(libs)

navierstokes.o: navierstokes.f90 NSadaptor.o
	$(fc) -c $(fcflags) $< -o $@ $(libs)

$(bin): navierstokes.o NSadaptor.o PointBasedDataSet.o
	$(fc) -o $@ $(fcflags) $(incdirs) $(libdirs) -Wl,-rpath $(pvinstalldir)/lib/paraview-4.0 $^ $(libs)

clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) $(BIN)
