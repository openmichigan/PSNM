#
#  Copyright 2012 NVIDIA Corporation
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

F90 = pgf90
OMPFLAGS = -O3 -Mcuda -Minfo 
CUDA_PATH = /usr/local/pgi/linux86-64/2012/cuda/4.1
WARNINGS = -Minform=inform
CCFLAGS  = -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcudart -lcufft

CC = g++
CXXFLAGS = -c

PARAVIEWSRCDIR = /usr/apps/vis/ParaView/ParaView-3.14.1-Source
PVFORTRANAPIDIR = $(PARAVIEWSRCDIR)/CoProcessing/Adaptors/FortranAdaptors

# OSMesa build, change these after X installed on Forge
PARAVIEWBUILDDIR = /usr/apps/vis/ParaView/ParaView-3.14.1-OSMesa-build
PARAVIEWDIR = /usr/apps/vis/ParaView/3.14.1-OSMesa/gnu/ParaViewServers-3.14.1-Linux-x86_64


INCDIRS = -I$(PARAVIEWSRCDIR)/CoProcessing/Adaptors/FortranAdaptors -I$(PARAVIEWSRCDIR)/CoProcessing/CoProcessor -I$(PARAVIEWSRCDIR)/VTK/Filtering -I$(PARAVIEWSRCDIR)/VTK/Common -I$(PARAVIEWBUILDDIR)/CoProcessing/Adaptors/FortranAdaptors -I$(PARAVIEWBUILDDIR)/VTK

LIBDIRS = -L$(PARAVIEWDIR)/lib/paraview-3.14 -L/usr/apps/vis/Mesa/7.11.2/lib

LIBS = -lvtkCoProcessorPythonD -lOSMesa
 
BIN =  ns2dcn

all: $(BIN)

FortranAdaptorAPI.o: $(PVFORTRANAPIDIR)/FortranAdaptorAPI.cxx
	$(CC) $(CXXFLAGS) $(INCDIRS) -o $@ $<

ns2dcnVTKDataSet.o: ns2dcnVTKDataSet.cxx
	$(CC) $(CXXFLAGS) $(INCDIRS) -o $@ $<

ns2dcnadaptor.o: ns2dcnadaptor.f90 ns2dcnVTKDataSet.o
	$(F90) $(CXXFLAGS) -o $@ $^

ns2dcn: navierstokes.cuf ns2dcnadaptor.o FortranAdaptorAPI.o ns2dcnVTKDataSet.o
	$(F90) $(WARNINGS) $(OMPFLAGS) $(CCFLAGS) $(LIBDIRS) -o $@ $^ $(LIBS)


clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) $(BIN)
