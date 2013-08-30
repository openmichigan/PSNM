// 20130816 Mark Van Moer, NCSA
// Adaptor for ParaView/Catalyst Coprocessing.
// Interfaces with a Fortran glue module to get data from a Fortran
// simulation.

// Fortran specific header for API called from Fortran glue module
#include "vtkCPPythonAdaptorAPI.h"

// Adaptor specific headers, will need to change if simulation data structure
// changes. This should available from the ParaView install include directory
// if you chose to install development files during configuration.
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkCellData.h"

#include <string>
#include <iostream>
// Names are mangled with trailing underscore for linker visibility, along with
// using extern "C"

// Creates the vtkDataSet container for the CoProcessor. The simulation data
// must be compatible.
extern "C" void createcpimagedata_(int* nx, int* ny, int* nz)
{
	if (!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
		vtkGenericWarningMacro("Unable to access coprocessor data.");
		return;
	}

	vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();
	
	// Note: expects Fortran indexing
	img->SetExtent(0, *nx - 1, 0, *ny - 1, 0, *nz - 1);

	vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(img);
}

// Add data to the container.
extern "C" void addfield_(double* simdata, char* name)
{
	vtkSmartPointer<vtkCPInputDataDescription> idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input");
	vtkSmartPointer<vtkImageData> img = vtkImageData::SafeDownCast(idd->GetGrid());

	if (!img) {
		vtkGenericWarningMacro("No adaptor grid to attach field data to.");
		return;
	}

	std::string fieldPythonName = name;
	
	if (idd->IsFieldNeeded(name)) {
		vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
		field->SetName(name);
		field->SetArray(simdata, img->GetNumberOfPoints(), 1);
		img->GetPointData()->AddArray(field);
	}
}
