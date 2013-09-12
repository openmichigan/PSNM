// 20130912 Mark Van Moer, NCSA
// Adaptor for ParaView/Catalyst Coprocessing.
// Interfaces with a Fortran glue module to get data from a Fortran simulation.
// This works with the ParaView 4.0.1/Catalyst 1.0-alpha API.

// Fortran specific header, yes, even though it says Python.
#include "vtkCPPythonAdaptorAPI.h" 

// CoProcessor specific headers
// Routines that take the place of VTK dataset object creation.
// Called from Fortran code which also calls the Fortran Adaptor API
// supplied with ParaView source.
// Note: names mangled with trailing underscores for linker visibility.
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkCellData.h"
#include <string>

// These will be called from the Fortran "glue" code"
// Completely dependent on data layout, structured vs. unstructured, etc.
// since VTK/ParaView uses different internal layouts for each.

// Creates the data container for the CoProcessor.
// Takes the extents for both the global dataset and the particular subsection
// visible to the current MPI rank.
// Note: expects to receive Fortran base-1 indices.
extern "C" void createcpimagedata_(int* nx, int* ny, int* nz, int* xst, int* xen,
	int* yst, int* yen, int* zst, int* zen)
{
  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }

  // The simulation grid is a 3-dimensional topologically and geometrically 
  // regular grid. In VTK/ParaView, this is considered an image data set.
  vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();
  // Indexing based on change from F90 to C++, and also from nodal to cellular.
  img->SetExtent(*xst - 1, *xen, *yst - 1, *yen, *zst - 1, *zen); 
  
  // Setting spacing is important so that the camera position in the pipeline makes
  // sense if using different sized meshes between setting up the pipeline and running
  // the simulation. Origin can often be ignored.
  img->SetSpacing( 1.0 / *nx, 1.0 / *ny, 1.0 / *nz); // considering passing in as args. 
  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(img);
  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetWholeExtent(0, *nx, 0, *ny, 0, *nz);
}

// Add field to the data container.
// Separate from above because this could be dynamic, grid is static.
// Might be an issue, VTK assumes row major C storage.
// Underscore is by-hand name mangling for fortran linking.
// Note: Expects the address of the data, has no way of determining
// if the array is densely packed or not.
// Note 2: Assumes "name" points to null-terminated array of chars.
// Easiest way to do that is to concatenate in caller. 
extern "C" void addfield_(double* a, char* name)
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
    field->SetArray(a, img->GetNumberOfCells(), 1); 
    img->GetCellData()->AddArray(field);
  }
}
