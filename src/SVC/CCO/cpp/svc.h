// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdlib>
#include<vector>
#include<fstream>
#include<iostream>
#include<omp.h>

// VItA Libs
#include<structures/tree/AbstractCostEstimator.h>
#include<structures/domain/AbstractDomain.h>
#include<structures/domain/StagedDomain.h>
#include<core/GeneratorData.h>
#include<core/StagedFRROTreeGenerator.h>
#include<constrains/ConstantConstraintFunction.h>
#include<constrains/ConstantPiecewiseConstraintFunction.h>
#include<structures/tree/SingleVesselCCOOTree.h>
#include<structures/tree/AbstractObjectCCOTree.h>
#include<structures/vascularElements/AbstractVascularElement.h>
#include<io/VTKObjectTreeNodalWriter.h>
#include<structures/domain/DomainNVR.h>
#include<structures/domain/SimpleDomain.h>
#include<structures/domain/NormalDistributionGenerator.h>
#include<structures/vascularElements/SingleVessel.h>
#include<structures/tree/VolumetricCostEstimator.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>
#include<creators/ParallelepipedCreator.h>
#include<creators/CylinderCreator.h>
#include<structures/domain/AnnulusDistributionGenerator.h>
#include<structures/CCOCommonStructures.h>

// VTK libs
#include<vtkTransform.h>
#include<vtkCylinderSource.h>
#include<vtkDiskSource.h>
#include<vtkSmartPointer.h>
#include<vtkPolyDataWriter.h>
#include<vtkTriangleFilter.h>
#include<vtkTransformPolyDataFilter.h>
#include<vtkTransform.h>

/* int svc(int n, char* s); */
std::string SVC_coarse(const std::string configFilename, bool verbose, int nThreads = 1);
std::string SVC_macula(const std::string configFilename, bool verbose, int nThreads = 1);

