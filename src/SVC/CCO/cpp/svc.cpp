#include"svc.h"

void MakeDisk(const std::string filename,
	      const double innerRadius,
	      const double outerRadius
	      )
{
  vtkSmartPointer<vtkDiskSource> diskSource = vtkSmartPointer<vtkDiskSource>::New();
  // diskSource->SetCenter(0.0,0.0,-0.1); // Will be (0,0,0) by default.
  diskSource->SetInnerRadius(innerRadius);
  diskSource->SetOuterRadius(outerRadius);
  diskSource->SetCircumferentialResolution(40);
  diskSource->Update();

  // vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  // transform->RotateX(90.0);
  // vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  // transformFilter->SetInputConnection(diskSource->GetOutputPort());
  // transformFilter->SetTransform(transform);
  // transformFilter->Update();
   
  
  vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
  triangleFilter->SetInputData(diskSource->GetOutput());
  triangleFilter->Update();

  //Write in file
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(triangleFilter->GetOutput());
  writer->SetFileTypeToASCII();
  writer->Write();
}


void MakeCylinder(const std::string filename,
		  const double radius
		  )
{
  vtkSmartPointer<vtkCylinderSource> cylinderSource = vtkSmartPointer<vtkCylinderSource>::New();
  cylinderSource->SetCenter(0.0,0.0,-0.1);
  cylinderSource->SetHeight(0.1);
  cylinderSource->SetRadius(radius);
  cylinderSource->SetResolution(40);
  cylinderSource->Update();

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->RotateX(90.0);
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetInputConnection(cylinderSource->GetOutputPort());
  transformFilter->SetTransform(transform);
  transformFilter->Update();
  
  
  vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
  triangleFilter->SetInputData(transformFilter->GetOutput());
  triangleFilter->Update();

  //Write in file
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(triangleFilter->GetOutput());
  writer->SetFileTypeToASCII();
  writer->Write();

  // delete cylinderSource;
  // delete writer;
  // delete triangleFilter;
  // delete transformFilter;
  // delete transform;
}

std::string SVC_coarse(const std::string configFilename, bool verbose, int nThreads)
{

  // Set number of threads
  if (nThreads > 0)
    {
      omp_set_dynamic(0);
      omp_set_num_threads(5);
      // cout << omp_get_num_threads() << " threads running with ";
      // cout << omp_get_max_threads() << " threads available." << endl; 
    } 
  
  std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
  std::ofstream   fout("/dev/null");
  if (!verbose)
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
  
  std::string line(configFilename);
  std::cout << "Reading file: " << line << std::endl;

  std::ifstream config;
  config.open(configFilename);
  if (!config){
    std::cerr << "Error: parameter file " << configFilename
	      << " could not be opened." << std::endl;
    std::exit(1);
  }

  // Output file
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  std::string outputFilename {line};

  // Copying the config file in the output folder
  std::string root_results {outputFilename.substr(0, outputFilename.find_last_of("/"))};
  std::string command {"mkdir -p " + root_results};
  std::cout << command << std::endl;
 
  const int dir_err = system(command.c_str());
  if (-1 == dir_err)
    {
      std::cout << "Error creating directory!" << std::endl;
      std::exit(1);
    }

  command = {"cp " + configFilename + " " + outputFilename + ".conf"};
  std::cout << "COPY COMMAND " << command << std::endl;
  const int copy_err = system(command.c_str());
  if (-1 == copy_err)
    {
      std::cout << "Error copying the config file!" << std::endl;
      std::exit(1);
    }


  // Root tree (.cco)
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  std::string rootTreeFilename {line};

  // Simulation parameters
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int n {std::stoi(line)}; // Terminal vessels

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double lLimFr {stod(line)}; // Correction step factor

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *gamma {new ConstantConstraintFunction<double, int>(stod(line))}; // Murray's law exponent

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *delta {new ConstantConstraintFunction<double, int>(stod(line))}; // Symmetry ratio

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *eta {new ConstantConstraintFunction<double, int>(stod(line))}; // viscosity in cP

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double viscTolerance {stod(line)}; // Viscosity tolerance

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double thetaMin {stod(line) * M_PI}; // Minimum bifurcation angle

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double perfAreaFr {stod(line)}; // Correction step factor

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double closeNeighFr {stod(line)}; // Close neighborhood factor  

  // Consecutive attempts to generate a point - nFail
  int nFail = 20;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int DeltaNu = 7;
  // Buffer size for random point generation
  int nDraw {1000};
  // Random seed
  long long int seed {std::time(nullptr)};

  // // // // Running the CCO
  StagedDomain *stagedDomain = new StagedDomain();
  AdimSproutingVolumetricCostEstimator *FSprout = new AdimSproutingVolumetricCostEstimator(50, 0.5, 1e+4, 40, 70*1e-6);
  AbstractCostEstimator *costEstimator = FSprout;
  GeneratorData *generatorData = new GeneratorData(100, nFail, lLimFr, perfAreaFr, closeNeighFr, 0.1, DeltaNu, 0, false, costEstimator);

  ifstream is_root_tree_correct {rootTreeFilename};
  if (!is_root_tree_correct){
    std::cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << std::endl;
    std::exit(1);
  }
  is_root_tree_correct.close();

  // Creates a temporary domain (in vtk format), centered around the fovea.
  SingleVesselCCOOTree *rootTree = new SingleVesselCCOOTree(rootTreeFilename, generatorData, gamma, delta, eta);

  point opticDisc = rootTree->getRoot()->getVessels()[0]->xDist;
  AnnulusDistributionGenerator *dist = new AnnulusDistributionGenerator(0.02, opticDisc.p, -0.5,4);

  double lb[3] {-2,-2,-0.1}, ub[3] {2.0,2.0,0.0}; // Assume (0,0) is the fovea.
  // ParallelepipedCreator *para = new ParallelepipedCreator(lb, ub);
  // std::vector<double> fovea{0.0,0.0,-0.1};
  // CylinderCreator *para = new CylinderCreator(fovea, 1.5, 0.2, 50);
  // para->create(outputFilename + "_tmp.vtp");
  // delete para;
  // MakeCylinder(outputFilename + "_tmp.vtp", 1.5);
  MakeDisk(outputFilename + "_tmp.vtp", 0.0, 1.5); // WILL IT WORK IN 2D? 
  
  // Load the temporary domain
  SimpleDomain *domain = new SimpleDomain(outputFilename + "_tmp.vtp", nDraw, seed, generatorData, dist);
  domain->setIsConvexDomain(true);
  domain->setIsBifPlaneContrained(false);
  domain->setMinBifurcationAngle(thetaMin);
  domain->setOpticDisc(opticDisc);
  stagedDomain->addStage(n+1, domain);
  stagedDomain->setOpticDisc(opticDisc);

  std::ofstream f;
  f.open(outputFilename+"_RandomPoints.dat");
  std::vector<point> points = dist->getNPoints(10000);
  for (point p : points)
    {
      f << p.p[0] << ' ' << p.p[1] << " " << p.p[2] << endl;
    }
  f.close();

  VTKObjectTreeNodalWriter *treeWriter = new VTKObjectTreeNodalWriter();
  rootTree->save(outputFilename + ".cco.root");
  treeWriter->write(outputFilename + ".vtp.root", rootTree);

  rootTree->setIsInCm(true);
  rootTree->setCurrentStage(1);

  long long int nTermTotal = rootTree->getNTerms()+n;

  StagedFRROTreeGenerator *tree_generator = new StagedFRROTreeGenerator(stagedDomain, rootTree, nTermTotal, {gamma}, {delta}, {eta});

  SingleVesselCCOOTree *tree = static_cast<SingleVesselCCOOTree *>(tree_generator->resume(200, "./"));

  tree->save(outputFilename + ".cco");
  treeWriter->write(outputFilename + ".vtp", tree);


  
  if (!verbose)
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer

  
  
  delete rootTree;
  // delete tree; /// Causes segfault because tree may not be a dynamic pointer? IDK 
  delete tree_generator;
  delete domain;
  delete treeWriter;
  delete dist;
  delete stagedDomain;
  delete FSprout;
  delete generatorData;
  delete gamma;
  delete eta;
  delete delta;

  return outputFilename + ".cco";
}
  
		 
std::string SVC_macula(const std::string configFilename, bool verbose, int nThreads)
{

  // Set number of threads
  if (nThreads > 0)
    {
      omp_set_dynamic(0);
      omp_set_num_threads(5);
      cout << omp_get_num_threads() << " threads running with ";
      cout << omp_get_max_threads() << " threads available." << endl; 
    } 

  std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
  std::ofstream   fout("/dev/null");		 // Empty buffer
  if (!verbose)
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
  
  std::string line(configFilename);
  std::cout << "Reading file: " << line << std::endl;

  std::ifstream config;
  config.open(configFilename);
  if (!config){
    std::cerr << "Error: parameter file " << configFilename
	      << " could not be opened." << std::endl;
    std::exit(1);
  }

  // Output file
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  std::string outputFilename {line};

  // Copying the config file in the output folder
  std::string root_results {outputFilename.substr(0, outputFilename.find_last_of("/"))};
  std::string command {"mkdir -p " + root_results};
  std::cout << command << std::endl;
 
  const int dir_err = system(command.c_str());
  if (-1 == dir_err)
    {
      std::cout << "Error creating directory!" << std::endl;
      std::exit(1);
    }

  command = {"cp " + configFilename + " " + outputFilename + ".conf"};
  std::cout << "COPY COMMAND " << command << std::endl;
  const int copy_err = system(command.c_str());
  if (-1 == copy_err)
    {
      std::cout << "Error copying the config file!" << std::endl;
      std::exit(1);
    }

  // Root tree (.cco)
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  std::string rootTreeFilename {line};

  // Simulation parameters
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int n1 {std::stoi(line)}; // Terminal vessels stage 1
  getline(config, line);
  int n2 {std::stoi(line)}; // Terminal vessels stage 2

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double lLimFr {stod(line)}; // Correction step factor

  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *gamma {new ConstantConstraintFunction<double, int>(stod(line))}; // Murray's law exponent
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *delta {new ConstantConstraintFunction<double, int>(stod(line))}; // Symmetry ratio

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *eta {new ConstantConstraintFunction<double, int>(stod(line))}; // viscosity in cP
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double viscTolerance {stod(line)}; // Viscosity tolerance

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double thetaMin {stod(line) * M_PI}; // Minimum bifurcation angle

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double perfAreaFr {stod(line)}; // Correction step factor

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double closeNeighFr {stod(line)}; // Close neighborhood factor  
  
  // Consecutive attempts to generate a point - nFail
  int nFail = 20;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int DeltaNu = 20;
  // Buffer size for random point generation
  int nDraw {1000};
  // Random seed
  long long int seed {std::time(nullptr)};

  // // // // Running the CCO
  StagedDomain *stagedDomain = new StagedDomain();
  AdimSproutingVolumetricCostEstimator *FSprout = new AdimSproutingVolumetricCostEstimator(50, 0.5, 1e+4, 40, 70*1e-6);
  AbstractCostEstimator *costEstimator = FSprout;
  GeneratorData *generatorData = new GeneratorData(20, nFail, lLimFr, perfAreaFr, closeNeighFr, 0.1, DeltaNu, 0, false, costEstimator);

  ifstream is_root_tree_correct {rootTreeFilename};
  if (!is_root_tree_correct){
    std::cerr << "Error: file '" << rootTreeFilename
	      << " could not be opened. The root tree's .cco file could not be found." << std::endl;
    std::exit(1);
  }
  is_root_tree_correct.close();

  // Creates a temporary domain (in vtk format), centered around the fovea.
  SingleVesselCCOOTree *rootTree = new SingleVesselCCOOTree(rootTreeFilename, generatorData, gamma, delta, eta);

  std::cout << "DP = " << rootTree->getDp() << " --- ROOT RADIUS = " << rootTree->getRoot()->getDistalRadius() << " vs " << rootTree->getRootRadius() << std::endl;
  double lb[3] {-2.0,-1.5,-0.1}, ub[3] {1.5,1.5,0.0}; // Assume (0,0) is the fovea.
  point opticDisc = rootTree->getRoot()->getVessels()[0]->xDist;
  opticDisc.p[2] = 0.0;
  AnnulusDistributionGenerator *dist = new AnnulusDistributionGenerator(0.02, opticDisc.p, -0.5,4);

  // Make the domain for stage 1
  // MakeCylinder(outputFilename+"_tmp1.vtp", 0.4);
  MakeDisk(outputFilename+"_tmp1.vtp", 0.0, 0.3);

  SimpleDomain *domain1 = new SimpleDomain(outputFilename + "_tmp1.vtp", nDraw, seed, generatorData, dist);
  domain1->setIsConvexDomain(true);
  domain1->setIsBifPlaneContrained(false);
  domain1->setMinBifurcationAngle(thetaMin);
  domain1->setOpticDisc(opticDisc);
  stagedDomain->addStage(n1, domain1);
  stagedDomain->setOpticDisc(opticDisc);
  
  // Make the domain for stage 2
  // MakeCylinder(outputFilename+"_tmp2.vtp", 0.2);
  MakeDisk(outputFilename+"_tmp2.vtp", 0.3, 0.6);

  SimpleDomain *domain2 = new SimpleDomain(outputFilename + "_tmp2.vtp", nDraw, seed, generatorData, dist);
  // domain2->setIsConvexDomain(true);
  domain2->setIsBifPlaneContrained(false);
  domain2->setMinBifurcationAngle(thetaMin);
  domain2->setOpticDisc(opticDisc);
  stagedDomain->addStage(n2, domain2);
  stagedDomain->setOpticDisc(opticDisc);

  std::ofstream f;
  f.open(outputFilename+"_RandomPoints.dat");
  for (int i = 0; i < 10000; i++)
    {
      point p = domain2->getRandomPoint();
      f << p.p[0] << ' ' << p.p[1] << " " << p.p[2] << endl;
    }
  f.close();

  VTKObjectTreeNodalWriter *treeWriter = new VTKObjectTreeNodalWriter();
  rootTree->setCurrentStage(2);
  rootTree->save(outputFilename + ".cco.root");
  treeWriter->write(outputFilename + ".vtp.root", rootTree);

  rootTree->setIsInCm(true);

  long long int nTermTotal = rootTree->getNTerms()+n1+n2;
  std::cout << "DP = " << rootTree->getDp() << " --- " << rootTree->getRoot()->getDistalRadius() << std::endl;

  // stagedDomain->setInitialStage(2);


  StagedFRROTreeGenerator *tree_generator = new StagedFRROTreeGenerator(stagedDomain, rootTree, nTermTotal, {gamma, gamma}, {delta, delta}, {eta, eta});
  tree_generator->setDLim(0.01);

  if (!verbose)
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
  SingleVesselCCOOTree *tree = static_cast<SingleVesselCCOOTree *>(tree_generator->resume(200, "./"));
  if (!verbose)
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
  

  tree->save(outputFilename + ".cco");
  treeWriter->write(outputFilename + ".vtp", tree);

  if (!verbose)
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer

  delete rootTree;
  // delete tree; /// Causes segfault because tree may not be a dynamic pointer? IDK 
  delete tree_generator;
  delete domain1;
  delete domain2;
  delete treeWriter;
  delete dist;
  delete stagedDomain;
  delete FSprout;
  delete generatorData;
  delete gamma;
  delete eta;
  delete delta;

  return outputFilename + ".cco";
}
