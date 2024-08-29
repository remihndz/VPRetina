#include "transitTimeDistribution.hpp"
#include "argumentParser.hpp"

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  
  size_t cutoff;
  std::vector<std::string> graphFiles;
  commandLineParser(argc, argv, cutoff, graphFiles);
  
  std::ofstream fout, foutLengths;
  std::ifstream graphFile;
  float tt, pr, oef;
  size_t len;
  for (auto graphFileName: graphFiles){

    std::cout << "Analysing " << graphFileName << " ... " << std::endl;
    Graph graph;
    graphFile.open(graphFileName);

    ReadGraph(graphFile, graph);
    graphFile.close();
    std::cout << "Graph read with " <<  graph.size() << " edges. CRA=" << graph.CRA
	 << "; CRV=" << graph.CRV << std::endl;
    GetEdgeData(graph);         
    auto pathsData = dfs(graph, graph.CRA, graph.CRV, cutoff);

    // // Write the graph data in binary (to save space, these arrays are pretty large)
    std::cout << std::flush;
    std::string fileName = graphFileName,
      fileNameLengths = graphFileName;

    // NOTE: This can be read using numpy.fromfile(fileName.pathdata.bin, dtype=numpy.float32)
    // By default, numpy.fromfile looks for DOUBLE which causes undefined behaviour!!
    fileName.replace(fileName.find(".graph"), std::string(".graph").length(), ".pathdata.bin"); // Replace .graph with .pathdata.bin
    fileNameLengths.replace(fileNameLengths.find(".graph"), std::string(".graph").length(), ".pathlength.dat"); // Replace .graph with .pathlength.dat
    fout.open(fileName);
    foutLengths.open(fileNameLengths);
    
    struct pathInfo{
      float tt, pr, oef;
    } p;
    BOOST_FOREACH(boost::tie(tt, pr, oef, len), boost::combine(pathsData.pathsTransitTimes, pathsData.pathsProbabilities, pathsData.pathsOEF, pathsData.pathsLength)){
      if (true){ // (!((std::isnan(tt)) || (std::isinf(tt)) || (std::isnan(pr)) || (std::isinf(pr)))){
	p.tt  = tt;
	p.pr  = pr;
	p.oef = oef;
	fout.write(reinterpret_cast<char *>(&p), sizeof(p));
	foutLengths << len << std::endl;
      }
    }
    fout.close();
    foutLengths.close();
  }  
  return 0;
}
