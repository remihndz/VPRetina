#include "transitTimeDistribution.hpp"
#include "argumentParser.hpp"

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  
  size_t cutoff;
  std::vector<std::string> graphFiles;
  commandLineParser(argc, argv, cutoff, graphFiles);
  
  std::ofstream fout;
  std::ifstream graphFile;
  double tt, pr, len, oef;
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
    std::string fileName = graphFileName;

    // NOTE: This can be read using numpy.fromfile(fileName.pathdata.bin, dtype=numpy.float32)
    // By default, numpy.fromfile looks for DOUBLE which causes undefined behaviour!!
    fileName.replace(fileName.find(".graph"), std::string(".graph").length(), ".pathdata.bin"); // Replace .graph with .pathdata.bin
    fout.open(fileName);
    struct pathInfo{
      float tt, pr, oef, len;
    } p;
    BOOST_FOREACH(boost::tie(tt, pr, oef, len), boost::combine(pathsData.pathsTransitTimes, pathsData.pathsProbabilities, pathsData.pathsOEF, pathsData.pathsLength)){
      if (true){ // (!((std::isnan(tt)) || (std::isinf(tt)) || (std::isnan(pr)) || (std::isinf(pr)))){
	p.tt  = tt;
	p.pr  = pr;
	p.oef = oef;
	p.len = len;
	fout.write(reinterpret_cast<char *>(&p), sizeof(p));
      }
    }
    fout.close();

    // // This writes as ascii (up to 6 significant digits?)
    // fileName = graphFileName;

    // fileName.replace(fileName.find(".graph"), std::string(".graph").length(), ".pathdata");
    // fout.open(fileName);
    
    // BOOST_FOREACH(boost::tie(tt, pr), boost::combine(pathsData.pathsTransitTimes, pathsData.pathsProbabilities)){
    //   if (!((std::isnan(tt)) || (std::isinf(tt)))){
    // 	fout << tt << "," << pr << '\n';
    //   }
    // }
    // fout.close();
  }  
  return 0;
}
