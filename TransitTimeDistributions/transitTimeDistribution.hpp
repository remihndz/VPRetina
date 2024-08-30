#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <iostream>
#include <vector>
#include <stack>
#include <utility>
#include <fstream>
#include <map>
#include <deque>
#include <boost/lexical_cast.hpp>
#include <boost/range/combine.hpp>
#include <boost/foreach.hpp>
#include <cmath> // For isnan function and exponential
#include <boost/program_options.hpp> // Argument parsing

const double PI = 3.1415926535897932384626433832795028841971693993751058209;
const double CompartmentVesselsLength = 0.36; // In cm the estimated length of the compartment based on a radius of 20um and resistance of 1e6:
                                              // L = (1e6*pi*r^4)/(8*mu(r,hd=0.45)) 


typedef std::pair<size_t, size_t> State;
struct Vertex{
  std::vector<size_t> adj; // Nodes it is connected to
  std::vector<double> flow, rad, ht, len,
    tt, proba; // Edge properties
};

struct Graph{
  std::vector<Vertex> vertices;
  size_t size(){return vertices.size();}
  size_t CRA, CRV; 
};

struct pathAnalysis {
  std::vector<double> pathsTransitTimes;
  std::vector<double> pathsProbabilities;
  std::vector<size_t> pathsLength;
  size_t numberOfPaths(){return pathsTransitTimes.size();}
};

pathAnalysis dfs(Graph  &graph, size_t start, size_t end, size_t cutoff);
std::vector<size_t> getLeaves(const Graph &graph);
void GetEdgeData(Graph &graph);
std::vector<size_t> getRoots(const Graph &graph);
void print(const std::vector<size_t> &path);
void print(const std::deque<size_t> &path);
void print(const std::deque<double> &path);
void print(const std::deque<float> &path);
void ReadGraph(std::ifstream& graphFile, Graph& graph);

#endif
