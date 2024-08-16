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
#include <cmath> // For isnan function

const double PI = 3.1415926535897932384626433832795028841971693993751058209;
const size_t CUTOFF = 300;
// TODO: give path vectors a fixed size (cutoff) and keep track of depth to save some time on checks from pop_back

using namespace std;

typedef vector<vector<double>> doubleMat;
typedef vector<vector<size_t>> intMat;

struct pathAnalysis {
  vector<double> pathsTransitTimes;
  vector<double> pathsProbabilities;
  
  size_t numberOfPaths(){return pathsTransitTimes.size();}
};

void print(const vector<size_t> &path){
  for (auto i: path)
    cout << i << "->";
  cout << endl;
};

void print(const deque<size_t> &path){
  for (auto i: path)
    cout << i << "->";
  cout <<endl;
};

enum State {Unvisited, Visited, Visiting};
struct Vertex{
  vector<size_t> adj; // Nodes it is connected to
  State state=Unvisited;
  vector<double> flow, rad, ht, len,
    tt, proba; // Edge properties
};

// typedef vector<Vertex> Graph;
struct Graph{
  vector<Vertex> vertices;
  size_t size(){return vertices.size();}
  size_t CRA, CRV; 
};

vector<size_t> getRoots(const Graph &graph){
  vector<size_t> roots;
  vector<bool> isRoot(graph.vertices.size(), true);
  for (size_t u=0; u<graph.vertices.size(); u++){
    for (size_t v: graph.vertices[u].adj)
      isRoot[v] = false;
  }
  for (size_t u=0; u<graph.vertices.size(); u++)
    if (isRoot[u]){
      roots.push_back(u);
    }      
  return roots;
}

vector<size_t> getLeaves(const Graph &graph){
  vector<size_t> leaves;
  size_t u = 0; // u will be the vertex index in adjacency matrix
  for (const Vertex &v: graph.vertices){
    if (v.adj.empty()) // No neighbors
      leaves.push_back(u); 
    u++;
  }
  return leaves;
}

void GetEdgeData(Graph &graph){
  // Compute the probability of RBC at node u going to node v, where v is a downstream neighbor of u.
  
  // Compute the probability of a RBC going from vertex u to v
  // and the transit time between u and v
  double totalRBC, rbc;
  int u=0;
  for (Vertex &v: graph.vertices){ // Iterate through vertices v
    v.proba.resize(v.adj.size());
    v.tt.resize(v.adj.size());
    totalRBC = 0;
    u++;
    for (size_t k=0; k<v.adj.size(); k++){ // Iterate through neighbors number k
      rbc = v.flow[k]*v.ht[k]; // Flow of RBC from u to v
      v.tt[k]    = PI*v.rad[k]*v.rad[k]*v.len[k]/v.flow[k];
      v.proba[k] = rbc;
      totalRBC  += rbc;
    }
    for (size_t k=0; k<v.adj.size(); k++)
      v.proba[k] /= totalRBC; // Scale to a probability
  }      
}

// pathAnalysis BFS(Graph &graph, size_t start, size_t end, size_t cutoff){
//   intMat &adj = graph.vertices; // Graph adjacency matrix
//   doubleMat &edgeProbabilities = graph.edgeProbabilities,
//     &edgeTransitTimes = graph.edgeTransitTimes;
//   //initialize:
//   pathAnalysis pathData;
//   vector<bool> visited(graph.size(), false);
// }

// //graph[i][j] stores the j-th neighbour of the node i
pathAnalysis dfs(Graph &graph, size_t start, size_t end, size_t cutoff) 
{
  //initialize:
  pathAnalysis pathData;
  //remember the node (first) and the index of the next neighbour (second)
  typedef pair<size_t, size_t> State;
  stack<State> to_do_stack;
  deque<size_t> path; //remembering the way
  vector<bool> visited(graph.size(), false); //caching visited - no need for searching in the path-vector
  deque<double> pathTransitTime(cutoff), pathProbability(cutoff); // Remembering the probability and transit time of the path. Use the same way as `path` but stores (cumulating) times/probabilities instead of nodes.

  // // Pre-allocate max size. Note: this keeps size() at 0.
  // path.reserve(cutoff);
  // pathTransitTime.reserve(cutoff);
  // pathProbability.reserve(cutoff);
  
  //start in start!
  to_do_stack.push(make_pair(start, 0));
  visited[start]=true;    
  path.push_back(start);
  pathTransitTime.push_back(0);
  pathProbability.push_back(1);

  size_t pathCount = 0, oneMil = static_cast<size_t>(1e6);
  while (!to_do_stack.empty())
    {
      State &current = to_do_stack.top();//current stays on the stack for the time being...
      
      if (current.first == end || current.second == graph.vertices[current.first].adj.size() || path.size()>cutoff)//goal reached or done with neighbours?
	{
          if (current.first == end)
	    {
	      pathCount+=1;
	      pathData.pathsProbabilities.push_back(pathProbability.back());
	      pathData.pathsTransitTimes.push_back(pathTransitTime.back());
	      if (pathCount % oneMil == 0){
		cout << "\rNumber of paths found (in millions): " << pathCount / oneMil << flush;
	      }
	      // print(path);//found a way!      
	    }
          //backtrack:
          visited[current.first]=false;//no longer considered visited
	  //go a step back	
          path.pop_back();
	  pathProbability.pop_back();
	  pathTransitTime.pop_back();
          to_do_stack.pop();//no need to explore further neighbours   
	}
      else{//normal case: explore neighbours
	size_t next=graph.vertices[current.first].adj[current.second];
	current.second++;//update the next neighbour in the stack!
	if(!visited[next]){
	  //putting the neighbour on the todo-list
	  to_do_stack.push(make_pair(next, 0));
	  visited[next]=true;
	  path.push_back(next);
	  pathTransitTime.push_back(pathTransitTime.back()+graph.vertices[current.first].tt[current.second-1]);
	  pathProbability.push_back(pathProbability.back()*graph.vertices[current.first].proba[current.second-1]);
	  // print(path);
	}      
      }
    }
  cout << endl;
  return pathData;
}

void ReadGraph(ifstream& graphFile, Graph& graph)
{
  // The nodes in the text file might be strings. We will map them to consecutive integers
  map<string, size_t> nodesToInt;
  string line;
  size_t nv;
  string delimiter;
  string sCRA, sCRV; // The node names as strings
  
  while (getline(graphFile, line)){
    // Be careful that the reader doesn't read another
    // graph attribute that ends in "CRA:" or "CRV:"
    // by keeping a blank space before, e.g., " CRA:"
    if (line.find(string(" CRA:")) != string::npos){
      delimiter = ": ";
      line.erase(0, line.find(delimiter)+delimiter.length());
      sCRA = line;
    }
    else if (line.find(string(" CRV:")) != string::npos){
      delimiter = ": ";
      line.erase(0, line.find(delimiter)+delimiter.length());
      sCRV = line;
    }
    else if (line.find(string("# Nodes")) != string::npos){ // If found the '# Nodes: nNodes' line
      delimiter = ": ";
      line.erase(0, line.find(delimiter)+delimiter.length()); // Removes the '# Nodes: ' before the number of nodes
      // cout << "Number of vertices " << line << endl;
      nv = boost::lexical_cast<size_t>(line);
      getline(graphFile, line); // Header line we could sparse
      break;
    }
  }

  graph.vertices.resize(nv); // Initialize with nv vertices
  
  // Make the map of nodes to integer
  delimiter = ",";
  for (int i=0; i<nv; i++){
    getline(graphFile, line);
    string node = line.substr(0, line.find(delimiter)); // Node name
    nodesToInt[node] = i;    
  }
  graph.CRA = nodesToInt[sCRA];
  graph.CRV = nodesToInt[sCRV];
  
  // Get the edges
  size_t ne;
  string edgeAttributeNames;
  int flowAttributePos=0, htAttributePos=0,
    radAttributePos=0, lenAttributePos=0; // That's how many comas we need to look for to find the column wiht flow values
  while (getline(graphFile, line)){
      if (line.find(string("# Edges")) != string::npos){
        delimiter = ": ";
	line.erase(0, line.find(delimiter)+delimiter.length()); // Removes the '# Edges: ' before the number of edges
	// cout << "Number of edges " << line << endl;
	ne = stoi(line);      
	getline(graphFile, edgeAttributeNames); // header line we could sparse
	string edgeAttribute;
	delimiter = ',';
	edgeAttributeNames.append(delimiter);
	int k =0;
	while(edgeAttributeNames.find(delimiter) != string::npos){
	  edgeAttribute = edgeAttributeNames.substr(0, edgeAttributeNames.find(delimiter));
	  //cout << "Edge attribute: " << edgeAttribute << endl;
	  edgeAttributeNames.erase(0, edgeAttributeNames.find(delimiter) + delimiter.length());
	  if (edgeAttribute.compare("flow")==0)
	    {
	      //cout << "Found flow columns at " << k << endl;
	      flowAttributePos = k;
	    }
	  else if (edgeAttribute.compare("hd")==0)
	    {
	      //cout << "Found ht columns at " << k << endl;
	      htAttributePos = k;
	    }
	  else if (edgeAttribute.compare("radius")==0)	    
	    {
	      //cout << "Found radius columns at " << k << endl;
	      radAttributePos = k;
	    }
	  else if (edgeAttribute.compare("length")==0){
	    // cout << "Found length columns at " << k << endl;
	    lenAttributePos = k;
	  }
	  k++;
	}
	break;
      }
  }

  // Read edge data
  size_t u,v;
  double f, h, r, l;
  delimiter = ',';
  string nodeName;
  for (int j=0; j<ne; j++){
    getline(graphFile, line);
    nodeName = line.substr(0, line.find(delimiter));
    u = nodesToInt[nodeName]; // Parent node
    line.erase(0, line.find(delimiter)+delimiter.length()); // Remove from string
    nodeName = line.substr(0, line.find(delimiter));
    v = nodesToInt[nodeName]; // Branch node
    line.erase(0, line.find(delimiter)+delimiter.length()); // Remove from string
    
    // Find the column that has flow value for the vessel
    int k=2; // Starts at one because we remove two columns (u and v the node names)
    line.append(delimiter); // In case flow or hematocrit is the last column
    while (line.find(delimiter) != string::npos){
      if (k==flowAttributePos)
	f = boost::lexical_cast<double>(line.substr(0, line.find(delimiter)));
      else if (k==htAttributePos){
	try {h = boost::lexical_cast<double>(line.substr(0, line.find(delimiter)));}
	catch (...) {h = 0.45;}
      }
      else if (k==radAttributePos){
	r = boost::lexical_cast<double>(line.substr(0, line.find(delimiter)));
      }
      else if (k==lenAttributePos){
	try {	  
	  l = boost::lexical_cast<double>(line.substr(0, line.find(delimiter)));}
	catch (...) {l = 1e20;}
      }
      line.erase(0, line.find(delimiter)+delimiter.length());
      k++;
    }

    if ((u==graph.CRA || v==graph.CRV) && f<0)
      cout << "CRA or CRV orientation is being switched!" << endl;
    
    if (f>=0){
      graph.vertices[u].adj.push_back(v); // Add the edge
      graph.vertices[u].ht.push_back(0.45);
      graph.vertices[u].flow.push_back(f);
      graph.vertices[u].rad.push_back(r);
      graph.vertices[u].len.push_back(l);
      // cout << "Edge (" << u << "," << v << ") with flow " << graph[u].flow.back() << " and ht " << graph[u].ht.back() << endl;
    }
    else{ // Reverse the edge to follow flowx
      // graph.vertices[v].adj.push_back(u);
      graph.vertices[u].adj.push_back(v); // Add the edge
      graph.vertices[v].ht.push_back(0.45);
      graph.vertices[v].flow.push_back(-1.0*f);
      graph.vertices[v].rad.push_back(r);
      graph.vertices[v].len.push_back(l);
      // cout << "Original edge was reversed to become (" << v << "," << u << ") with flow " << graph.vertices[v].flow.back() << " and ht " << graph.vertices[v].ht.back() << endl;
    }
  }
};

int main(int argc, char *argv[])
{
  
  ofstream fout;
  ifstream graphFile;
  double tt, pr;
  for (int i=1; i<argc; i++){

    cout << "Analysing " << argv[i] << " ... " << endl;
    Graph graph;
    graphFile.open(argv[i]);

    ReadGraph(graphFile, graph);
    graphFile.close();
	    
    GetEdgeData(graph);         
    auto pathsData = dfs(graph, graph.CRA, graph.CRV, CUTOFF);

    // // Write the graph data in binary (to save space, these arrays are pretty large)
    cout << flush;
    string fileName = argv[i];


    // NOTE: This can be read using numpy.fromfile(fileName.pathdata.bin, dtype=numpy.float32)
    // By default, numpy.fromfile looks for DOUBLE which causes undefined behaviour!!
    fileName.replace(fileName.find(".graph"), string(".graph").length(), ".pathdata.bin");
    fout.open(fileName);
    struct pathInfo{
      float tt, pr;
    } p;
    BOOST_FOREACH(boost::tie(tt, pr), boost::combine(pathsData.pathsTransitTimes, pathsData.pathsProbabilities)){
      if (!((std::isnan(tt)) || (std::isinf(tt)))){
	p.tt = tt;
	p.pr = pr;
	fout.write(reinterpret_cast<char *>(&p), sizeof(p));
      }
    }
    fout.close();

    // // This writes as ascii (up to 6 significant digits?)
    // fileName = argv[i];

    // fileName.replace(fileName.find(".graph"), string(".graph").length(), ".pathdata");
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
