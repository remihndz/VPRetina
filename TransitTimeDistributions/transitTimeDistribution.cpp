#include <iostream>
#include <vector>
#include <stack>
#include <utility>
#include <fstream>
#include <map>
#include <deque>

const double PI = 3.1415926535897932384626433832795028841971693993751058209;

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

enum State {Unvisited, Visited, Visiting};
struct Vertex{
  vector<size_t> adj; // Nodes it is connected to
  State state=Unvisited;
  vector<double> flow, rad, ht, len,
    tt, proba; // Edge properties
};

//typedef vector<Vertex> Graph;
struct Graph{
  vector<Vertex> adjacencyMatrix;
  size_t size(){return adjacencyMatrix.size();}
};

vector<size_t> getRoots(Graph &graph){
  vector<size_t> roots;
  vector<Vertex> &adj = graph.adjacencyMatrix;
  vector<bool> isRoot(adj.size(), true);
  for (size_t u=0; u<adj.size(); u++){
    for (auto v: adj[u].adj)
      isRoot[v] = false;
  }
  for (size_t u=0; u<adj.size(); u++)
    if (isRoot[u]){
      roots.push_back(u);
    }      
  return roots;
}

vector<size_t> getLeaves(Graph &graph){
  vector<size_t> leaves;
  size_t u = 0; // u will be the vertex index in adjacency matrix
  for (Vertex &v: graph.adjacencyMatrix){
    if (v.adj.empty()) // No neighbors
      leaves.push_back(u); 
    u++;
  }
  return leaves;
}

void GetEdgeData(Graph &graph){
  // Compute the probability of RBC at node u going to node v, where v is a downstream neighbor of u.
  vector<Vertex> &adj = graph.adjacencyMatrix;
    
  // Compute the probability of a RBC going from vertex u to v
  // and the transit time between u and v
  double totalRBC, rbc;
  int u=0;
  for (Vertex &v: adj){ // Iterate through vertices v
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
//   intMat &adj = graph.adjacencyMatrix; // Graph adjacency matrix
//   doubleMat &edgeProbabilities = graph.edgeProbabilities,
//     &edgeTransitTimes = graph.edgeTransitTimes;
//   //initialize:
//   pathAnalysis pathData;
//   vector<bool> visited(graph.size(), false);
// }

// //graph[i][j] stores the j-th neighbour of the node i
pathAnalysis dfs(Graph &graph, size_t start, size_t end, size_t cutoff) 
{
  vector<Vertex> &adj = graph.adjacencyMatrix; // Graph adjacency matrix shortcut
  //initialize:
  pathAnalysis pathData;
  //remember the node (first) and the index of the next neighbour (second)
  typedef pair<size_t, size_t> State;
  stack<State> to_do_stack;
  deque<int> path(cutoff); //remembering the way
  vector<bool> visited(adj.size(), false); //caching visited - no need for searching in the path-vector
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
  while(!to_do_stack.empty())
    {
      State &current = to_do_stack.top();//current stays on the stack for the time being...
      
      if (current.first == end || current.second == adj[current.first].adj.size() || path.size()>=cutoff)//goal reached or done with neighbours?
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
	size_t next=adj[current.first].adj[current.second];
	current.second++;//update the next neighbour in the stack!
	if(!visited[next]){
	  //putting the neighbour on the todo-list
	  to_do_stack.push(make_pair(next, 0));
	  visited[next]=true;
	  path.push_back(next);
	  pathTransitTime.push_back(adj[current.first].tt[current.second-1]);
	  pathProbability.push_back(adj[current.first].proba[current.second-1]);
	  // print(path);
	}      
      }
    }
  return pathData;
}

void ReadGraph(ifstream& graphFile, Graph& graph)
{
  // The nodes in the text file might be strings. We will map them to consecutive integers
  map<string, size_t> nodesToInt;
  string line;
  size_t nv;
  string delimiter;
  
  while (getline(graphFile, line)){
    if (line.find(string("# Nodes")) != string::npos){ // If found the '# Nodes: nNodes' line
      delimiter = ": ";
      line.erase(0, line.find(delimiter)+delimiter.length()); // Removes the '# Nodes: ' before the number of nodes
      cout << "Number of vertices " << line << endl;
      nv = stoi(line);
      getline(graphFile, line); // Header line we could sparse
      break;
    }
  }

  graph.adjacencyMatrix.reserve(nv); // Initialize with nv vertices
  
  // Make the map of nodes to integer
  delimiter = ",";
  for (int i=0; i<nv; i++){
    getline(graphFile, line);
    string node = line.substr(0, line.find(delimiter)); // Node name
    nodesToInt[node] = i;    
  }

  // Get the edges
  size_t ne;
  string edgeAttributeNames;
  int flowAttributePos=0, htAttributePos=0,
    radAttributePos=0, lenAttributePos=0; // That's how many comas we need to look for to find the column wiht flow values
  while (getline(graphFile, line)){
      if (line.find(string("# Edges")) != string::npos){
        delimiter = ": ";
	line.erase(0, line.find(delimiter)+delimiter.length()); // Removes the '# Edges: ' before the number of edges
	cout << "Number of edges " << line << endl;
	ne = stoi(line);      
	getline(graphFile, edgeAttributeNames); // header line we could sparse
	string edgeAttribute;
	delimiter = ',';
	int k =0;
	while(edgeAttributeNames.find(delimiter) != string::npos){
	  edgeAttribute = edgeAttributeNames.substr(0, edgeAttributeNames.find(delimiter));
	  edgeAttributeNames.erase(0, edgeAttributeNames.find(delimiter) + delimiter.length());
	  if (edgeAttribute.compare("flow"))
	    flowAttributePos = k;
	  else if (edgeAttribute.compare("hd"))
	    htAttributePos = k;
	  else if (edgeAttribute.compare("radius"))
	    radAttributePos = k;
	  else if (edgeAttribute.compare("length"))
	    lenAttributePos = k;
	  k++;
	}
	break;
      }
  }

  // Read edge data
  size_t u,v;
  float f, h, r, l;
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
    int k=1; // Starts at one because we remove two columns (u and v the node names)
    line.append(delimiter); // In case flow or hematocrit is the last column
    while (line.find(delimiter) != string::npos){
      if (k==flowAttributePos)
	  f = stof(line.substr(0, line.find(delimiter)));
      else if (k==htAttributePos)       
	  h = stof(line.substr(0, line.find(delimiter)));
      else if (k==radAttributePos)       
	  r = stof(line.substr(0, line.find(delimiter)));
      else if (k==lenAttributePos)       
	  l = stof(line.substr(0, line.find(delimiter)));
      line.erase(0, line.find(delimiter)+delimiter.length());
      k++;
    }
    
    if (f>=0){
      graph.adjacencyMatrix[u].adj.push_back(v); // Add the edge
      graph.adjacencyMatrix[u].ht.push_back(h);
      graph.adjacencyMatrix[u].flow.push_back(f);
      graph.adjacencyMatrix[u].rad.push_back(r);
      graph.adjacencyMatrix[u].len.push_back(l);
      // cout << "Edge (" << u << "," << v << ") with flow " << flow[u].back() << " and ht " << ht[u].back() << endl;
    }
    else{ // Reverse the edge to follow flowx
      graph.adjacencyMatrix[v].adj.push_back(u); 
      graph.adjacencyMatrix[v].ht.push_back(h);
      graph.adjacencyMatrix[v].flow.push_back(-1.0*f);
      graph.adjacencyMatrix[v].rad.push_back(r);
      graph.adjacencyMatrix[v].len.push_back(l);
      // cout << "Original edge was reversed to become (" << v << "," << u << ") with flow " << flow[v].back() << " and ht " << ht[v].back() << endl;
    }
  }
}

int main(void)
{

  Graph graph;

  // ifstream graphFile;
  // graphFile.open("./sim_58_AV.graph"); //"./sim_0_AV.graph");
  // ReadGraph(graphFile, graph);
  // cout << "Got graph" << endl;
  // graphFile.close();

  graph.adjacencyMatrix.resize(5);  
  graph.adjacencyMatrix[0].adj.push_back(1);
  graph.adjacencyMatrix[2].adj.push_back(1);
  graph.adjacencyMatrix[1].adj.push_back(3);
  graph.adjacencyMatrix[1].adj.push_back(4);

  int ne=0, nv=0;
  for (Vertex& v: graph.adjacencyMatrix){
    for (size_t k=0; k<v.adj.size(); k++){
      v.ht.push_back(0.45);
      v.flow.push_back(1.0);
      v.rad.push_back(1.0);
      v.len.push_back(1.0);
      v.proba.push_back(1.0);
      v.tt.push_back(1.0);
      ne+=1;
    }
    nv++;
  }
      
  vector<size_t> roots = getRoots(graph),
    leaves = getLeaves(graph);
  
  print(roots);
  print(leaves);
  GetEdgeData(graph);
  cout << "Got edge data" << endl;
  auto pathsData = dfs(graph, roots[0], leaves[0], 300);
  cout << "Found " << pathsData.numberOfPaths() << " paths." << endl;
  
  // // vector<vector<size_t> > graph(5);

  // // graph[0].push_back(1);
  // // graph[2].push_back(1);
  // // graph[1].push_back(3);
  // // graph[1].push_back(4);

  // vector<vector<size_t>> graph(4);
  // for (int i=0; i<graph.size(); i++){
  //   for (int j=0; j<graph.size(); j++)
  //     if (i!=j)
  // 	graph[i].push_back(j);
  // }
  
  // vector<vector<double>> edgeProbas(graph.size());
  // for (int i=0; i<graph.size(); i++){
  //   for (int j=0; j<graph[i].size(); j++)
  //     edgeProbas[i].push_back(1.0);
  // }
  
  // auto pathData = dfs(fgraph, 0, 2, edgeProbas, edgeProbas, 10);

  // vector<size_t> roots = {0}, leaves = {3,2};
  // for (auto root: roots)
  //   for (auto leaf: leaves)
  //     dfs(graph, root, leaf, edgeProbas, edgeProbas, 10);
  return 0;
}
