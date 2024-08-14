#include <iostream>
#include <vector>
#include <stack>
#include <utility>
#include <fstream>
#include <map>

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

struct Graph{
  intMat adjacencyMatrix;
  doubleMat flow, radius, ht, length,
    edgeProbabilities, edgeTransitTimes;
  size_t size(){return adjacencyMatrix.size();}
};

class pathVectorInt {
public:
  size_t size;
  vector<size_t> path;
  pathVectorInt(size_t maxSize){
    size = 0;
    path.reserve(maxSize);
  }
  size_t pop_back(){
    size_t val = path[size]; // The last value
    size -= 1; // We removed the last value so size = size-1
    return val;
  }
  void push_back(size_t val){
    path[size] = val; // Write the value in our vector
    size+=1; // Move the cursor to the back of the vector
  }
};

class pathVectorDouble {
public:
  size_t size;
  vector<double> path;
  pathVectorDouble(size_t maxSize){
    size = 0;
    path.reserve(maxSize);
  }
  double pop_back(){
    size_t val = path[size]; // The last value
    size -= 1; // We removed the last value so size = size-1
    return val;
  }
  void push_back(double val){
    path[size] = val; // Write the value in our vector
    size+=1; // Move the cursor to the back of the vector
  }
};


vector<size_t> getRoots(Graph &graph){
  vector<size_t> roots;
  intMat &adj = graph.adjacencyMatrix;
  vector<bool> isRoot(adj.size(), true);
  for (size_t u=0; u<adj.size(); u++){
    for (auto v: adj[u])
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
  int u = 0;
  for (auto vec: graph.adjacencyMatrix){
    if (vec.empty())
      leaves.push_back(u);
    u++;
  }
  return leaves;
}

void GetEdgeData(Graph &graph){
  // Compute the probability of RBC at node u going to node v, where v is a downstream neighbor of u.
  doubleMat &tt = graph.edgeTransitTimes,
    &pr   = graph.edgeProbabilities,
    &flow = graph.flow,
    &ht   = graph.ht,
    &rad  = graph.radius,
    &len  = graph.length; // References to the graph's data vectors for cleaner code
  intMat &adj = graph.adjacencyMatrix;

  // Allocate memory and size
  pr.resize(graph.size());
  tt.resize(graph.size());
  for (size_t u=0; u<graph.size(); u++){
    tt[u].resize(flow[u].size());
    pr[u].resize(flow[u].size());
  }
  
  // Compute the probability of a RBC going from vertex u to v
  // and the transit time between u and v
  double totalRBC, rbc;
  for (size_t u=0; u<adj.size(); u++){ // Iterate through vertex u
    totalRBC = 0;
    for (size_t k=0; k<adj[u].size(); k++){ // Iterate through neighbors number k
      rbc = flow[u][k]*ht[u][k]; // Flow of RBC from u to v
      tt[u][k] = PI*rad[u][k]*rad[u][k]*len[u][k]/flow[u][k];
      pr[u][k] = rbc;
      totalRBC += rbc;
    }
    for (size_t k=0; k<pr[u].size(); k++)
      pr[u][k] /= totalRBC; // Scale to a probability
  }      
}
  
//graph[i][j] stores the j-th neighbour of the node i
pathAnalysis dfs(Graph &graph, size_t start, size_t end, size_t cutoff) 
{
  intMat &adj = graph.adjacencyMatrix; // Graph adjacency matrix
  doubleMat &edgeProbabilities = graph.edgeProbabilities,
    &edgeTransitTimes = graph.edgeTransitTimes;
  //initialize:
  pathAnalysis pathData;
  //remember the node (first) and the index of the next neighbour (second)
  typedef pair<size_t, size_t> State;
  stack<State> to_do_stack;
  pathVectorInt path(cutoff); //remembering the way
  vector<bool> visited(adj.size(), false); //caching visited - no need for searching in the path-vector
  pathVectorDouble pathTransitTime(cutoff), pathProbability(cutoff); // Remembering the probability and transit time of the path. Use the same way as `path` but stores (cumulating) times/probabilities instead of nodes.

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
      
      if (current.first == end || current.second == adj[current.first].size() || path.size>=cutoff)//goal reached or done with neighbours?
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
	size_t next=adj[current.first][current.second];
	current.second++;//update the next neighbour in the stack!
	if(!visited[next]){
	  //putting the neighbour on the todo-list
	  to_do_stack.push(make_pair(next, 0));
	  visited[next]=true;
	  path.push_back(next);
	  pathTransitTime.push_back(edgeTransitTimes[current.first][current.second-1]);
	  pathProbability.push_back(edgeProbabilities[current.first][current.second-1]);
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
  bool readNodes = false;
  size_t nv;
  string delimiter;

  doubleMat &flow = graph.flow,
    &ht   = graph.ht,
    &rad  = graph.radius,
    &len  = graph.length; // References to the graph's data vectors for cleaner code
  intMat &adj = graph.adjacencyMatrix;
  
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
  // Make the map of nodes to integer
  delimiter = ",";
  for (int i=0; i<nv; i++){
    getline(graphFile, line);
    string node = line.substr(0, line.find(delimiter)); // Node name
    nodesToInt[node] = i;    
  }
  // Create the arrays
  adj.resize(nv);
  flow.resize(nv);
  ht.resize(nv);
  rad.resize(nv);
  len.resize(nv);

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
      adj[u].push_back(v); // Add the edge
      ht[u].push_back(h);
      flow[u].push_back(f);
      rad[u].push_back(r);
      len[u].push_back(l);
      // cout << "Edge (" << u << "," << v << ") with flow " << flow[u].back() << " and ht " << ht[u].back() << endl;
    }
    else{ // Reverse the edge to follow flowx
      adj[v].push_back(u); 
      ht[v].push_back(h);
      flow[v].push_back(-1.0*f);
      rad[v].push_back(r);
      len[v].push_back(l);
      // cout << "Original edge was reversed to become (" << v << "," << u << ") with flow " << flow[v].back() << " and ht " << ht[v].back() << endl;
    }
  }
}

int main(void)
{

  Graph graph;

  ifstream graphFile;
  graphFile.open("./sim_58_AV.graph"); //"./sim_0_AV.graph");
  ReadGraph(graphFile, graph);
  cout << "Got graph" << endl;
  graphFile.close();

  // graph.adjacencyMatrix.resize(5);
  // graph.flow.resize(5);
  // graph.radius.resize(5);
  // graph.length.resize(5);
  // graph.ht.resize(5);
  // graph.edgeProbabilities.resize(5);
  // graph.edgeTransitTimes.resize(5);
  
  // graph.adjacencyMatrix[0].push_back(1);
  // graph.adjacencyMatrix[2].push_back(1);
  // graph.adjacencyMatrix[1].push_back(3);
  // graph.adjacencyMatrix[1].push_back(4);

  // for (size_t u=0; u<graph.size(); u++){
  //   int k = 0;
  //   for (auto v: graph.adjacencyMatrix[u]){
  //     graph.flow[u].push_back(1.0);
  //     graph.radius[u].push_back(1.0);
  //     graph.length[u].push_back(1.0);
  //     graph.ht[u].push_back(1.0);
  //     graph.edgeProbabilities[u].push_back(1.0);
  //     graph.edgeTransitTimes[u].push_back(1.0);
  //   }
  // }
    
  vector<size_t> roots = getRoots(graph),
    leaves = getLeaves(graph);
  
  // print(roots);
  //print(leaves);
  GetEdgeData(graph);
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
