#include <iostream>
#include <vector>
#include <stack>
#include <utility>
#include <fstream>
#include <map>

// TODO: give path vectors a fixed size (cutoff) and keep track of depth to save some time on checks from pop_back

using namespace std;
struct pathAnalysis {
  vector<double> pathsTransitTimes;
  vector<double> pathsProbabilities;
};

struct edgeData {
  vector<double> edgeTransitTimes;
  vector<double> edgeProbabilities;
}

void print(const vector<size_t> &path){
  for (auto i: path)
    cout << i << "->";
  cout << endl;
}

  
//graph[i][j] stores the j-th neighbour of the node i
pathData dfs(const vector<vector<size_t> > &graph, size_t start, size_t end, const vector<vector<double>> &edgeProbabilities, const vector<vector<double>> &edgeTransitTimes, size_t cutoff) 
{
  //initialize:
  pathAnalysis pathData; 
  //remember the node (first) and the index of the next neighbour (second)
  typedef pair<size_t, size_t> State;
  stack<State> to_do_stack;
  vector<size_t> path; //remembering the way
  vector<bool> visited(graph.size(), false); //caching visited - no need for searching in the path-vector
  vector<double> pathTransitTime, pathProbability; // Remembering the probability and transit time of the path. Use the same way as `path` but stores (cumulating) times/probabilities instead of nodes.
  
  //start in start!
  to_do_stack.push(make_pair(start, 0));
  visited[start]=true;    
  path.push_back(start);
  pathTransitTime.push_back(0);
  pathProbability.push_back(1);
  
  while(!to_do_stack.empty())
    {
      State &current = to_do_stack.top();//current stays on the stack for the time being...
      
      if (current.first == end || current.second == graph[current.first].size() || path.size()>=cutoff)//goal reached or done with neighbours?
	{
          if (current.first == end)
	    {
	      pathData.pathsProbabilities.push_back(pathProbability.back());
	      pathData.pathsTransitTimes.push_back(pathTransitTime.back());
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
	size_t next=graph[current.first][current.second];
	current.second++;//update the next neighbour in the stack!
	if(!visited[next]){
	  //putting the neighbour on the todo-list
	  to_do_stack.push(make_pair(next, 0));
	  visited[next]=true;
	  path.push_back(next);
	  pathTransitTime.push_back(edgeTransitTimes[current.first][next]);
	  pathProbability.push_back(edgeProbabilities[current.first][next]);
	}      
      }
    }
}

void ReadGraph(ifstream& graphFile, vector<vector<size_t>>& graph,
	       vector<vector<double>>& flow, vector<vector<double>>& ht)
{
  // The nodes in the text file might be strings. We will map them to consecutive integers
  map<string, size_t> nodesToInt;
  string line;
  bool readNodes = false;
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
  // Make the map of nodes to integer
  delimiter = ",";
  for (int i=0; i<nv; i++){
    getline(graphFile, line);
    string node = line.substr(0, line.find(delimiter)); // Node name
    nodesToInt[node] = i;
    graph.push_back(vector<size_t>());
    flow.push_back(vector<double>());
    ht.push_back(vector<double>());
  }

  // Get the edges
  size_t ne;
  string edgeAttributeNames;
  int flowAttributePos=0, htAttributePos=0; // That's how many comas we need to look for to find the column wiht flow values
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
	  if (edgeAttribute.compare("hd"))
	    htAttributePos = k;
	  k++;
	}
	break;
      }
  }
  // Read edge data
  size_t u,v;
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
    graph[u].push_back(v); // Add the edge
    // Find the column that has flow value for the vessel
    int k=1; // Starts at one because we remove two columns (u and v the node names)
    line.append(delimiter); // In case flow or hematocrit is the last column
    while (line.find(delimiter) != string::npos){
      if (k==flowAttributePos)
	{
	  string f = line.substr(0, line.find(delimiter));
	  flow[u].push_back(stof(f));
	}
      else if (k==htAttributePos)
	{
	  string h = line.substr(0, line.find(delimiter));
	  ht[u].push_back(stof(h));
	}
      line.erase(0, line.find(delimiter)+delimiter.length());
      k++;
    }
    // cout << "Edge (" << u << "," << v << ") with flow " << flow[u].back() << " and ht " << ht[u].back() << endl;
  }
}
   	     
int main(void)
{

  vector<vector<size_t>> graph;
  vector<vector<double>> flow, ht;

  ifstream graphFile;
  graphFile.open("./sim_0_AV.graph");
  ReadGraph(graphFile, graph, flow, ht);
  graphFile.close();

  vector<bool> isRoot(graph.size(), true);
  vector<size_t> leaves, roots;
  for (int u=0; u<graph.size(); u++){
    if(graph[u].empty())
      leaves.push_back(u);
    for (auto v: graph[u])
      isRoot[v] = false;
  }
  for (int i=0; i<graph.size(); i++)
    if (isRoot[i])
      roots.push_back(i);

  auto edgeData = getEdgeProbabilitiesAndTransitTimes(graph, flow, ht);

  auto pathData = dfs(graph, roots[0], leaves[0], edgeData.edgeProbabilities, edgeData.edgeTransitTimes);
  
  
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
