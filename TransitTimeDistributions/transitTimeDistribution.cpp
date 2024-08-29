#include "transitTimeDistribution.hpp"

// //graph[i][j] stores the j-th neighbour of the node i
pathAnalysis dfs(Graph &graph, size_t start, size_t end, size_t cutoff) 
{
  //initialize:
  pathAnalysis pathData;
  //remember the node (first) and the index of the next neighbour (second)
  std::stack<State> to_do_stack;
  std::deque<size_t> path; //remembering the way
  std::vector<bool> visited(graph.size(), false); //caching visited - no need for searching in the path-vector
  std::deque<float> pathTransitTime(cutoff), pathProbability(cutoff),
    pathOEF(cutoff); // Remembering the probability and transit time of the path. Used the same way as `path` but stores (cumulating) times/probabilities instead of nodes.

  // // Pre-allocate max size. Note: this keeps size() at 0.
  // path.reserve(cutoff);
  // pathTransitTime.reserve(cutoff);
  // pathProbability.reserve(cutoff);
  
  //start in start!
  to_do_stack.push(std::make_pair(start, 0));
  visited[start]=true;    
  path.push_back(start);
  pathTransitTime.push_back(0);
  pathProbability.push_back(1);
  pathOEF.push_back(1);

  size_t pathCount = 0, oneMil = static_cast<size_t>(1e6);
  std::cout << "\rNumber of paths found (in millions): " << pathCount / oneMil << std::flush;
  
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
	      pathData.pathsLength.push_back(path.size());
	      pathData.pathsOEF.push_back(1-pathOEF.back());
	      
	      if (pathCount % oneMil== 0){
		std::cout << "\rNumber of paths found (in millions): " << pathCount / oneMil << std::flush;
	      }
	      // print(path);//found a way!      
	    }
          //backtrack:
          visited[current.first]=false;//no longer considered visited
	  //go a step back	
          path.pop_back();
	  pathProbability.pop_back();
	  pathTransitTime.pop_back();
	  pathOEF.pop_back();
          to_do_stack.pop();//no need to explore further neighbours   
	}
      else{//normal case: explore neighbours
	size_t next=graph.vertices[current.first].adj[current.second];
	current.second++;//update the next neighbour in the stack!
	if(!visited[next]){
	  //putting the neighbour on the todo-list
	  to_do_stack.push(std::make_pair(next, 0));
	  visited[next]=true;
	  path.push_back(next);
	  pathTransitTime.push_back(pathTransitTime.back()+graph.vertices[current.first].tt[current.second-1]);
	  pathProbability.push_back(pathProbability.back()*graph.vertices[current.first].proba[current.second-1]);
	  pathOEF.push_back(pathOEF.back()*(std::exp(-K*graph.vertices[current.first].tt[current.second-1])));
	  // print(path);
	}      
      }
    }
  std::cout << std::endl;
  return pathData;
}

std::vector<size_t> getLeaves(const Graph &graph){
  std::vector<size_t> leaves;
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

std::vector<size_t> getRoots(const Graph &graph){
  std::vector<size_t> roots;
  std::vector<bool> isRoot(graph.vertices.size(), true);
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


// Printing helpers
void print(const std::vector<size_t> &path){
  for (auto i: path)
    std::cout << i << "->";
  std::cout << std::endl;
};

void print(const std::deque<size_t> &path){
  for (auto i: path)
    std::cout << i << "->";
  std::cout << std::endl;
};

void print(const std::deque<double> &path){
  for (auto i: path)
    std::cout << i << "->";
  std::cout << std::endl;
};

void print(const std::deque<float> &path){
  for (auto i: path)
    std::cout << i << "->";
  std::cout << std::endl;
};


void ReadGraph(std::ifstream& graphFile, Graph& graph)
{
  // The nodes in the text file might be strings. We will map them to consecutive integers
  std::map<std::string, size_t> nodesToInt;
  std::string line;
  size_t nv;
  std::string delimiter;
  std::string sCRA, sCRV; // The node names as strings
  
  while (getline(graphFile, line)){
    // Be careful that the reader doesn't read another
    // graph attribute that ends in "CRA:" or "CRV:"
    // by keeping a blank space before the string
    // to compare 'line' to, e.g., " CRA:"
    if (line.find(std::string(" CRA:")) != std::string::npos){
      delimiter = ": ";
      line.erase(0, line.find(delimiter)+delimiter.length());
      sCRA = line;
    }
    else if (line.find(std::string(" CRV:")) != std::string::npos){
      delimiter = ": ";
      line.erase(0, line.find(delimiter)+delimiter.length());
      sCRV = line;
    }
    else if (line.find(std::string("# Nodes")) != std::string::npos){ // If found the '# Nodes: nNodes' line
      delimiter = ": ";
      line.erase(0, line.find(delimiter)+delimiter.length()); // Removes the '# Nodes: ' before the number of nodes
      // std::cout << "Number of vertices " << line << std::endl;
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
    std::string node = line.substr(0, line.find(delimiter)); // Node name
    nodesToInt[node] = i;    
  }
  graph.CRA = nodesToInt[sCRA];
  graph.CRV = nodesToInt[sCRV];
  
  // Get the edges
  size_t ne;
  std::string edgeAttributeNames;
  int flowAttributePos=0, htAttributePos=0,
    radAttributePos=0, lenAttributePos=0; // That's how many comas we need to look for to find the column wiht flow values
  while (getline(graphFile, line)){
      if (line.find(std::string("# Edges")) != std::string::npos){
        delimiter = ": ";
	line.erase(0, line.find(delimiter)+delimiter.length()); // Removes the '# Edges: ' before the number of edges
	// std::cout << "Number of edges " << line << std::endl;
	ne = stoi(line);      
	getline(graphFile, edgeAttributeNames); // header line we could sparse
	std::string edgeAttribute;
	delimiter = ',';
	edgeAttributeNames.append(delimiter);
	int k =0;
	while(edgeAttributeNames.find(delimiter) != std::string::npos){
	  edgeAttribute = edgeAttributeNames.substr(0, edgeAttributeNames.find(delimiter));
	  //std::cout << "Edge attribute: " << edgeAttribute << std::endl;
	  edgeAttributeNames.erase(0, edgeAttributeNames.find(delimiter) + delimiter.length());
	  if (edgeAttribute.compare("flow")==0)
	    {
	      //std::cout << "Found flow columns at " << k << std::endl;
	      flowAttributePos = k;
	    }
	  else if (edgeAttribute.compare("hd")==0)
	    {
	      //std::cout << "Found ht columns at " << k << std::endl;
	      htAttributePos = k;
	    }
	  else if (edgeAttribute.compare("radius")==0)	    
	    {
	      //std::cout << "Found radius columns at " << k << std::endl;
	      radAttributePos = k;
	    }
	  else if (edgeAttribute.compare("length")==0){
	    // std::cout << "Found length columns at " << k << std::endl;
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
  std::string nodeName;
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
    while (line.find(delimiter) != std::string::npos){
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
	l = boost::lexical_cast<double>(line.substr(0, line.find(delimiter)));
	if (l>1e20)
	  l = CompartmentVesselsLength; // Either u or v is the compartment.
      }
      line.erase(0, line.find(delimiter)+delimiter.length());
      k++;
    }

    if ((u==graph.CRA || v==graph.CRV) && f<0)
      std::cout << "CRA or CRV orientation is being switched!" << std::endl;
    
    if (f>=0){
      graph.vertices[u].adj.push_back(v); // Add the edge
      graph.vertices[u].ht.push_back(0.45);
      graph.vertices[u].flow.push_back(f);
      graph.vertices[u].rad.push_back(r);
      graph.vertices[u].len.push_back(l);
      // std::cout << "Edge (" << u << "," << v << ") with flow " << graph[u].flow.back() << " and ht " << graph[u].ht.back() << std::endl;
    }
    else{ // Reverse the edge to follow flowx
      // graph.vertices[v].adj.push_back(u);
      graph.vertices[u].adj.push_back(v); // Add the edge
      graph.vertices[v].ht.push_back(0.45);
      graph.vertices[v].flow.push_back(-1.0*f);
      graph.vertices[v].rad.push_back(r);
      graph.vertices[v].len.push_back(l);
      // std::cout << "Original edge was reversed to become (" << v << "," << u << ") with flow " << graph.vertices[v].flow.back() << " and ht " << graph.vertices[v].ht.back() << std::endl;
    }
  }
};

