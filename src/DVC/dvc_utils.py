import networkx as nx
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.spatial import Delaunay, Voronoi
import matplotlib.pyplot as plt
from statistics import StatisticsError, mean

def CapillaryBed(nPoints:int=500,
                 tortuosityIndex:float=0.9,
                 capillaryRadius:float=2e-4,
                 radiusROI:float=0.3,
                 radiusFAZ:float=2.5e-2,
                 homogeneous:bool=True)->nx.DiGraph:
    '''
    Assumes length scale is cm.
    Creates a Voronoi diagram to be used as a capillary.
    
    ### Parameters
    - 
    '''
    ## Unused for now because it creates unrealistic behaviors (intersections of capillaries)
    alpha = 1-tortuosityIndex # Lower values of alpha creates straighter lines
        
    # nPoints = nTheta x nRadius
    # we split nTheta=2*nRadius/3
    nRadius = int((3.*nPoints/2.)**0.5)
    points = GetCentroidsOfDelaunayTriangulation(nTheta=int(2.*nRadius/3),
                                                 nR=nRadius,
                                                 r=radiusROI,
                                                 rMin=radiusFAZ,
                                                 homogeneous=homogeneous)
    G = MakeGraphCapillaryBed(points, radiusFAZ,
                              outerRadius=radiusROI,
                              capillaryRadius=capillaryRadius)
        
    return G

def GenerateRandomPointOnDisk(n:int, r:float, rMin:float=0) -> np.ndarray:
    """
    Generates random points on a disk or annulus.
    
    ### Parameters
    n: int
        The number of points to sample.
    r: float
        The outer radius of the disk.
    rMin: float
        The inner radius of the annulus. Default 0.
    
    ### Returns
    points: ndarray
        A (n,2) array of points.
    """
    X = (rMin + (r-rMin)*np.sqrt(np.random.random(n))) * np.cos(np.random.random(n)*2*np.pi)
    Y = (rMin + (r-rMin)*np.sqrt(np.random.random(n))) * np.cos(np.random.random(n)*2*np.pi)
    return np.vstack((X,Y)).T

def GetCentroidsOfDelaunayTriangulation(nTheta:int, nR:int,
                                        r:float, rMin:float,
                                        homogeneous:bool=True) -> np.ndarray:
    """
    Create a Delaunay triangulation of the annulus and returns
    its centroids.
    
    ### Parameters
    nTheta: int
        The number of \theta samples if homogeneous is True.
    nR: int
        The number of radial samples if homogeneous is True.
    r: float
        The outer radius of the annulus. 
    rMin: float
        The inner radius of the annulus.
    homogeneous: bool
        Whether or not to sample the points homogeneously. Default is True.

    ### Returns
    points: ndarray
        A (nTheta*nR, 2) array of points.
    """
    if homogeneous:
        u,v = np.meshgrid(np.linspace(0,2*np.pi, nTheta),
                        np.linspace(rMin, r, nR, endpoint=True))
        u,v = u.flatten(), v.flatten()
        points = np.vstack((v*np.cos(u), v*np.sin(u))).T
    else:
        points = GenerateRandomPointOnDisk(nTheta*nR, r, rMin)
    points = np.append(points, np.zeros((1,2)), axis=0)
    tri = Delaunay(points)
    return tri.points[tri.simplices].mean(1) # The centroids

def MakeGraphCapillaryBed(points:np.ndarray, radiusFAZ:float=2.5e-2,
                          outerRadius:float=0.4,
                          capillaryRadius:float=2e-4) -> nx.Graph:
    """
    Creates an undirected graph of capillaries
    from a set of points using a Voronoi tessellation.

    ### Parameters
    points: ndarray
        The points used a centroids for the Voronoi tessellation.
    radiusFAZ: float
        The radius of the inner circle which is void of vessels (FAZ).
    outerRadius: float
        The outer boundary of the capillary bed.
    capillaryRadius: float
        The radius of the vessels created vessels.
        
    ### Returns
    G: nx.Graph
        An undirected graph made of the vertices and edges of the Voronoi tessellation.
    """

    vor = Voronoi(points)
    G = nx.Graph()
    # Mask out the far away vertices
    maskEdges = ~(np.array(vor.ridge_vertices)==-1).any(axis=1)
    # Get rid of points in the FAZ or outside the boundaries we are interested in
    # TODO: Find a way to vectorize this.
    capillaries = []
    for e in np.array(vor.ridge_vertices)[maskEdges,:]:
        dists = np.linalg.norm(vor.vertices[e], axis=1)
        if ( (dists>outerRadius).any() or (dists<radiusFAZ).any() ):
            continue
        # plt.plot(*np.array(vor.vertices[e]).T,  linewidth=2,  c='k')
        capillaries.append(e)
    # plt.show()
    G.add_nodes_from((n, {'pos':vor.vertices[n], 'nodeType':'cap'}) for n in np.unique(capillaries))
    G.add_edges_from(capillaries, radius=capillaryRadius)
    nx.convert_node_labels_to_integers(G)
    print("Capillary bed:",G)
    return G
    

    
def MakeAVGraph(capillaryGraph:nx.Graph, 
                arteries:np.ndarray, 
                veins:np.ndarray,
                z:float=0
                ) -> nx.DiGraph:
    """
    Make an arterio-venous graph by linking the pre-capillary arterioles to
    the post-capillary venules through the capillary bed in G.
    This was used for development/testing. Use the methods vp_utils's
    VirtualRetinalVasculature instead.

    ### Parameters
    capillaryGraph: Graph
        The graph of the capillary bed.
    arteries: ndarray
        A list of vessel end points from an arterial tree to connect to the capillary bed.
    veins: ndarray
        A list of vessel end points from a venous tree to connect to the capillary bed.

    ### Returns
    G: DiGraph
        The graph of the arterio-venous tree.      
    """

    # We don't need DiGraph at this points, but for some reason, 
    # it gives weird behavior with a non-directed graph when adding the edges.
    G = nx.DiGraph()
    G.add_nodes_from(n for n in capillaryGraph.nodes(data=True))
    G.add_nodes_from(((f'art{i}', {'pos':arteries[i], 'nodeType':'art'}) for i in range(arteries.shape[0])))
    G.add_nodes_from(((f'vei{i}', {'pos':veins[i], 'nodeType':'vei'}) for i in range(veins.shape[0])))
    G.add_edges_from(e for e in capillaryGraph.edges(data=True))
    
    # Link the arteries and veins to closest capillary.
    xyCap = np.array([[n, *(capillaryGraph.nodes[n]['pos'])] for n,t in capillaryGraph.nodes.data('nodeType') if t=='cap'])
    closestArt = np.argmin(np.linalg.norm(xyCap[:, 1:]-arteries[:,np.newaxis,:], axis=2), axis=1)
    closestVei = np.argmin(np.linalg.norm(xyCap[:,1:]-veins[:,np.newaxis,:], axis=2), axis=1)
    G.add_edges_from(((f'art{i}', closestArt[i]) for i in range(closestArt.size)), radius=0.1)
    G.add_edges_from(((closestVei[i], f'vei{i}') for i in range(closestVei.size)), radius=0.1) 
    G = nx.convert_node_labels_to_integers(G, label_attribute='OldLabel')
    del xyCap, closestArt, closestVei

    # Find a topological ordering of the nodes from their distance to/from veins/arteries
    undirectedCopyOfG = nx.Graph()
    undirectedCopyOfG.add_nodes_from(G.nodes(data=True))
    undirectedCopyOfG.add_edges_from(G.edges(data=True))

    L = nx.laplacian_matrix(undirectedCopyOfG)
    phi = np.zeros(L.shape[0])+0.5 # Node temperature
    sources = [n for i,(n,t) in enumerate(undirectedCopyOfG.nodes.data('nodeType')) if t=='art']
    sinks = [n for i,(n,t) in enumerate(undirectedCopyOfG.nodes.data('nodeType')) if t=='vei']
    for t in range(2): # Transient heat equation on graph
        phi[sources] = 1.0 
        phi[sinks] = 0.0 
        phi = phi-0.1*L.dot(phi)
    topologicalOrder = np.flip(np.argsort(phi))
    topologicalOrder = sorted(topologicalOrder, key=lambda n:undirectedCopyOfG.nodes[n]['nodeType']) # Should put the arteries first, then capillaries then veins. Hopefully without changing the order
    del undirectedCopyOfG

    # Relabel according to that topological order
    G = nx.relabel_nodes(G, {topologicalOrder[i]:i for i in range(len(topologicalOrder))}, copy=True)
    edges = [(u,v) if u<v else (v,u) for u,v in G.edges]
    G.clear_edges()
    G.add_edges_from(edges)

    # Fix the dangling nodes
    edges = []
    ## No parents capillaries
    danglingNodes = [n for n,t in G.nodes.data('nodeType') if G.in_degree(n)==0 and t=='cap']
    for n in danglingNodes:
        pos = G.nodes[n]['pos']
        candidates = np.array([np.append(G.nodes[u]['pos'], u) for u in range(n) if G.nodes[u]['nodeType']!='vei'])
        u = int(candidates[np.argmin(np.linalg.norm(pos[np.newaxis,...]-candidates[:,:-1],axis=1), axis=0),-1])
        edges.append((u,n))
    ## No children capillaries
    danglingNodes = [n for n,t in G.nodes.data('nodeType') if G.out_degree(n)==0 and t=='cap']
    for n in danglingNodes:
        pos = G.nodes[n]['pos']
        candidates = np.array([np.append(G.nodes[u]['pos'], u) for u in range(n+1,G.number_of_nodes()) if G.nodes[u]['nodeType']!='vei']+
                                [np.append(data['pos'], u) for u,data in G.nodes.data() if data['nodeType']=='vei'])
        u = int(candidates[np.argmin(np.linalg.norm(pos[np.newaxis,...]-candidates[:,:-1],axis=1), axis=0),-1])
        edges.append((n,u))
    G.add_edges_from(edges)

    # Add a 'position' node dictionary with (x,y,z) components (z constant, default 0)
    # TODO: delete 'pos' and make the notaiton and size of position vector consistent across modules
    # Posibility: overwrite the G.nodes[n]['pos'] to return a view of G.nodes[n]['position'] (showing only (x,y))
    # in order to avoid storing redundant information.
    G.add_nodes_from([(n,{'position':np.array([*p,0.0])}) for n,p in G.nodes.data('pos')])
    
    assert nx.is_directed_acyclic_graph(G), "Something went wrong and loops have been created. Assertion nx.is_directed_acyclic_graph(G) failed."
    return G

def BezierCurve(points:list, t:"list[float]") -> list:
    """
    Return points from the cubic Bezier curve defined
    by the 4 points in points sampled at times t.

    ### Parameters
    points: list
        A list of 4 points (2D or 3D).
    t: float or np.ndarray
        The time(s) at which to sample new points.
    
    ### Returns
    newPoints: list
        A list of points on the Bezier curve.
    """
    assert len(points)==4, f"Cubic Bezier curves require 4 points. {len(points)} were given."
    P0, C0, C1, P1 = points
    return [P0*(1-ti)**3 + C0*3*ti*(1-ti)**2 + 
            3*ti*ti*(1-ti)*C1 + ti*ti*ti*P1*ti 
            for ti in t]

def Smoothing(G:nx.DiGraph, alpha:float=0.5, nPoints:int=5):
    """
    Smoothes the vessels by adding tortuosity (controlled by alpha)
    with Bezier curves. Method from Linninger et al. 2013.
    TODO: Diameter smoothing.

    ### Parameters
    G: nx.DiGraph
        The directed graph of the vasculature.
    alpha: float
        The level of tortuosity to add. 
        alpha=0 to keep straight lines, alpha=1 for maximum tortuosity.
    nPoints: int
        The number of points to sample the Bezier curve (regular intervals).
    """
    def FindControlPoints(n0, n1):
        A0 = [G.nodes[n]['position'] for n in G.predecessors(n0) if n!=n1]
        A0 = np.mean(A0, axis=0)
        A1 = [G.nodes[n]['position'] for n in G.successors(n1) if n!=n0]
        p0,p1 = G.nodes[n0]['position'], G.nodes[n1]['position']
        
        if not A1: # Do differently for the terminal vessels: it gives strange results
            #raise StatisticsError
            points = (p0, p0 + alpha * (p0-A0), p1 + alpha * (p1-np.random.normal(p1, [0.01,0.01,0])), p1)
            # t = [i/(nPoints+1) for i in range(0,nPoints+2)]
            # b = BezierCurve(points, t)
            # plt.plot(*np.array(b).T, linewidth=0.6)
            return points
        A1 = np.mean(A1, axis=0)
        
        if np.isnan(A1).any() or np.isnan(A0).any() or any([A0[-1]!=0.0, A1[-1]!=0.0]):
            raise StatisticsError
            #return (p0, 0.01*(p1+p0)/4, 3.99*(p1+p0)/4, p1)
        return  (p0, p0 + alpha * (p0-A0), p1 + alpha * (p1-A1), p1)

    nx.set_node_attributes(G, True, name='Original node')
    t = [i/(nPoints+1) for i in range(1,nPoints+1)]
    nodeName, newNodes, newPaths = sorted(list(G.nodes))[-1]+1, [], []

    for n0,n1, d in G.edges(data=True):
        try:
            newPts = BezierCurve(FindControlPoints(n0,n1), t)
            newNodes.extend([(nodeName+i, {'position':p, 'Original node':False}) 
                            for i,p in enumerate(newPts)])
            newPaths.append([(n0,*list(range(nodeName, nodeName+len(newPts))), n1), d])
            nodeName += len(newPts)
        except StatisticsError:
            newPaths.append(((n0,n1),d))
        
    G.clear_edges()
    G.add_nodes_from(newNodes)
    for newPath in newPaths:
        nx.add_path(G, newPath[0], **newPath[1]) # Add the new segments, sharing the original edge's data.
    return

def _findBifurcationPoints(G:nx.DiGraph,
                        arteryRadiusCondition:callable,
                        veinRadiusCondition:callable):
    """
    NOT USED. 

    Finds the bifurcation points in the graph for each type of vessel, where the radius of the parent vessel is greater than the
    child vessel.

    Parameters
    ----------
    G : networkx.DiGraph
        The graph of the SVC to find the bifurcation points in.
    arteryRadiusCondition : function
        A function that takes an edge tuple (n1,n2,data) as input and returns True 
        if it satisfies the condition for arteries, False otherwise.
    veinRadiusCondition : function
        A function that takes an edge tuple (n1,n2,data) as input and returns True
        if it satisfies the condition for veins, False otherwise.

    Returns
    -------
    tuple
        Two lists of edges, one for the arteries and one for the veins, which satisfy the conditions.
    """
    artery_edges = []
    vein_edges = []
    for n1,n2,d in G.edges.data(True):
        nodeType = G.nodes[n1]['nodeType'] # Vein or artery
        if (nodeType=='art' or nodeType=='artery') and arteryRadiusCondition(n1,n2,d):
            artery_edges.append((n1,n2))
        elif (nodeType=='vei' or nodeType=='vein') and veinRadiusCondition(n1,n2,d):
            vein_edges.append((n1,n2))
    return (artery_edges, vein_edges)

        
def FindClosestNodes(C:np.ndarray):
    '''
    Find the nodes assignment that minimize the edges distance.
    Note that, for the time being, it DOES allow for one node to be
    assigned to multiple others. 
    
    ### Parameters
    C: numpy array
        An n-by-m array of costs for each edge (here, distance between nodes).
    
    ### Returns
    ind_row: list
        A list of n nodes.
    ind_col: list
        A list of the closest node for each node.
    '''
    return list(range(C.shape[0])), list(np.argmin(C, axis=1))
