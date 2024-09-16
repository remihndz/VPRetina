import os
os.environ["OPENBLAS_NUM_THREADS"] = "5"

from typing import Dict
import networkx as nx
import numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.interpolate import interp1d
from scipy.spatial import Voronoi 
from random import shuffle
from tqdm import tqdm
from itertools import cycle

from DVC.dvc_utils import GenerateRandomPointOnDisk, FindClosestNodes, _findBifurcationPoints
from SVC.CCO.cco_utils import AVGraph

# To find intersection between lines
def ccw(A, B, C):
    return (C[:, 1] - A[:, 1]) * (B[:, 0] - A[:, 0]) > (B[:, 1] - A[:, 1]) * (C[:, 0] - A[:, 0])

def intersect(A, B, C, D):
    return np.logical_and(ccw(A, C, D) != ccw(B, C, D), ccw(A, B, C) != ccw(A, B, D))



def GenerateParameters(resultsFolder:str,
                       n:int,
                       baseTerm:list=[100,600,500],
                       **kwargs) -> dict:
    params = {**kwargs}
    params['ssmFileName'] = f"{resultsFolder}/SSM/sim_{n}"

    ## Patient parameters
    params['rCRA'] = kwargs.get('rCRA',
                                np.random.normal(163/2.0, 17/2.0))*1e-4 # Dorner et al. 2009
    params['vCRA'] = kwargs.get('vCRA',
                                np.random.normal(6.3, 1.2)) # Dorner et al. 2009
    params['MAP']  = kwargs.get('MAP',
                                np.random.normal(84, 6)) # Dorner et al. 2009
    params['IOP']  = kwargs.get('IOP',
                                np.random.normal(14.7, 2.8)) # Wang et al. 2018
    params['pCRV'] = kwargs.get('pCRV',
                                np.random.normal(15, 2.7)) # Stodmeister et al. 2013
    params['capPressure'] = kwargs.get('capPressure', 23) # Takahashi et al. 2009
    
    # File names
    params['confFileNameVein'] = f"{resultsFolder}/Coarse/tmp_vein.conf"
    params['confFileNameArtery'] = f"{resultsFolder}/Coarse/tmp_artery.conf"
    params['veinFileName'] = f"{resultsFolder}/Coarse/sim_{n}_vein"
    params['arteryFileName'] = f"{resultsFolder}/Coarse/sim_{n}_artery"
    params['AVTreeFileName'] = f"{resultsFolder}/AVTrees/sim_{n}_AV.cco"
    params['confFileNameMacula'] = f"{resultsFolder}/Macula/tmp_macula.conf"
    params['SVCFileName'] = f"{resultsFolder}/Macula/sim_{n}"

    # CCO parameters
    params['nTerm'] =  int(np.random.normal(baseTerm[0], 3))
    params['nTerm0'] = int(np.random.normal(baseTerm[1], baseTerm[1]/15.))
    params['nTerm1'] = int(np.random.normal(baseTerm[2], baseTerm[2]/15.))

    return params


class VirtualRetinalVasculature(AVGraph):
    '''
    Assumes the incoming graph has a nodeType attribute.
    The attribute is either 'art' or 'vei'.
    '''

    zICP = -125e-4
    zDCP = zICP-20e-4
    
    capillaryRadius = 2.5e-4 # Minimum capillary radius. Also the initial radius of the ICP/DCP vessels.
    plexusNames = {0:'SVP', 1:'ICP', 2:'DCP'}

    nodeDataToSave = ['position', 'plexus', 'nodeType', 'stage', 'pressure']
    edgeDataToSave = ['radius', 'length','hd', 'flow']

    _rDummy, _lDummy = 50e-4, 500e-4 # In cm, the parameters of the dummy vessels
    _RDummy = 1e6

    def __init__(self, ccoFile:str=None,
                 pCRA:float=None,
                 pCRV:float=None,
                 QCRA:float=None,
                 smoothing=False,
                 **kwargs) -> None:
        if not ccoFile:
            super().__init__()
            return
        
        super().__init__(ccoFile)
        self._isSVCConnected = False
        

        self._ccoFile = os.path.abspath(ccoFile) # From where the SVC tree has been read

        self.pCRA = pCRA # *133.322
        self.pCRV = pCRV # *133.322
        self.QCRA = QCRA

        # Connect the end vessels in the arterial tree to the venous tree
        assert nx.is_directed_acyclic_graph(self), "Not a DAG before adding capillaries."
        
        if kwargs.get('nPointsSVC', 3000):
            self._AddCapillaries(kwargs.get('ROI',0.2), kwargs.get('nPointsSVC', 3000), kwargs.get('alpha', 0.6))

        nx.set_node_attributes(self, values=0 , name='plexus')

        # self._DiameterSmoothing()

        assert nx.is_directed_acyclic_graph(self), "Not a DAG after adding SVC capillaries."
        if kwargs.get('nPointsICP', 1000) and kwargs.get('nPointsDCP', 1000):
            self._Connect_ICP_DCP(**kwargs)
        assert nx.is_directed_acyclic_graph(self), "Not a DAG after adding DVC capillaries."

        if smoothing:
            nx.relabel_nodes(self, {n:i for i,n in enumerate(self.nodes, start=self.number_of_nodes())}, copy=False)
            Smoothing(self, kwargs.get("alpha", 0.3), kwargs.get("nPointsSmoothing", 5))

        # self._DiameterSmoothing()

    @property
    def SVCVessels(self)->list:
        return [(n1,n2) for n1,n2 in self.edges
                if self.nodes[n1]['plexus']==self.nodes[n2]['plexus']==0]
    @property
    def ICPVessels(self)->list:
        return [(n1,n2) for n1,n2 in self.edges
                if self.nodes[n1]['plexus']==self.nodes[n2]['plexus']==1]
    @property
    def DCPVessels(self)->list:
        return [(n1,n2) for n1,n2 in self.edges
                if self.nodes[n1]['plexus']==self.nodes[n2]['plexus']==2]
    @property
    def ConnectingVessels(self)->list:
        return [(n1,n2) for n1,n2 in self.edges
                if self.nodes[n1]['plexus']!=self.nodes[n2]['plexus']]
    @property
    def _OutletsNotCRV(self)->list:
        return [n for n,d in self.out_degree if d==0 and not self.nodes[n].get('stage',1)==-2]
    @property
    def _InletsNotCRA(self)->list:
        return [n for n,d in self.in_degree if d==0 and not self.nodes[n].get('stage',1)==-3]
    @property
    def CRA(self)->int:
        inletNodes = [n for n,d in self.in_degree if d==0 and self.nodes[n].get('stage',1)==-3]
        assert len(inletNodes)==1, f"Found {len(inletNodes)} inlets."
        return inletNodes[0]
    @property
    def CRV(self)->int:
        outletNodes = [n for n,d in self.out_degree if d==0 and self.nodes[n].get('stage',1)==-2]
        assert len(outletNodes)==1, f"Found {len(outletNodes)} outlets."
        return outletNodes[0]
    
    def SplitFOVNotFOVVessels(self, FOV:float=0.4,
                              shape='square' # square or circle
                              )->tuple:
        """
        Return a list of vessels inside the FOV and a list
        of the vessels outside the FOV. 
        """
        fov, other = [], []
        fov2 = FOV/2.
        order = np.inf if shape=='square' else 2 # Inf norm if we want a squared FOV
        for n1,n2 in self.edges:
            if (np.linalg.norm(self.nodes[n1]['position'], order)<fov2
                and np.linalg.norm(self.nodes[n2]['position'], order)<fov2):
                fov.append((n1,n2))
            else:
                other.append((n1,n2))
        return fov, other
                
    def _flowInCompartment(self):
        try:
            return sum(f for u,v,f in self.edges.data('flow', None)
                       if v==-1 and self.nodes[u]['nodeType']=='art')
        except TypeError:
            raise TypeError("Flow has not been computed yet.")

    def _flowOutCompartment(self):
        try:
            return sum(f for u,v,f in self.edges.data('flow', None)
                       if u==-1 and self.nodes[v]['nodeType']=='vei')
        except TypeError:
            raise TypeError("Flow has not been computed yet.")

    def MaculaFlow(self):
        # return 1 - self._flowInCompartment()/self.TRBF()
        subG = self.subgraph(n for n,p in self.nodes.data('position') if sum(p[:2]**2)**.5<=.2
                             and self.nodes[n]['nodeType']=='art'
                             and self.nodes[n]['plexus']==0)
        return sum(f for u,v,f in subG.edges.data('flow') if subG.in_degree(u)==0)/self.TRBF()
    def TRBF(self):
        return max(f for u,v,f in self.edges.data('flow'))
    
    def _BifurcationLevel(self) -> dict:
        '''
        Computes the bifurcation level of all segments in the tree.
        This assumes the graph has not been connected

        ### Returns
        levels: dict[(n1,n2)]
            A dictionnary with edges as keys and their bifurcation levels as values.
        '''
        # if self._isSVCConnected:
        #     raise ValueError("Graph is already connected. Can't update bifurcation levels with this method.")

        levels = {}
        
        # Levels in the arterial tree
        paths = nx.all_simple_paths(self, self.CRA, [-1, self.CRV])
        
        for path in paths:
            level = 0
            n1 = path[0]
            for n2 in path[1:]:
                level = level + int(self.out_degree(n1)>1)
                levels[(n1,n2)] = level
                n1 = n2

        # Levels in the arterial tree
        reversedTree = nx.reverse(self, copy=True) # We need to reverse the tree to be able to do the same as above
        paths = nx.all_simple_paths(self, self.CRV, [-1])
        for path in paths:
            level = 0
            n1 = path[0]
            for n2 in path[1:]:
                level = level + int(reversedTree.out_degree(n1)>1)
                levels[(n2,n1)] = level # Reversed in our initial tree
                n1 = n2
        del reversedTree 

        return levels     
            
    
    def __str__(self):
        description = f'''Graph of retinal vasculature with {self.number_of_edges()} vessels.
        {len(self.SVCVessels)} SVC vessels
        {len(self.ICPVessels)} ICP vessels
        {len(self.DCPVessels)} DCP vessels
        pCRA: {getattr(self, "pCRA",np.nan)}
        pCRV: {getattr(self, "pCRV",np.nan)}
        {len(self.ConnectingVessels)} inter-plexi vessels
        ICP/DCP depth: {self.zICP}/{self.zDCP} [cm]
        Compartment resistance: {self._RDummy}
        SVC file: {self._ccoFile}
        Node attributes: {self.nodeDataToSave}
        Edge attributes: {self.edgeDataToSave}
        CRA: {self.CRA}
        CRV: {self.CRV}
        '''
        ## TODO: add some metrics e.g., vessel density
        return description

    @property
    def isSVCConnected(self):
        return _isSVCConnected

    def Write(self, fileName:str, mode:str='w'):
        '''
        Saves the retinal vascular network into a text file.

        ### Parameters
        fileName: str
            File to save the graph in.
        '''
        self._UpdateLength()
        with open(fileName, mode) as f:

            f.write(str(self) + '\n')
            # Node data
            f.write(f"\n# Nodes: {self.number_of_nodes()}\n")
            f.write(','.join(attr for attr in ["node"]+self.nodeDataToSave)+'\n')
            f.write('\n'.join(str(n)+','+','.join(str(d.get(attr, '')) for attr in self.nodeDataToSave) for n, d in self.nodes.data()))

            # Edge data
            f.write(f"\n# Edges: {self.number_of_edges()}\n")
            f.write(','.join(attr for attr in ["ProxNode", "DistNode"]+self.edgeDataToSave)+'\n')
                        
            f.write('\n'.join(
                str(n1)+','+str(n2)+','+
                ','.join(str(d.get(attr, None)) for attr in self.edgeDataToSave)
                    for n1,n2,d in self.edges.data()))
                
    def Read(self, fileName:str, smoothing:bool=False):
        '''
        Reads graph data from a file.

        ### Parameters
        fileName: str
            File to read the graph from.
        '''

        def ifFloat(x):
            try:
                return float(x)
            except ValueError:
                return x

        self.clear()
        self._ccoFile = os.path.realpath(fileName)
        
        with open(fileName, 'r') as f:
            line = f.readline()
            while not '# Nodes' in line:
                line = f.readline()
            nNodes = int(line.split(' ')[-1])
            attributes = f.readline().strip('\n').split(',')[1:] # First item is Name
            self.nodeDataToSave.extend([n for n in attributes[1:] if n not in self.nodeDataToSave])
            
            nodes = []
            for n in range(nNodes):
                node = f.readline().strip('\n').split(',')
                nodes.append((node[0], {attributes[i]:ifFloat(node[i+1]) if attributes[i]!='position'
                                        else np.array(list(map(float,node[i+1].strip('[').strip(']').split()))) for i in range(len(attributes))}))
            self.add_nodes_from(nodes)

            line = f.readline()
            while not '# Edges' in line:
                line = f.readline()
            nEdges = int(line.split(' ')[-1])
            attributes = f.readline().strip('\n').split(',')[2:] # First two items are names
            self.edgeDataToSave.extend([n for n in attributes[2:] if n not in self.edgeDataToSave])
            edges = []
            for n in range(nEdges):
                edge = f.readline().split(',')
                edges.append((edge[0], edge[1], {attributes[i]:ifFloat(edge[i+2]) for i in range(len(attributes))}))
            self.add_edges_from(edges)

        nx.relabel_nodes(self, {n:i for i,n in enumerate(self.nodes)}, copy=False)


    def _SplitVessels(self, vessels:list)->list:
        '''
        Split the specified vessels into two, adding an
        intermediate node in the middle.

        ### Parameters
        vessels: list
            A list of edges (n1,n2) or (n1,n2,d).

        ### Returns
        A list of the new nodes.
        '''

        nodeCount = sum((1 for n in self.nodes # Count the nodes to make sure we don't repeat names.
                         if (isinstance(n,str) and n.split('_')[0]=='intermediate') ))
        nodes = []
        for n1,n2 in vessels:

            intermediateNode = f"intermediate_{nodeCount}"
            p = (self.nodes[n1]['position']+self.nodes[n2]['position'])/2.0
            data = {**self.nodes[n1]}
            data.pop('pos', None);
            data['position'] = p
            self.add_node(intermediateNode, **{**data, 'plexus':0})

            data = self[n1][n2]
            self.remove_edge(n1,n2)
            self.add_edges_from([(n1,intermediateNode,data),
                                 (intermediateNode, n2, data)])

            nodes.append(intermediateNode)
            nodeCount+=1
        return nodes

    def _MakeSVCGraph(self,
                capillaryGraph:nx.Graph, 
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
        G = nx.DiGraph(capillaryGraph)
        #G.add_nodes_from(n for n in capillaryGraph.nodes(data=True))
        G.add_nodes_from(((f'art{i}', {'pos':arteries[i], 'nodeType':'art','original':i}) for i in range(arteries.shape[0])))
        G.add_nodes_from(((f'vei{i}', {'pos':veins[i], 'nodeType':'vei','original':i}) for i in range(veins.shape[0])))
        #G.add_edges_from(e for e in capillaryGraph.edges(data=True))

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
        allNodes = np.array([G.nodes[u]['pos'] for u in range(G.number_of_nodes())])
        dists = np.linalg.norm(np.array([G.nodes[u]['pos'] for u in danglingNodes])[:,np.newaxis,:] - allNodes, axis = 2)
        for i,n in enumerate(danglingNodes):
            u = np.argmin(dists[i, :n])
            #print((u,n), dists[i,u])
            edges.append((u,n))
        
        ## No children capillaries
        danglingNodes = [n for n,t in G.nodes.data('nodeType') if G.out_degree(n)==0 and t=='cap']
        dists = np.linalg.norm(np.array([G.nodes[u]['pos'] for u in danglingNodes])[:,np.newaxis,:] - allNodes, axis = 2)
        for i,n in enumerate(danglingNodes):
            u = np.argmin(dists[i, n+1:]) + n+1
            #print((n,u), dists[i,u])
            edges.append((n,u))

        G.add_edges_from(edges)

        # Add a 'position' node dictionary with (x,y,z) components (z constant, default 0)
        # TODO: delete 'pos' and make the notaiton and size of position vector consistent across modules
        # Posibility: overwrite the G.nodes[n]['pos'] to return a view of G.nodes[n]['position'] (showing only (x,y))
        # in order to avoid storing redundant information.
        G.add_nodes_from([(n,{'position':np.array([*p,0.0])}) for n,p in G.nodes.data('pos')])

        print("Dangling nodes left:", len([(n,t) for n,t in G.nodes.data('nodeType') if ((G.in_degree(n)==0 or G.out_degree(n)==0) and t=='cap')]))

        nx.set_edge_attributes(G, {e:{'radius':self.capillaryRadius if all(G.nodes[n]['nodeType']=='cap' for n in e) else 2*self.capillaryRadius} for e in G.edges})

        assert nx.is_directed_acyclic_graph(G), "Something went wrong and loops have been created. Assertion nx.is_directed_acyclic_graph(G) failed."
        return G

    def _MakeDVCGraph(self,
                     ICP:nx.Graph,
                     DCP:nx.Graph,
                     arteries:list,
                     veins:list,
                     ):
        """
        Update the graph to add the ICP and DCP.
        
        ### Parameters
        ICP, DCP: networkx.Graph
        Undirected graphs generated by the CapillaryBed function.
        arteries, veins: list
        Lists of nodes to be used to connect the capillary with a main graph.
        """

        nx.set_node_attributes(ICP, values='cap', name='nodeType')
        ICP.add_nodes_from([(n,{**self.nodes[n], 'nodeType':'art'}) for n in arteries])
        ICP.add_nodes_from([(n,{**self.nodes[n], 'nodeType':'vei'}) for n in veins])
        
        # Link arteries and veins to the closest capillaries in the ICP
        capillaries = [n for n,d in ICP.nodes.data() if d['nodeType']=='cap' and d['plexus']==1]
        mapping = {i:n for i,n in enumerate(capillaries)}
        xyCap = np.array([ICP.nodes[n]['position'] for n in capillaries])
        xyArteries = np.array([self.nodes[n]['position'] for n in arteries])
        xyVeins    = np.array([self.nodes[n]['position'] for n in veins])
        
        closestArt = np.argmin(np.linalg.norm(xyCap[:, :]-xyArteries[:,np.newaxis,:], axis=2), axis=1)
        closestVei = np.argmin(np.linalg.norm(xyCap[:,:]-xyVeins[:,np.newaxis,:], axis=2), axis=1)

        # We add an intermediate node that will serve as a connection points
        # with DCP (A1/V1 branch of the SVP-ICP connecting vessel, see An et al. 2020)
        intermediateArteries = [(f'intermediateICPArt_{i}', {**ICP.nodes[mapping[j]], 'plexus':1,
                                                             'position':(xyArteries[i]+xyCap[j])/2.0})
                                for i,j in enumerate(closestArt)] 
        intermediateVeins = [(f'intermediateICPVei_{i}', {**ICP.nodes[mapping[j]], 'plexus':1,
                                                          'position':(xyVeins[i]+xyCap[j])/2.0})
                             for i,j in enumerate(closestVei)]

        ICP.add_nodes_from(intermediateArteries+intermediateVeins)
        ICP.add_edges_from(((arteries[i], f'intermediateICPArt_{i}') for i in range(closestArt.size)))
        ICP.add_edges_from(((f'intermediateICPArt_{i}', mapping[closestArt[i]]) for i in range(closestArt.size)))

        ICP.add_edges_from(((mapping[closestVei[i]], f'intermediateICPVei_{i}') for i in range(closestVei.size)))
        ICP.add_edges_from(((f'intermediateICPVei_{i}', veins[i]) for i in range(closestVei.size)))        

        # Link arteries and veins to the closest capillaries in the DCP
        DCP.add_nodes_from(intermediateArteries+intermediateVeins)
        capillaries = [n for n,d in DCP.nodes.data() if d['nodeType']=='cap' and d['plexus']==2]
        mapping = {i:n for i,n in enumerate(capillaries)}
        xyCap = np.array([DCP.nodes[n]['position'] for n in capillaries])
        xyArt = np.array([n[1]['position'] for n in intermediateArteries])
        xyVei = np.array([n[1]['position'] for n in intermediateVeins])

        closestArt = np.argmin(np.linalg.norm(xyCap[:, :]-xyArteries[:,np.newaxis,:], axis=2), axis=1)
        closestVei = np.argmin(np.linalg.norm(xyCap[:,:]-xyVeins[:,np.newaxis,:], axis=2), axis=1)

        try:
            DCP.add_edges_from(((intermediateArteries[i][0], mapping[closestArt[i]]) for i in range(closestArt.size)))
            DCP.add_edges_from(((mapping[closestVei[i]], intermediateVeins[i][0]) for i in range(closestVei.size)))
        except KeyError:
            raise KeyError(f"{closestArt.shape=}, {closestVei.shape=}, {len(capillaries)=}, {DCP.number_of_nodes()=}, {xyArteries.shape=}")

        ## Combine ICP-DCP
        G = nx.DiGraph()
        G.add_nodes_from((n for n in ICP.nodes.data()))
        G.add_nodes_from((n for n in DCP.nodes.data()))
        G.add_edges_from((e for e in ICP.edges.data()))
        G.add_edges_from((e for e in DCP.edges.data()))

        ## Orient the graph
        undirectedCopyOfG = nx.Graph()
        undirectedCopyOfG.add_nodes_from([(n,{'nodeType':d}) for n,d in G.nodes(data='nodeType')])
        undirectedCopyOfG.add_edges_from(G.edges)
        
        L = nx.laplacian_matrix(undirectedCopyOfG)
        phi = np.zeros(L.shape[0])+0.5 # Node temperature
        sources = [i for i,(n,t) in enumerate(undirectedCopyOfG.nodes.data('nodeType')) if t=='art']
        sinks = [i for i,(n,t) in enumerate(undirectedCopyOfG.nodes.data('nodeType')) if t=='vei']
        
        for t in range(200): # Transient heat equation on graph
            phi[sources] = 1.0 
            phi[sinks] = 0.0 
            phi = phi-0.01*L.dot(phi)        
        topologicalOrder = {n:phi[i] for i,n in enumerate(undirectedCopyOfG.nodes)}
        del undirectedCopyOfG
        
        # Re-orient the edges according to topological order
        newEdges = [(n1,n2,d) if topologicalOrder[n1]>topologicalOrder[n2] else (n2,n1,d)
                    for n1,n2,d in G.edges(data=True)]
        G = G.to_directed()
        G.clear_edges()
        G.add_edges_from(newEdges)

        #TODO: vectorize this. Most time is spent on generating the 'candidates' array.
        # Fix dangling nodes
        edges = []
        # No-parent capillaries
        danglingNodes = [n for n,t in G.nodes.data('nodeType') if (t=='cap' and G.in_degree(n)==0)]        
        for n in danglingNodes:
            pos = G.nodes[n]['position']
            temp = topologicalOrder[n]
            plexus = G.nodes[n]['plexus']
            
            candidates = [u for u in G.nodes # Upstream nodes list
                          if (G.nodes[u]['nodeType']=='cap'
                              and G.nodes[u]['plexus']==plexus
                              and topologicalOrder[u]>temp)]
            candidatesPos = np.array([G.nodes[u]['position'] for u in candidates])
            
            u = candidates[np.argmin(np.linalg.norm(pos[np.newaxis,...]-candidatesPos[:,:], axis=1), axis=0)]
            edges.append((u,n))

        ## No-children capillaries
        danglingNodes = [n for n,t in G.nodes.data('nodeType') if (t=='cap' and G.out_degree(n)==0)]
        for n in danglingNodes:
            pos = G.nodes[n]['position']
            temp = topologicalOrder[n]
            plexus = G.nodes[n]['plexus']
            
            candidates = [u for u in G.nodes # Downstream nodes list
                          if (G.nodes[u]['nodeType']=='cap'
                              and G.nodes[u]['plexus']==plexus
                              and topologicalOrder[u]<temp)]
            candidatesPos = np.array([G.nodes[u]['position'] for u in candidates])
            u = candidates[np.argmin(np.linalg.norm(pos[np.newaxis,...]-candidatesPos[:,:], axis=1), axis=0)]
            edges.append((n,u))

        G.add_edges_from(edges)
        assert nx.is_directed_acyclic_graph(G), "Something went wrong and loops have been created. Assertion nx.is_directed_acyclic_graph(G) failed."

        # Add radius to the capillaries if they don't have one already
        nx.set_edge_attributes(G, {e:{'radius':self.capillaryRadius if all(G.nodes[n]['nodeType']=='cap' for n in e) else 2*self.capillaryRadius} for e in G.edges})

        self.add_nodes_from(n for n in G.nodes.data())
        self.add_edges_from(e for e in G.edges.data())

        # self._DiameterSmoothing()
        assert nx.is_directed_acyclic_graph(self), "Something went wrong and loops have been created. Assertion nx.is_directed_acyclic_graph(self) failed."
        nx.convert_node_labels_to_integers(self)
        dummyNode = [n for n,t in self.nodes.data('nodeType') if t=='dummy']
        assert len(dummyNode)<=1, f"More than one dummy nodes found. Found {dummyNode}"
        if dummyNode:
            nx.relabel_nodes(self, {dummyNode[0]:-1}, copy=False)
          
    def _Connect_ICP_DCP(self,
                         ROI:float=0.3, # In cm the FOV (circle)
                         arteryRadiusConditionSVC:callable=None,
                         veinRadiusConditionSVC:callable=None,
                         probaToKeepConnectingVessels:float=.3,
                         **kwargs
                         ):

        # SVC-ICP connections. We can introduce some condition to constrain the
        # number of connections or which vessels are allowed to connect the plexi
        # by defining these two functions.
        arteryCond = veinCond = lambda n1,n2,d:True
        # arteryCond = ((lambda n1,n2,d: 18e-4<2*d['radius']<4*22e-4
        #                #and np.linalg.norm([self.nodes[n1]['position'], self.nodes[n2]['position']], axis=1).min()<radiusROI
        #                )
        #               if not arteryRadiusConditionSVC
        #               else arteryRadiusConditionSVC)
        # veinCond = ((lambda n1,n2,d: 18.8e-4<2*d['radius']<4*22e-4
        #              #and np.linalg.norm([self.nodes[n1]['position'], self.nodes[n2]['position']], axis=1).min()<radiusROI
        #              )
        #             if not veinRadiusConditionSVC
        #             else veinRadiusConditionSVC)

        subgraph = self.subgraph(n for n,p in self.nodes.data('position') if np.linalg.norm(p[:2])<ROI/2)

        bifurcatingArterioles = self._SplitVessels([(n1,n2) for n1,n2,d in subgraph.edges(data=True)
                                                    if (subgraph.nodes[n1]['nodeType']==subgraph.nodes[n2]['nodeType']=='art')
                                                    and arteryCond(n1,n2,d)
                                                    and np.random.random()>probaToKeepConnectingVessels])
        bifurcatingVenules = self._SplitVessels([(n1,n2) for n1,n2,d in subgraph.edges(data=True)
                                                 if (subgraph.nodes[n1]['nodeType']==subgraph.nodes[n2]['nodeType']=='vei')
                                                 and veinCond(n1,n2,d)
                                                 and np.random.random()>probaToKeepConnectingVessels])

        print(f"{len(bifurcatingArterioles)=}, {len(bifurcatingVenules)=}")

        def _MakeCapBed(n, radiusFOV):
            vor, verticesToKeep, verticesToRemove = self._MakePoints(n, radiusFOV)

            edges = np.array([e for e in vor.ridge_vertices if all(n not in verticesToRemove for n in e)])
            newEdges = []
            vessels = np.array([[self.nodes[u]['position'][:2], self.nodes[v]['position'][:2]] for u,v in self.edges()])

            capBed = nx.Graph()
            capBed.add_nodes_from((i, {'pos':vor.vertices[i], 'nodeType':'cap'}) for i in edges.flat)
            capBed.add_edges_from((e[0],e[1], {'length':np.linalg.norm(np.diff(vor.vertices[e])), 'radius':self.capillaryRadius}) for e in edges)
            capBed = nx.convert_node_labels_to_integers(capBed)
            return capBed
        
        ICP = _MakeCapBed(n=kwargs.get('nPointsICP', 3000), radiusFOV=ROI/2)

        nx.relabel_nodes(ICP, {n:f'ICP{n}' for n in ICP.nodes}, copy=False)
        nx.set_node_attributes(ICP, values={n:np.append(p,self.zICP) for n,p in ICP.nodes.data('pos')}, name='position')
        nx.set_node_attributes(ICP, values=1, name='plexus')

        DCP = _MakeCapBed(n=kwargs.get('nPointsDCP', 3000), radiusFOV=ROI/2)

        
        # print('\n'.join(str(x) for x in np.linalg.norm(np.array([DCP.nodes[n]['pos'][:2] for n in DCP]), axis=1) if x>ROI))

        nx.relabel_nodes(DCP, {n:f'DCP{n}' for n in DCP.nodes}, copy=False)        
        nx.set_node_attributes(DCP, values={n:np.append(p,self.zDCP) for n,p in DCP.nodes.data('pos')}, name='position')
        nx.set_node_attributes(DCP, values=2, name='plexus')
        
        self._MakeDVCGraph(ICP, DCP, arteries=bifurcatingArterioles, veins=bifurcatingVenules)
        self._UpdateLength()

    def _MakePoints(self, n, radiusFOV):

        points = GenerateRandomPointOnDisk(n, radiusFOV) 
        vor = Voronoi(points)

        dists = np.min(np.linalg.norm(np.array([p[:2] for n,p in self.nodes.data('position')])[...] - vor.vertices[:,np.newaxis,:], axis=2), axis=1)
        mask = ((np.linalg.norm(vor.vertices, axis=1)>0.02) & (dists<1000*1e-4) & (np.linalg.norm(vor.vertices, axis=1)<radiusFOV*1.1))
        verticesToKeep = np.arange(vor.vertices.shape[0])[mask]
        verticesToRemove = np.append( np.arange(vor.vertices.shape[0])[~mask], -1)

        vertices = vor.vertices[verticesToKeep]
        return vor, verticesToKeep, verticesToRemove

    def _MakeSVCCapBed(self, n:int, radiusFOV,):
        
        vor, verticesToKeep, verticesToRemove = self._MakePoints(n, radiusFOV)

        edges = np.array([e for e in vor.ridge_vertices if all(n not in verticesToRemove for n in e)])
        newEdges = []
        vessels = np.array([[self.nodes[u]['position'][:2], self.nodes[v]['position'][:2]] for u,v in self.edges()])

        # Remove intersections
        for edge in edges:
            line = vor.vertices[edge]
            if np.any(intersect(line[0, np.newaxis], line[1,np.newaxis], vessels[:,0,:], vessels[:,1,:])):
                continue
            newEdges.append(edge)

        edges = np.array(newEdges)

        capBed = nx.Graph()
        capBed.add_nodes_from((i, {'pos':vor.vertices[i], 'nodeType':'cap'}) for i in edges.flat)
        capBed.add_edges_from((e[0],e[1], {'length':np.linalg.norm(np.diff(vor.vertices[e])), 'radius':self.capillaryRadius}) for e in edges)
        capBed = nx.convert_node_labels_to_integers(capBed)

        return capBed
                
    def _AddCapillaries(self,
                        ROI:float=0.3,
                        n:int=3000,
                        alpha:float=0.6) -> None:
        '''
        Add capillary connections between end arteriolar and
        end venous vessels. Within the disk of radius ROI, the
        vessels are linked to closest neighbour. 

        ### Parameters
        peripheralCapillariesTotalResistance: float
            The total resistance of the capillary compartment in the periphery (not ROI).
        ROI: float
            Region of interest, a disk within which we explicitly
            connect end vessels. The other vessels are connected
            via a dummmy node representing a 'capillary' compartment.
        n: int
            Number of seeds to use for the Voronoi diagram.
        alpha: float
            The fraction of venule/arteriole nodes that branch out to connect to the capillary bed.
        '''

        print(f"Adding capillary connections. Initial number of nodes/edges: {self.number_of_nodes()}/{self.number_of_edges()}")
        radiusFOV = ROI/2.
        macula = self.subgraph(n for n,p in self.nodes.data('position') if sum(p[:2]**2)<radiusFOV**2)

        # Create the initial capillary bed, not connected to self yet.  
        capBed = self._MakeSVCCapBed(n, radiusFOV) 

        # Define nodes to connect to the capillary bed
        random = np.random.random(size=macula.number_of_nodes()).tolist() # Used to randomly select nodes connecting to capillaries
        inletVeins = [n for n,d in macula.in_degree if d==0 and macula.nodes[n]['nodeType'] == 'vei']
        outletArteries = [n for n,d in macula.out_degree if d==0 and macula.nodes[n]['nodeType'] == 'art']
        arteries = [n for n,t in macula.nodes.data('nodeType') if (t=='art' and random.pop()<alpha)] + outletArteries
        veins = [n for n,t in macula.nodes.data('nodeType') if (t=='vei' and random.pop()<alpha)] + inletVeins
        
        # Make the connections
        capBed = self._MakeSVCGraph(capBed,
            np.array([macula.nodes[n]['position'][:2] for n in arteries]), 
            np.array([macula.nodes[n]['position'][:2] for n in veins]), 
            z=0
        )
        
        # Add the connections to self
        self.add_nodes_from((f'cap_{n}',d) for n,d in capBed.nodes.data() if d['nodeType']=='cap')
        self.add_edges_from(
            (arteries[capBed.nodes[u]['original']],f'cap_{v}', d) if (capBed.nodes[u]['nodeType']=='art' and capBed.nodes[v]['nodeType']=='cap')
            else (f'cap_{u}', veins[capBed.nodes[v]['original']], d) if (capBed.nodes[v]['nodeType']=='vei' and capBed.nodes[u]['nodeType']=='cap')
            else (f'cap_{u}', f'cap_{v}', d)
            for u,v,d in capBed.edges.data())

        self._UpdateLength()
        self._DiameterSmoothing()
        self._isSVCConnected = nx.is_directed_acyclic_graph(self) # True if it is acyclic
        
    def _UpdateLength(self, returnLengths=False):
        """
        Update the 'length' attribute of edges. 
        """
        lengths = np.linalg.norm(np.array([self.nodes[u]['position']-self.nodes[v]['position'] for u,v in self.edges]), axis=1)
        nx.set_edge_attributes(self, values={e:l for e,l in zip(self.edges, lengths) if not any(n==-1 for n in e)}, name='length')
        if returnLengths:
            return {e:l for e,l in zip(self.edges, lengths) if not any(n==-1 for n in e)}

    def _DiameterSmoothing(self, capillaryRadius:float=2.5e-4) -> None:
        '''
        Smooth the diameter transitions at bifurcations. 
        Before smoothing, all vessels with 
        radius <capillaryRadius are assigned capillaryRadius.

        ### Parameters
        capillaryRadius: float
            Minimal radius of a capillary. Default 2.5 microns 
            (assumes unit of radius is cm).
        '''
        CRA, CRV = self.CRA, self.CRV
        newRadii = {(n1,n2): np.mean([r if r>capillaryRadius else capillaryRadius
                                      for r in (self[n2][n]['radius'] for n in self.successors(n2))]
                                     +[r if r>capillaryRadius else capillaryRadius
                                      for r in (self[n][n1]['radius'] for n in self.predecessors(n1)) 
                                       ] + [self[n1][n2]['radius']])
                    for n1,n2 in self.edges if n1!=CRA and n2!=CRV}
        nx.set_edge_attributes(self, values=newRadii, name='radius')                
                
    def ComputeFlow(self) -> None:
        ''' 
        Compute blood flow in all vessels. 
        Pressures must be given in mmHg and flow in cm^3/s.

        ### Parameters
        inletBC, outletBC: dictionnary
            Specifiy inlet/outlet boundary condition as {'BCType':value}. 
            Valid boundary conditions types are 'flow' and 'pressure'.
        '''

        CRA,CRV = self.CRA, self.CRV # This checks we only have one CRA/CRV

        # Add connections to the peripheral vascular compartment
        self.remove_nodes_from([n for n,t in self.nodes.data('nodeType') if t=='dummy']) # Remove the compartment if existing
        originalHangingArteries, originalHangingVeins = self._OutletsNotCRV, self._InletsNotCRA 
        dummyParameters = {'dummy':True, 'resistance':self._RDummy, 'radius':self._rDummy, 'length':self._lDummy, 'hd':0.45}
        self.add_node(-1, nodeType='dummy', plexus=0, position=np.array(3*[np.inf]))        
        self.add_edges_from([(a,-1,dummyParameters) for a in originalHangingArteries])
        self.add_edges_from([(-1,v,dummyParameters) for v in originalHangingVeins])
        
        # Sanity check
        assert nx.is_directed_acyclic_graph(self), "Something went wrong and loops have been created. Assertion nx.is_directed_acyclic_graph(self) failed."

        nodeToInt = {n:i for i,n in enumerate(self.nodes)} ## Save the ordering. Allows for node names not being positive integers

        # Set up the boundary conditions
        ### Decision matrix and rhs            
        D = np.zeros((self.number_of_nodes(),))
        pBar = np.zeros_like(D)
        D[nodeToInt[CRA]] = 1.0
        pBar[nodeToInt[CRA]] = self.pCRA
        D[nodeToInt[CRV]] = 1.0
        pBar[nodeToInt[CRV]] = self.pCRV
 
        D = sp.dia_matrix(([D],[0]), shape = (self.number_of_nodes(), self.number_of_nodes()),
                            dtype=np.float32)
        I = sp.eye(D.shape[0], dtype=np.float32)
        R = sp.dia_matrix(([self.Resistances()], [0]), shape=(self.number_of_edges(), self.number_of_edges()))
        C = nx.incidence_matrix(self, oriented=True) * (-1.0)
        A = sp.vstack([sp.hstack([R, -C.T]),
                        sp.hstack([(I-D).dot(C), D])],
                        format='csr', dtype=np.float32)
        b = np.concatenate([np.zeros(self.number_of_edges()),
                            D.dot(pBar)])
        del R, D, I

        x = spl.spsolve(A,b,) # Find flow and pressure for each segment
        f, p = x[:self.number_of_edges()], x[self.number_of_edges():] # /133.322 # Edge flow, node pressures
        dp = C.T.dot(p) # Pressure

        nx.set_node_attributes(self, values={n:p[i] for n,i in nodeToInt.items()},
                               name='nodal pressure') # Add nodal pressure to graph data
        nx.set_edge_attributes(self, values={e:{'flow':f[i], 'pressure drop':dp[i]} for i,e in enumerate(self.edges())})

        # flowLoss = self[CRA][next(self.successors(CRA))]['flow'] - self[next(self.predecessors(CRV))][CRV]['flow']
        flowLoss = (sum(self[n][next(self.successors(n))]['flow'] for n,d in self.in_degree if d==0) -
                    sum(self[next(self.predecessors(n))][n]['flow'] for n,d in self.out_degree if d==0) )/ sum(self[n][next(self.successors(n))]['flow'] for n,d in self.in_degree if d==0)
        
        v = np.abs(np.array([self.GetVelocity(e) for e in self.edges() if self.nodes[e[0]]['nodeType']!='dummy' and self.nodes[e[1]]['nodeType']!='dummy']))
        print(f"Flow loss with max/min velocity={v.max()}/{v.min()}, max/min flow={abs(f).max()}/{abs(f).min()}, max/min pressure {p.max()}/{p.min()}:")
        print(f"\tFlow loss (flow in - flow out)/flow in = {flowLoss}")
        print(f"\tFlow loss (C*f).sum()/flow in = {(C@f).sum()/sum(self[n][next(self.successors(n))]['flow'] for n,d in self.in_degree if d==0)}")

    def GetVelocity(self, e:tuple) -> float:
        '''
        Given a tuple of nodes e, returns blood velocity in cm/s.
        If flow has not been calculated, returns -1.
        
        ### Parameters
        e: tuple
            Tuple of nodes forming an edge.
        '''
        try:
            return self[e[0]][e[1]]['flow']/((self[e[0]][e[1]]['radius']**2)*np.pi)
        except KeyError: # If flow has not been calculated yet
            return -1e12

    def Viscosity(self, radius, hd=0.45):
        ## In vitro law from Pries and Secomb
        d = 2*radius*1e4 # Convert to micron
        # mu045 = 220*np.exp(-1.3*d) + 3.2 - 2.44*np.exp(-0.06*d*0.645)
        C = (0.8 + np.exp(-0.075*d))*(-1+(1+(10**-11)*(d**12))**-1) + (1+(10**-11)*(d**12))**-1
        # visc = ( 1 + (mu045-1)*((1-hd)**C-1)/((1-0.45)**C-1) ) # In cP

        ## In vivo law from Pries and Secomb
        mu045 = 6*np.exp(-0.085*d) + 3.2 - 2.44*np.exp(-0.06*(d**0.645))
        visc = (1 + (mu045-1) * ((1-hd)**C-1)/(0.55**C-1) * ((d/(d-1.1))**2)) * ((d/(d-1.1))**2)
        
        #return 3.6 * 1e-3 / 133.332

        # delta = 4.29e-4 # Convert to cm
        # visc = 1.09 * np.exp(0.024 * hd)/(1+delta/radius)
        
        return visc * 1e-3/133.332 # Convert to mmHg.s

    def Resistances(self):
        '''
        Returns a vector of the resistances for each vessel.
        '''
        _res = self.Resistance
        mu = self.Viscosity
        resistance = lambda n1,n2,d: (self._RDummy if d.get('dummy', False) else _res((n1,n2)))

        R = [resistance(*e) for e in self.edges.data()]
        
        nx.set_edge_attributes(self, values={(n1,n2):mu(r) for n1,n2,r in self.edges.data('radius')}, name='viscosity')
        nx.set_edge_attributes(self, values={e:Ri for e,Ri in zip(self.edges, R)}, name='resistance')

        # print(f"Resistance ranges: {max(R)}-{min(R)}")
        # print(f"Viscosity ranges: {max(v for _,_,v in self.edges.data('viscosity'))}-{min(v for _,_,v in self.edges.data('viscosity'))}")

        return R

    def Resistance(self, e:tuple) -> float:
        '''
        Returns the resistance of vessel e=(n1,n2).
        
        ### Parameters
        e: tuple
            Tuple of nodes forming the edge.
        
        ### Returns
        The resistance of the given vessel
        '''
        mu = self.Viscosity # Might speed up by not having to look up the function         
        resistance = lambda d: 8*mu(d['radius'])*d['length']/(np.pi*(d['radius']**4))
        return resistance(self[e[0]][e[1]]) 


'''
A collection of functions to read a .cco into a
graph format, write a graph of vessels into .cco,
smooth a graph of vessels and plot a graph.
'''

import networkx as nx
import numpy as np
from statistics import StatisticsError, mean
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def PlotGraph(G:nx.DiGraph, alpha=None, nPoints:int=5, ax=None, **kwargs):
    '''
    Plot the graph edges.
    If alpha is given, the graph is first smoothed
    using Bezier curves.

    ### Parameters
    G: networkx.DiGraph
        The graph to plot. Nodes must have attribute 'position'.
    alpha: float
        The parameter controlling the tortuosity added by smoothing.
    nPoints: int
        Number of points to add during the smoothing process.
    ax: matplotlib.axis
        Axis on which to plot the graph. Default create a new figure.
    '''
    
    newG = G.copy(as_view=False)
    if alpha:
        Smoothing(newG, alpha, nPoints=nPoints)
                
    p = np.array([[newG.nodes[n1]['position'][:2], newG.nodes[n2]['position'][:2]] for n1,n2 in newG.edges])
    if not ax:
        fig, ax = plt.subplots()
        fig.suptitle(f'Plot of graph with smoothing parameter {alpha=}')
    ax.scatter(*p[:,1,:].T, s=0)
    lines = LineCollection(p, linewidths=[kwargs.get('radiusFactor',10)*r for n1,n2,r in newG.edges.data('radius', default=0.1)], color='black')
    ax.add_collection(lines)

def CreateGraph(ccoFile : str, convertUnitsTo) -> nx.DiGraph:

    lengthConversionDict = {'mm':1e1,
                            'cm':1.0,
                            'micron':1e4,
                            'mum':1e4,
                            'um':1e4,
                            'microns':1e4}
    try:
        lengthConversion = lengthConversionDict[str(convertUnitsTo)]
        print(f"{lengthConversion=}")
    except:
        raise KeyError(f"Wrong key '{convertUnitsTo}' for unit conversion. Valid keys are {lengthConversionDict.keys()}")

    G = nx.DiGraph()

    with open(ccoFile, 'r') as f:

        token = f.readline()
        token = f.readline().split() # Tree information

        # G.add_node(0, position=np.array([float(xi) for xi in token[:3]]) * lengthConversion)

        f.readline() # Blank line
        f.readline() # *Vessels
        nVessels = int(f.readline())
        print(f'The tree has {nVessels} vessels.')

        edges = dict()

        for i in range(nVessels):

            vessel = f.readline().split()
            vesselId = int(vessel[0])
            x1,x2 = np.array([float(x) for x in vessel[1:4]])*lengthConversion, np.array([float(x) for x in vessel[4:7]])*lengthConversion
            r = float(vessel[12]) * lengthConversion
            l = np.linalg.norm(x1-x2) * lengthConversion

            edges[vesselId] = {'radius':r, 'length':l,'start':x1,'end':x2,'stage':vessel[-1]}

            # G.add_node(vesselId+1, position=np.array([float(xi) for xi in vessel[4:7]]) * lengthConversion)

        f.readline() # Blank line
        f.readline() # *Connectivity

        rootId = None
        for i in range(nVessels):

            vessel = f.readline().split()
            vesselId = int(vessel[0])

            edges[vesselId]['parent'] = int(vessel[1])
            if int(vessel[1])==-1:
                rootId = vesselId
            edges[vesselId]['descendants'] = [int(descendant) for descendant in vessel[2:]]

            # if int(vessel[1]) == -1: # If root
            #     edges.append((0, 1, {'radius':vesselRadii[0], 'length':vesselLength[0], 'hd':0.45}))

            # else:
            #     vessel, parent = int(vessel[0]), int(vessel[1])
            #     edges.append((parent+1, vessel+1, {'radius':vesselRadii[vessel], 'length':vesselLength[vessel], 'hd':0.45}))

        # G.add_edges_from(edges)
        vesselId, node, nodep = rootId, 0, 1 # Start with the root


        def AddVesselToGraph(vesselId, startNode):
            endNode = G.number_of_nodes()
            vessel = edges.pop(vesselId)
            G.add_node(endNode, position=vessel['end'], stage=int(vessel.pop('stage')))
            G.add_edge(startNode, endNode, radius=vessel['radius'], length=vessel['length'], hd=vessel.pop('hd',0.45))

            for descendant in vessel['descendants']:
                AddVesselToGraph(descendant, endNode)

        G.add_node(0, position=edges[rootId]['start'], stage=-2)
        AddVesselToGraph(rootId, 0)                

        # nodeToRemove = [n for n,stage in G.nodes(data='stage')
        #                 if stage==-2 ]
        # for node in nodeToRemove:
        #     G.remove_node(node)

    return nx.convert_node_labels_to_integers(G)


def Graph2CCO(G:nx.DiGraph, ccoFileName:str, **kwargs):

    dp = kwargs.get('dp', 16.6*133.3)
    refPressure = kwargs.get('refPressure', 60*133.3)
    
    with open(ccoFileName, 'w') as f:

        nVessels = G.number_of_edges()
        root = [p for n,p in G.nodes.data('position') if G.in_degree(n)==0]
        d = [f for n1,n2,f in G.edges(data=True) if G.in_degree(n1)==0][0]
        assert len(root)==1, "More than one root found."
        root = root[0] # Unpack
        
        f.write("*Tree\n")
        treeInfo = f"{root[0]} {root[1]} {root[2]} {d.get('flow', 0.0)} {kwargs.get('psiFactor',9e-6)}"
        treeInfo+= f" {dp} 2 {refPressure} {G.number_of_nodes()} {d['radius']} {kwargs.get('tol', 1e-6)}\n\n"
        f.write(treeInfo)
        
        f.write("*Vessels\n")
        f.write(f"{nVessels}\n")
        branchingMode = 2
        edgeID = {(n1,n2):i for i,(n1,n2) in enumerate(G.edges)}
        vesselsInfo = [f"{edgeID[(n1,n2)]} {G.nodes[n1]['position'][0]} {G.nodes[n1]['position'][1]} {G.nodes[n1]['position'][2]}"
                       + f" {G.nodes[n2]['position'][0]} {G.nodes[n2]['position'][1]} {G.nodes[n2]['position'][2]}"
                       + f" 0.0 0.0 0.0 0.0 {branchingMode} {d['radius']} {d.get('flow',0.0)} 0.0 0.0 0.0 0 0.0 0.0 {d.get('stage', 0)}"
                       for i,(n1,n2,d) in enumerate(G.edges(data=True))]
        f.write('\n'.join(vesselsInfo))
        
        f.write('\n\n*Connectivity\n')
        vesselsConn = []
        for n1,n2 in G.edges:
            ID = edgeID[(n1,n2)]
            parent = list(G.predecessors(n1))
            if parent:
                parent = str(edgeID[(parent[0],n1)])
            else:
                parent = -1    
            children = ' '.join(str(edgeID[(n2,suc)]) for suc in G.successors(n2))
            vesselsConn.append(f"{ID} {parent} {children}")
            
        f.write('\n'.join(vesselsConn))

    return


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
    return [P0*((1-ti)**3) + C0*3*ti*((1-ti)**2) + 
            3*ti*ti*(1-ti)*C1 + ti*ti*ti*P1 
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

    CRA, CRV = G.CRA, G.CRV
    for n0,n1, d in G.edges(data=True):
        if n0 == CRA or n1 == CRV:
            continue
        try:
            newPts = BezierCurve(FindControlPoints(n0,n1), t)
            newNodes.extend([(nodeName+i, {**G.nodes[n0], 'position':p, 'Original node':False}) 
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

def RandomSmoothing(G:nx.DiGraph, alpha:float=0.5, nPoints:int=5):
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

        p0,p1 = G.nodes[n0]['position'], G.nodes[n1]['position']
        v = p1-p0
        v /= np.linalg.norm(v)
        normal = np.array([-v[1],v[0],0])
        k0,k1 = np.random.normal(0.25,0.02,2), np.random.normal(0.01, 0.002,2)
        if np.isnan(normal).any():
            return (p0, (p0+p1)/4, 3*(p0+p1)/4, p1)

        return (p0, p0 + k0[0]*v + np.random.choice([-1,1])*k1[0]*normal,
                p0 + (0.5+k0[1])*v + np.random.choice([-1,1])*k1[1]*normal, p1)
                
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

def SplitEdgesToLength(G:nx.DiGraph, l:float=None):
    '''
    Splits the edges of the graph until all have length <l.
    By default l is the length of the smallest edge.
    '''
    if not l:
        l = np.linalg.norm(np.array([G.nodes[n1]['position']-G.nodes[n2]['position'] for n1,n2 in G.edges]), axis=1).min()
        print(l, G)

    newPaths = []
    newNodes = []
    newNode = max([n for n in G.nodes])+1
    for n1,n2,d in G.edges(data=True):
               
        v = (G.nodes[n2]['position']-G.nodes[n1]['position'])
        dn1 = {key: value for key, value in G.nodes[n1].items() if key!='position'}
        p0 = G.nodes[n1]['position']
        li = np.linalg.norm(v)
        v/=li
        
        n = int(li/l)
        newDict = d
        dx = li/(2*n)
        newDict['length'] = dx
        
        newNodes.extend([(newNode+i, {'position':p0 + i*dx*v, **dn1}) for i in range(1,2*n)])
        newPaths.append([(n1,*[newNode+i for i in range(1,2*n)], n2), newDict])
        newNode = newNode+2*n-1

        # newNode = newPaths[-1][0][-2]+1
        # ns, path = SplitNTimes(G, (n1,n2), n, newNode)
        # newNode = ns[-1][0]+1
        # newNodes.extend(ns)
        # newPaths.append([path,newDict])
        
    G.clear_edges()
    G.add_nodes_from(newNodes)
    for path in newPaths:
        nx.add_path(G, path[0], **path[1])
        

def SplitNTimes(G:nx.DiGraph, e:tuple, n:int, newNodeName:int)->tuple:

    n1,n2 = e
    v = (G.nodes[n2]['position']-G.nodes[n1]['position'])
    dn1 = {key: value for key, value in G.nodes[n1].items()}
    _ = dn1.pop('position')
    p0 = G.nodes[n1]['position']
    li = np.linalg.norm(v)
    v/=li
    dx = li/(2**n)
    
    newNodes = [(newNodeName+i, {'position':p0+i*dx*v, **dn1}) for i in range(1,2*n)]
    newPath = [n1, *[newNodeName+i for i in range(1,2*n)], n2]
    return newNodes, newPath

def ComputeFlow(net:VirtualRetinalVasculature):
    """Computes the flow in the virtual retinal vasculature.
    This function can be adapted to more complex haemodynamics models.

    Args:
        net (VirtualRetinalVasculature): The vasculature. Should have the attributes pCRA, pCRV, _RDummy.
    """
    CRA, CRV = net.CRA, net.CRV
    
    # Add the connections to the peripheral vascular compartment
    terminalArteries, terminalVeins = self._OutletsNotCRV, self._InletsNotCRA 
    dummyParameters = {'dummy':True, 'resistance':net._RDummy, 'radius':net._rDummy, 'length':net._lDummy, 'hd':0.45}
    net.remove_nodes_from([n for n,t in self.nodes.data('nodeType') if t=='dummy']) # Remove if existing to avoid duplicates
    net.add_node(-1, nodeType='dummy', plexus=0, position=np.array(3*[np.inf]))     # Compartment is at [inf,inf,inf]        
    net.add_edges_from((a,-1,dummyParameters) for a in terminalArteries)   
    net.add_edges_from((-1,v,dummyParameters) for v in terminalVeins)
    assert nx.is_directed_acyclic_graph(self), "Something went wrong and loops have been created. Assertion nx.is_directed_acyclic_graph(self) failed."

    inlets = np.array([n for n,d in self.in_degree if d==0])
    outlets = np.array([n for n,d in self.out_degree if d==0])

    D = np.zeros(nn)
    pBar = np.zeros(nn)
    qBar = np.zeros(nn)
    for node, pressure in pressureConditions.items():
        D[node] = 1.0
        pBar[node] = pressure
    for node, flow in flowConditions.items():
        qBar[node] = flow
    
    D = sp.dia_matrix(([D],[0]), shape=(nn,nn))
    I = sp.eye(nn)
    R = sp.dia_matrix(([mu(d) for _,_,d in self.edges.data()], [0]), shape=(ne,ne))
    C = nx.incidence_matrix(self, oriented=True, nodelist=list(range(nn))) * (-1.0)
    A = sp.vstack([sp.hstack([R, -C.T]), 
                    sp.hstack([(I-D).dot(C), D])],
                    format='csc')
    rhs = np.concatenate([np.zeros((self.number_of_edges())),
                            (I-D).dot(qBar) + D.dot(pBar)])
    return A, rhs

