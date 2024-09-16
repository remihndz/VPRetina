import sys
sys.path.append("../src")
from vp_utils  import VirtualRetinalVasculature
import pandas as pd
import numpy as np 
from pathlib import Path 
import networkx as nx
import scipy.sparse as sp
from tqdm.auto import tqdm

outputGraphsDir = Path("/home/remi/VPRetina/TransitTimeDistributions/GraphsForExperimentOcclusion/")
outputGraphsDir.mkdir(exist_ok=True)

inputGraphsDir = Path('/media/Storage3.6TB/VariedPopulationOfRetinas/OCTAMetrics' )
# inputGraphsDir = Path('/home/remi/MATLAB_3DOxygen/RetinaWith90percentOccludedVein')
populationParameters = pd.concat((
    pd.read_csv('/media/Storage3.6TB/VariedPopulationOfRetinas/PopulationParameters.csv', index_col=0,
                usecols=["Unnamed: 0","rCRA", "vCRA"]),
    pd.read_csv('/media/Storage3.6TB/VariedPopulationOfRetinas/outputs.csv', index_col=0)), axis=1)

def ComputeFlow(G, qIn):
    '''
    Compute the flow using a target inflow to the retina (in cm^3/s).
    '''
    CRA, CRV = G.CRA, G.CRV

    # Add a compartment for hanging nodes (that way we only have one inlet, one outlet)
    G.remove_nodes_from([n for n,d in G.nodes.data('dummy') if d])
    hangingArteries, hangingVeins = G._OutletsNotCRV, G._InletsNotCRA
    if hangingArteries or hangingVeins:
        dummyParams = {'resistance':G._RDummy, 'radius':G._rDummy,
                       'length':np.inf, 'hd':0.45, 'dummy':True}
        G.add_node(-1, nodeType='dummy', plexus=0, position=np.array(3*[np.inf]))
        G.add_edges_from((a, -1, dummyParams) for a in hangingArteries)
        G.add_edges_from((-1, v, dummyParams) for v in hangingVeins)
     
    nv, ne = G.number_of_nodes(), G.number_of_edges()
    R = sp.dia_matrix((G.Resistances(), 0), shape=(ne,ne))

    C = nx.incidence_matrix(G, oriented=True) # nv x ne matrix
    C_T = C.T.copy()    
    D_p = sp.dia_matrix(([1 if (v==CRV) else 0 for v in G], [0]), shape=(nv,nv)) # Decision matrix
    C = (sp.eye(nv)-D_p)@C
    bcLocs = {'CRA' if v==CRA else 'CRV': i+ne for i,v in enumerate(G) if (v==CRA or v==CRV)} # Which row they correspond to
    
    A = sp.vstack([
    sp.hstack([R, -C_T]),
    sp.hstack([C, D_p])]).tocsr()
    b = np.zeros(nv+ne)
    b[bcLocs['CRA']] = -qIn
    b[bcLocs['CRV']] = 1 # "Ground" pressure
    sol = sp.linalg.spsolve(A, b)
    q = sol[:ne]/qIn # Flow relative to inlet flow
    for (_,_,d), qi in zip(G.edges.data(), q):
        d['flow'] = qi
        d['ht']   = 0.45
    # Delete the compartment before saving?
    # G.remove_node(-1)
    #G.Write(outputGraphsDir / Path(G._fileName).name)
        
graphFiles = list(inputGraphsDir.glob("*.graph"))[:20]
print(f"Found {len(graphFiles)} graph files.")

for graphFile in tqdm(graphFiles):
    simIndex = int(graphFile.name.split("_")[1]) # Simulation number
    G = VirtualRetinalVasculature()
    G.Read(graphFile)
    G._fileName = graphFile
    # rCRA, vCRA = populationParameters.loc[simIndex]
    # qIn = np.pi*rCRA*rCRA*vCRA
    qIn = 1.# populationParameters.loc[simIndex, 'TRBF']
    ComputeFlow(G, qIn)            
    G.Write("TransientTransitTimes/20NormalSimsWithNormalizedFlow/" + graphFile.name.replace("AV", "normalizedFlow"))
