import sys, os
import numpy as np
import pandas as pd
import vp_utils as vp
import networkx as nx
from tqdm import tqdm


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


DVC = 'on'

csvOutFileFolder = sys.argv[1][:-4]+ '-Haemodynamics'
# csvOutFileFolder = os.path.join(*sys.argv[1].split('/')[:-1], "Haemodynamics-DVCOn" if DVC=='on' else "Haemodynamics-DVCOff")

print("Output folder:", csvOutFileFolder)
os.system("mkdir -p " + csvOutFileFolder)
population = pd.read_csv(sys.argv[1])

for i,row in tqdm(list(population.iterrows())):
    
    MAP, IOP, CRVP = row['MAP'], row['IOP'], row['CRVP']
    QCRA = row['vCRA']*np.pi*(row['rCRA']**2)
    pCRA = 2./3.*MAP

    # Skip if we have already simulated haemodynamics for this graph. 
    # if os.path.isfile(os.path.join(csvOutFileFolder, row['AVTreeFileName'].split('/')[-1][:-4]+'.graph')): 
    #     print("Skipping file:", row['AVTreeFileName'])
    #     continue

    try:
        with HiddenPrints():
            G = vp.VirtualRetinalVasculature(row['AVTreeFileName'], nPointsICP = 8000 if DVC=='on' else 0 ,
                                             nPointsDCP = 8000 if DVC=='on' else 0,
                                             nPointsSVC = 3000,
                                             pCRA=pCRA, pCRV=CRVP, QCRA=QCRA)
            G._RDummy = 1e6
            G.Write(os.path.join(csvOutFileFolder, row['AVTreeFileName'].split('/')[-1][:-4]+'.graph'))

    except FileNotFoundError:
        print("Did not find file:", row['AVTreeFileName'])
        continue


    #G.remove_nodes_from([n for n,p in G.nodes.data('plexus') if p!=0])
    CRA, CRV = G.CRA, G.CRV

    with HiddenPrints():
        G.ComputeFlow({'pressure':G.pCRA}, {'pressure':G.pCRV})

    # Writing the results
    inflow = lambda: sum(G[CRA][u]['flow'] for u in iter(G[CRA]))
    maculaFlow = G.MaculaFlow

    fileName = csvOutFileFolder + '/' + row['SVCFileName'].split('/')[-1] + f"_{100*maculaFlow():.0f}.csv"
    edgeData = ['radius', 'length', 'viscosity', 'hd', 'flow', 'pressure drop']
    print(fileName)    

    with open(fileName, 'w') as f:

        f.write(','.join(['idx']+edgeData+['inlet pressure',
                                           'outlet pressure',
                                           'inlet node',
                                           'outlet node',
                                           'inlet node type',
                                           'outlet node type',
                                           'inlet node plexus',
                                           'outlet node plexus',
                                           'isMacula'])) # Header
        f.write('\n')
        f.write('\n'.join(f'{idx},'+','.join(str(d.get(data, None)) for data in edgeData)
                          +f",{G.nodes[n1]['nodal pressure']},{G.nodes[n2]['nodal pressure']}"
                          +f",{n1},{n2}"
                          +f",{G.nodes[n1]['nodeType']},{G.nodes[n2]['nodeType']}"
                          +f",{G.nodes[n1]['plexus']},{G.nodes[n2]['plexus']}"
                          +f",{int((0.15>np.linalg.norm([G.nodes[n1]['position'], G.nodes[n2]['position']], axis=1)).all())}"                          
                          for idx,(n1,n2,d) in enumerate(G.edges.data()))) # Edges
