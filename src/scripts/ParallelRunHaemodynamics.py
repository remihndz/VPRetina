import sys, os
import numpy as np
import pandas as pd
import vp_utils as vp
import networkx as nx
from tqdm import tqdm
from parallelbar import progress_map

import metrics
import multiprocessing
from functools import partial

from multiprocessing import set_start_method, get_context, Lock
import traceback

edgeData = ['radius', 'length', 'viscosity', 'hd', 'flow', 'pressure drop']
dirname = os.path.dirname(os.path.realpath(sys.argv[1]))

csvOutFileFolder = os.path.join(dirname, "Haemodynamics-120OPP-1e6R/")
RDummy = 1e6
OPPFolds = 1.2

nPointsSVC=5500
nPointsICP=8000
nPointsDCP=7000

print_lock = Lock()

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def func_imap(iterrow):
    return (iterrow[0], func(iterrow))        

def func(iterrow):

    i, row = iterrow

    MAP, IOP = row['MAP'], row['IOP']
    QCRA = row['vCRA']*np.pi*(row['rCRA']**2)
    pCRA = 2*MAP/3
    pCRV = IOP # row['CRVP']    

    if os.path.isfile(os.path.join(csvOutFileFolder, os.path.basename(row['AVTreeFileName'])[:-4]+'.graph')): 
        G = vp.VirtualRetinalVasculature()
        G.Read(os.path.join(csvOutFileFolder, os.path.basename(row['AVTreeFileName'])[:-4]+'.graph'))
        G.remove_nodes_from(n for n,t in G.nodes.data('nodeType') if t=='dummy')

    elif os.path.isfile(os.path.join(dirname, 'OCTAMetrics', os.path.basename(row['AVTreeFileName'])[:-4]+'.graph')):
        # print("Reading from ", os.path.join(dirname, 'OCTAMetrics', os.path.basename(row['AVTreeFileName'])[:-4]+'.graph'))
        G = vp.VirtualRetinalVasculature()
        G.Read(os.path.join(dirname, 'OCTAMetrics', os.path.basename(row['AVTreeFileName'])[:-4]+'.graph'))
        G.remove_nodes_from(n for n,t in G.nodes.data('nodeType') if t=='dummy')
    
    else:
        try:
            with HiddenPrints():
                G = vp.VirtualRetinalVasculature(row['AVTreeFileName'],
                                                 nPointsSVC=nPointsSVC,
                                                 nPointsICP=nPointsICP,
                                                 nPointsDCP=nPointsDCP,
                                                 )

        except FileNotFoundError:
            print("Did not find file:", row['AVTreeFileName'])
            return 9*[None]

    G._RDummy = RDummy

    with HiddenPrints():
        G.ComputeFlow({'pressure':pCRA}, {'pressure':pCRV})
    G.Write(os.path.join(csvOutFileFolder, os.path.basename(row['AVTreeFileName'])[:-4]+'.graph'))

    fileName = os.path.join(csvOutFileFolder , row['SVCFileName'].split('/')[-1] + f"_{G.MaculaFlow()*100:.0f}.csv")

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

    return G.TRBF(), QCRA, 100*G.MaculaFlow(), pCRA, pCRV, row['rCRA'], nPointsSVC, nPointsICP, nPointsDCP


def Func(data:pd.DataFrame):#, csvOutFileFolder:str):

    Haemodynamics_outcomes = pd.DataFrame(index=data.index, columns=['TRBF', 'Original flow', 'Macula flow', 'pCRA', 'pCRV', 'rCRA', 'nPointsSVC', 'nPointsICP','nPointsDCP'])

    text = f"Process {multiprocessing.Process()._identity}"
    print(' ', end='', flush=True)
    for iterrow in tqdm(list(data.iterrows()), desc=text):

        try:
            Haemodynamics_outcomes.loc[i, ['TRBF', 'Original flow', 'Macula flow', 'pCRA', 'pCRV', 'rCRA', 'nPointsSVC', 'nPointsICP','nPointsDCP']] = func(iterrow)
        except:
            with print_lock:
                print(f"{iterrow[0]}: an exception occured")
                print(traceback.format_exc())
                    
    return Haemodynamics_outcomes
        
if __name__=='__main__':        

    try:
        os.mkdir(csvOutFileFolder)
    except FileExistsError:
        pass
    
    print("Output folder:", csvOutFileFolder)
    os.system("mkdir -p " + csvOutFileFolder)
    population = pd.read_csv(sys.argv[1])
    population['pCRA'] *= OPPFolds
    population['IOP'] *= OPPFolds

    population.to_csv(os.path.join(csvOutFileFolder, "PopulationParameters.csv"))

    Haemodynamics = pd.DataFrame(index=population.index, columns=['TRBF', 'Original flow', 'Macula flow', 'pCRA', 'pCRV', 'rCRA', 'nPointsSVC', 'nPointsICP','nPointsDCP'])

    num_processes = 10

    if num_processes > 1:
        with get_context("spawn").Pool(processes=num_processes) as pool:
            print("Running the sims with", pool._processes, "processes.")
            # Haemodynamics = pd.concat(pool.map(Func, df_split))
            for i, _haemodynamics in tqdm(pool.imap_unordered(func_imap, population.iterrows()),
                                          total=population.shape[0]):
                Haemodynamics.loc[i] = _haemodynamics            

    else:
        print("Running the sims with", 1, "process.")
        for i, _haemodynamics in tqdm(map(func_imap, population.iterrows()),
                                      total=population.shape[0]):
            Haemodynamics.loc[i] = _haemodynamics            

    Haemodynamics.index.name = "sim"
    Haemodynamics.to_csv(os.path.join(csvOutFileFolder, "Haemodynamics.csv"))
    print("Terminated.")
    
