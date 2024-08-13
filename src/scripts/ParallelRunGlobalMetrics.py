import sys
from metrics import *
from math import pi
from tabulate import tabulate
import pandas as pd
import vp_utils as vp
import metrics

from multiprocessing import set_start_method, get_context, Lock
import traceback

FOV = 0.6

nPointsSVC=5500
nPointsICP=8000
nPointsDCP=7000

csvOutFileFolder = os.path.join(*sys.argv[1].split('/')[:-1], f'OCTAMetrics-{FOV}FOV/')
print_lock = Lock()

def func_imap(iterrow):

    metrics = func(iterrow)

    return (iterrow[0], metrics)

def func(iterrow, returnG:bool=False):

    i, params = iterrow
    
    if os.path.isfile(os.path.join(csvOutFileFolder, params['AVTreeFileName'].split('/')[-1][:-4]+'.graph')): 
        G = vp.VirtualRetinalVasculature()
        G.Read(os.path.join(csvOutFileFolder, params['AVTreeFileName'].split('/')[-1][:-4]+'.graph'))
        G.ccoFile = os.path.join(csvOutFileFolder, params['AVTreeFileName'].split('/')[-1][:-4]+'.graph')

    else:
        try:
            with HiddenPrints():
                G = vp.VirtualRetinalVasculature(os.path.join(csvOutFileFolder, params['AVTreeFileName']),
                                                 nPointsSVC=nPointsSVC,
                                                 nPointsICP=nPointsICP,
                                                 nPointsDCP=nPointsDCP,
                                                 ROI = FOV*1.5
                                                 )
                
                G.Write(os.path.join(csvOutFileFolder, params['AVTreeFileName'].split('/')[-1][:-4]+'.graph'))
                G.ccoFile = os.path.join(csvOutFileFolder, params['AVTreeFileName'].split('/')[-1][:-4]+'.graph')
                
        except FileNotFoundError:
            print("Did not find file:", params['AVTreeFileName'], "or", csvOutFileFolder, params['AVTreeFileName'].split('/')[-1][:-4]+'.graph')
            return 13*[None]

    imgFileName_root = os.path.join(csvOutFileFolder,'..', 'Images', params['AVTreeFileName'].split('/')[-1][:-4])
    FD = metrics.FractalDimension(G,FOV=FOV, imageFileName=imgFileName_root+"_FD_SVP.png")
    IVD = metrics.InterVesselDistance(G,FOV=FOV, imageFileName=imgFileName_root+"_IVD_SVP.png")

    rFAZ = min(np.linalg.norm(p[:2]) for n,p in G.nodes.data('position'))
    domainAreaInCm = (FOV**2)#*np.pi/4 # - (rFAZ**2)*pi

    SVP = G.subgraph(n for n,t in G.nodes.data('plexus') if t==0
                     and np.linalg.norm(G.nodes[n]['position'][:2], np.inf)<FOV/2)
    area = sum(d['radius']*d['length'] for u,v,d in SVP.edges.data())*2
    length = sum(l for u,v,l in SVP.edges.data('length'))
    VAD = area/domainAreaInCm
    VSD = length/domainAreaInCm
    VPI = 2e-4*length/domainAreaInCm # cm-1 to micron-1
    VDI = 1e4*area/length # cm to micron
    #VDI = 2e4*np.mean([r for _,_,r in SVP.edges.data('radius')])
    VCI = ((2*length)**2)/(4*np.pi*area)/1.5711 # Normalized as in Chu 2016
    VSD = 1e-4*length/domainAreaInCm
    del SVP
    
    ICP = G.subgraph(n for n,t in G.nodes.data('plexus') if t==1
                     and np.linalg.norm(G.nodes[n]['position'][:2])<FOV/2)
    VAD_ICP = (sum(d['radius']*d['length'] for u,v,d in ICP.edges.data())*2)/(domainAreaInCm)
    IVD_ICP  = metrics.InterVesselDistance(G, FOV=FOV, plexus=1, imageFileName=imgFileName_root+"_IVD_ICP.png")
    del ICP

    DCP = G.subgraph(n for n,t in G.nodes.data('plexus') if t==2
                     and np.linalg.norm(G.nodes[n]['position'][:2])<FOV/2)
    VAD_DCP = (sum(d['radius']*d['length'] for u,v,d in DCP.edges.data())*2)/(domainAreaInCm)
    IVD_DCP  = metrics.InterVesselDistance(G, FOV=FOV, plexus=2, imageFileName=imgFileName_root+"_IVD_DCP.png")
    del DCP

    vesselsICP_DCP, FOVArea = metrics.DVCVessels(G, FOV=FOV)    
    VAD_ICP_DCP = metrics.VesselAreaDensity(vesselsICP_DCP, FOVArea)
    FD_DVC  = metrics.FractalDimension(G, FOV=FOV, plexus=[1,2], imageFileName=imgFileName_root+"_FD_DVC.png")
    IVD_DVC  = metrics.InterVesselDistance(G, FOV=FOV, plexus=[1,2], imageFileName=imgFileName_root+"_IVD_DVC.png")

    if returnG:
        return G, VAD, VPI, VDI, VCI, VSD, FD, IVD, VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DVC, IVD_DVC, rFAZ

    return VAD, VPI, VDI, VCI, VSD, FD, IVD, IVD_ICP, IVD_DCP, VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DVC, IVD_DVC, rFAZ

def Func(data:pd.DataFrame):

    df = pd.DataFrame(index=data.index, columns=['VAD', 'VPI', 'VDI', 'VCI', 'VSD', 'FD', 'IVD', 'IVD_ICP', 'IVD_DCP', 'VAD_ICP_DCP', 'VAD_ICP', 'VAD_DCP', 'FD_DVC', 'IVD_DVC', 'rFAZ'])
    
    with metrics.HiddenPrints():

        for iterrow in data.iterrows():
            try:
                df.loc[i] = func(iterrow, returnG=False)
            except:
                with print_lock:
                  print(f"{iterrow[0]}: an exception occured")
                  print(traceback.format_exc())

    return df


if __name__=='__main__':

    assert(len(sys.argv) > 1)    
    os.system("mkdir -p " + os.path.join(csvOutFileFolder))

    population = pd.read_csv(sys.argv[1], index_col='sim', usecols=['sim', 'AVTreeFileName'], nrows=200)

    num_processes = 6
    # pop_split = np.array_split(population, num_processes)

    metrics = pd.DataFrame(index=population.index, columns=['VAD', 'VPI', 'VDI', 'VCI', 'VSD', 'FD', 'IVD', 'IVD_ICP', 'IVD_DCP', 'VAD_ICP_DCP', 'VAD_ICP', 'VAD_DCP', 'FD_DVC', 'IVD_DVC', 'rFAZ'])

    with get_context("spawn").Pool(processes=num_processes) as pool:
        print("Running the sims with", pool._processes, "processes.")
        with open(os.path.join(csvOutFileFolder, "OCTAMetrics.csv"), 'w') as f:
            f.write(','.join(['sim','VAD', 'VPI', 'VDI', 'VCI', 'VSD', 'FD', 'IVD', 'IVD_ICP', 'IVD_DCP', 'VAD_ICP_DCP', 'VAD_ICP', 'VAD_DCP', 'FD_DVC', 'IVD_DVC', 'rFAZ'])+'\n')
            for i, _metrics in tqdm(pool.imap_unordered(func_imap, population.iterrows()), total=population.shape[0]):
                metrics.loc[i] = _metrics
                f.write(','.join([str(i)]+[f"{m}" for m in metrics.loc[i].values]) + '\n')
        
        #metrics = pd.concat(pool.map(Func, pop_split))
#    del pop_split
    del population
            
    VAD, VPI, VDI, VCI, VSD, FD, IVD, IVD_ICP, IVD_DCP, VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DVC, IVD_DVC, rFAZ = metrics.T.values

    table = tabulate([['VAD', f"{VAD.mean():.4f} ({VAD.std():.4f})", 0.505],
                    ['VPI', f"{VPI.mean():.4f} ({VPI.std():.4f})", 0.05],
                    ['VDI', f"{VDI.mean():.4f} ({VDI.std():.4f})", 24.073],
                    ['VCI', f"{VCI.mean():.4f} ({VCI.std():.4f})", 17962],
                    ['VSD', f"{VSD.mean():.4f} ({VSD.std():.4f})", 0.021],
                    ['FD', f"{FD.mean():.4f} ({FD.std():.4f})", 1.42],
                    ['IVD', f"{IVD.mean():.4f} ({IVD.std():.4f})", 22*1e-4],
                    ['IVD_ICP', f"{IVD_ICP.mean():.4f} ({IVD_ICP.std():.4f})", 22*1e-4],
                    ['IVD_DCP', f"{IVD_DCP.mean():.4f} ({IVD_DCP.std():.4f})", 29*1e-4],
                    ['VAD_ICP', f"{VAD_ICP.mean():.4f} ({VAD_ICP.std():.4f})", 0.21],
                    ['VAD_DCP', f"{VAD_DCP.mean():.4f} ({VAD_DCP.std():.4f})", 0.17],
                    ['VAD_ICP_DCP', f"{VAD_ICP_DCP.mean():.4f} ({VAD_ICP_DCP.std():.4f})", 0.29],
                    ['FD_DVC', f"{FD_DVC.mean():.4f} ({FD_DVC.std():.4f})", 0.9]],
                   headers=['Metric', 'Mean (std)', 'Literature'],
                   tablefmt='simple_grid')
    print(table)

    csvOutFile = os.path.join(csvOutFileFolder,  "OCTAMetrics.csv")
    metrics.index.name = 'sim'
    print("Outcomes written in ", csvOutFile)
    metrics.to_csv(csvOutFile)
