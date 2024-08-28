import sys
sys.path.append("/home/remi/VPRetina/src/")
sys.path.append("/home/remi/VPRetina/src/SVC/CCO")
sys.path.append("/home/remi/VPRetina/src/SVC/SSM")
sys.path.append("/home/remi/VPRetina/src/SVC/CCO")

import vp_utils as vp
import ssm_utils
import cco_utils
import svc 

# from SVC.CCO.cpp import svc
# from SVC.CCO import cco_utils
# from SVC.SSM import ssm_utils
import matplotlib.pyplot as plt
import numpy as np
import time
from tabulate import tabulate
import os
import pandas as pd
from tqdm import tqdm
from SALib.sample import sobol
import seaborn as sns
from pathlib import Path

import metrics
import multiprocessing
from functools import partial

from multiprocessing import set_start_method
from multiprocessing import get_context

verbose = False
FOV = 0.3
resultsFolder = os.path.abspath('/media/Storage3.6TB/VariedPopulationOfRetinas')


def _func(params):
    return params['sim'],0,1,2,3,4,5,6

def computeMetrics(params, returnG:bool=False):

    #i, params = iterrow
    
    if 0==1 and os.path.isfile(os.path.join(resultsFolder, 'OCTAMetrics', params['AVTreeFileName'].split('/')[-1][:-4]+'.graph')): 
        G = vp.VirtualRetinalVasculature()
        G.Read(os.path.join(resultsFolder, 'OCTAMetrics', params['AVTreeFileName'].split('/')[-1][:-4]+'.graph'))
        G.ccoFile = os.path.join(resultsFolder, 'OCTAMetrics', params['AVTreeFileName'].split('/')[-1][:-4]+'.graph')

    else:
        try:
            with metrics.HiddenPrints():
                G = vp.VirtualRetinalVasculature(os.path.join(resultsFolder, 'OCTAMetrics', params['AVTreeFileName']),
                                                 nPointsSVC=params['nPointsSVC'],
                                                 nPointsICP=params['nPointsICP'],
                                                 nPointsDCP=params['nPointsDCP'],
                                                 ROI = FOV*1.5
                                                 )
                
                G.Write(os.path.join(resultsFolder, 'OCTAMetrics', params['AVTreeFileName'].split('/')[-1][:-4]+'.graph'))
                G.ccoFile = os.path.join(resultsFolder, 'OCTAMetrics', params['AVTreeFileName'].split('/')[-1][:-4]+'.graph')
                
        except FileNotFoundError:
            print("Did not find file:", params['AVTreeFileName'], "or", os.path.join(resultsFolder, 'OCTAMetrics', params['AVTreeFileName'].split('/')[-1][:-4]+'.graph'))
            return 13*[None]

    imgFileName_root = os.path.join(resultsFolder, 'Images', params['AVTreeFileName'].split('/')[-1][:-4])
    FD = metrics.FractalDimension(G,FOV=FOV, imageFileName=imgFileName_root+"_SVP.png")
    IVD = metrics.InterVesselDistance(G,FOV=FOV, imageFileName=imgFileName_root+"_SVP.png")

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
    del ICP

    DCP = G.subgraph(n for n,t in G.nodes.data('plexus') if t==2
                     and np.linalg.norm(G.nodes[n]['position'][:2])<FOV/2)
    VAD_DCP = (sum(d['radius']*d['length'] for u,v,d in DCP.edges.data())*2)/(domainAreaInCm)
    del DCP

    vesselsICP_DCP, FOVArea = metrics.DVCVessels(G, FOV=FOV)    
    VAD_ICP_DCP = metrics.VesselAreaDensity(vesselsICP_DCP, FOVArea)
    FD_DVC  = metrics.FractalDimension(G, FOV=FOV, plexus=[1,2], imageFileName=imgFileName_root+"_DVC.png")
    IVD_DVC  = metrics.InterVesselDistance(G, FOV=FOV, plexus=[1,2], imageFileName=imgFileName_root+"_DVC.png")

    if returnG:
        return G, VAD, VPI, VDI, VCI, VSD, FD, IVD, VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DVC, IVD_DVC, rFAZ
    
    return VAD, VPI, VDI, VCI, VSD, FD, IVD, VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DVC, IVD_DVC, rFAZ


def Func(params):

    t0 = time.time()
    cco_utils.ConfFileStage0(confFileName=params['confFileNameVein'],
                             backboneFileName=params['ssmFileName']+'_vein.cco',
                             outputFileName=params['veinFileName'],
                             # nTerm=params['nTerm'],
                             # thetaMin=0.2205,  # Vein bifurcation angles
                             # perfAreaFr=0.4,
                             **params)
    cco_utils.ConfFileStage0(confFileName=params['confFileNameArtery'],
                             backboneFileName=params['ssmFileName']+'_artery.cco',
                             outputFileName=params['arteryFileName'],
                             # nTerm=params['nTerm'],
                             # thetaMin=0.232, # Artery bifurcation angles
                             # perfAreaFr=0.4,
                             **params)

    svc.SVC_coarse(params['confFileNameArtery'], verbose)
    svc.SVC_coarse(params['confFileNameVein'], verbose)

    cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula_artery'],
                             backboneFileName=params['arteryFileName']+'.cco',
                             outputFileName=params['SVCFileName']+'_artery',
                             nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                             thetaMin=0.224) # Capillary bifurcation angles

    cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula_vein'],
                             backboneFileName=params['veinFileName']+'.cco',
                             outputFileName=params['SVCFileName']+'_vein',
                             nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                             thetaMin=0.224) # Capillary bifurcation angles

    svc.SVC_macula(params['confFileNameMacula_artery'], verbose)
    svc.SVC_macula(params['confFileNameMacula_vein'], verbose)

    mergedTree = cco_utils.MergeArteryVeinCCO(arteryFile=params['SVCFileName']+'_artery.cco',
                                              veinFile=params['SVCFileName']+'_vein.cco',
                                              outputFile=params['AVTreeFileName'])
    
    with metrics.HiddenPrints():

        G = vp.VirtualRetinalVasculature(params['AVTreeFileName'],
                                         nPointsSVC=params['nPointsSVC'],
                                         nPointsICP=params['nPointsICP'],
                                         nPointsDCP=params['nPointsDCP'],
                                         pCRA=params['pCRA'],
                                         pCRV=params['pCRV'],
                                         )
        G._RDummy = 1e6

        FD = metrics.FractalDimension(G,FOV, imageFile=os.path.join(resultsFolder, "Images", G.ccoFile[:-4]+"_SVP.png"))
        IVD = metrics.InterVesselDistance(G,FOV=FOV, imageFile=os.path.join(resultsFolder, "Images", G.ccoFile[:-4]+"_SVP.png"))
        vessels, _ = metrics.SVPVessels(G, FOV=FOV)
        vesselDensity = metrics.VesselAreaDensity(vessels, FOV)
        VAD = metrics.VesselAreaDensity(vessels, FOV)
        VPI = metrics.VesselPerimeterIndex(vessels, FOV)
        VDI = metrics.VesselDiameterIndex(vessels, FOV)
        VCI = metrics.VesselComplexityIndex(vessels, FOV)
        VSD = metrics.VesselSkeletonDensity(vessels, FOV)

        vesselsICP, _ = metrics.ICPVessels(G, FOV=FOV)
        vesselsDCP, _ = metrics.DCPVessels(G, FOV=FOV)
        vesselsICP_DCP, _ = metrics.DVCVessels(G, FOV=FOV)

        VAD_ICP_DCP = metrics.VesselAreaDensity(vesselsICP_DCP, FOV)
        VAD_DCP = metrics.VesselAreaDensity(vesselsDCP, FOV)
        VAD_ICP = metrics.VesselAreaDensity(vesselsICP, FOV)
        FD_DCP  = metrics.FractalDimension(G, FOV, plexus=2, imageFile=os.path.join(resultsFolder, "Images", G.ccoFile[:-4]+"_DCP.png"))

        G, VAD, VPI, VDI, VCI, VSD, FD, IVD, VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DVC, IVD_DVC, rFAZ = computeMetrics(params, returnG=True)
        G._RDummy = 1e6
        
        G.ComputeFlow()#{'pressure':params['pCRA']}, {'pressure':params['pCRV']})
        TRBF = G.TRBF()
        MaculaFlow = G.MaculaFlow()
        del G
        
    return (params['sim'], FD, IVD,
            VAD, VPI, VDI, VCI, VSD,
            VAD_ICP_DCP,VAD_ICP,
            VAD_DCP,FD_DVC,
            TRBF,MaculaFlow)
    # return (params['sim'], *13*[0])

if __name__=='__main__':

    # Run in multiprocessing
    num_processes = multiprocessing.cpu_count()//int(os.getenv('OMP_NUM_THREADS', 4))  # Adjust as needed

    BaselineParameters = {'rCRA':163/2.0, 'vCRA':6.3, 'MAP':84, 'IOP':14.7, 'pCRV':15, 'capPressure':23,
                          'nTerm':20, # Unused
                          'i':0}

    baseTerm = [200,150,75]
    capPressure=35

    times = []
    nTerms = []

    n = 500

    os.system("mkdir -p " + os.path.join(resultsFolder, "SSM"))
    os.system("mkdir -p " + os.path.join(resultsFolder, "Coarse"))
    os.system("mkdir -p " + os.path.join(resultsFolder, "Macula"))
    os.system("mkdir -p " + os.path.join(resultsFolder, "AVTrees"))
    os.system("mkdir -p " + os.path.join(resultsFolder, "Images"))
    os.system("mkdir -p " + os.path.join(resultsFolder, "OCTAMetrics"))

    with metrics.HiddenPrints():
        ssm = ssm_utils.CreateShapeModels('../../data/Shape_Dataset/')
    parameters = []
    fullListOfParams = {}

    # headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

    problem = {
        'num_vars':12,
        'names':["lLimFr", "delta", "eta", "gamma", "perfAreaFr", "thetaMin", "closeNeighFr",
                 "nTerm0", "nTerm1", 
                 "nPointsSVC", "nPointsICP", "nPointsDCP"],
        'bounds':[[0.1, 0.9],
                  [0.1, 0.9],
                  [0.2, 0.5],
                  [2., 3.0],
                  [0.1, 0.9],
                  [0.0, 0.4],
                  [0.0, 5.0],
                  [150*.8, 150*1.2], # nTerm0
                  [75*.8, 75*1.2], # nTerm1
                  [5500*.8, 5500*1.2], # nPointsSVC
                  [16000*.8, 16000*1.2], # nPointsICP
                  [10500*.8, 10500*1.2]] # nPointsDCP
    }
    params_value = sobol.sample(problem, n, calc_second_order=False)
    np.random.shuffle(params_value)
    params_value = params_value[:n]
    print(params_value.shape[0], "simulations running for an estimated", int(params_value.shape[0]*300/60/60/24/4), "days")
    capPressures = np.random.normal(capPressure, capPressure/15., size=n)

    # n = params_value.shape[0]
    
    ## Create the parameters for the sims.
    for i in tqdm(range(n), desc='Creating the parameter samples'):
        
        params = vp.GenerateParameters(resultsFolder, n=i, capPressure=capPressures[i],
                                       baseTerm=baseTerm) 
        params['sim'] = i
        params['pCRV'] = max(params['pCRV'], params['IOP'])
        params['pCRA'] = (2./3.)*params.get('MAP',84) - params.get('IOP', 15) # Guidoboni 2014
        # params['pCRA'] = np.random.normal(52.8, 9.3)

        params['lengthCRA'] = 0.1
        
        params['confFileNameMacula_artery'] = params['confFileNameMacula'] + f'sim_{i}_artery.conf'
        params['confFileNameMacula_vein'] = params['confFileNameMacula'] + f'sim_{i}_vein.conf'
        params['confFileNameVein'] = f"{resultsFolder}/Coarse/tmp_sim_{i}_vein.conf"
        params['confFileNameArtery'] = f"{resultsFolder}/Coarse/tmp_sim_{i}_artery.conf"

        varyingParams = {key:(val if ("nPoints" not in key
                                     and "nTerm" not in key)
                              else int(val))
                         for key, val in zip(problem['names'], params_value[i])}
        params = {**params, **varyingParams} # Add the parameters sampled using Sobol sampling

        if params['capPressure'] >= params['pCRA']:
            params['capPressure'] = params['pCRA'] - (params['pCRA']-params['pCRV'])
        
        fullListOfParams[i] = params
        
        ssm.Generate(c=0.8, saveIn=params['ssmFileName'], **params)

    df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
    df.to_csv(resultsFolder + '/PopulationParameters.csv')

    # sns.displot(data=df.melt(value_vars=['rCRA','vCRA','MAP','IOP','pCRV']),
    #             col='variable', kind='kde')
    # plt.show()

    del ssm, params, df,
    
    with get_context("spawn").Pool(processes=num_processes) as pool:
        print("Running the sims with", pool._processes, "processes.")
        # #pool.map(Func, list(fullListOfParams.values()))
        with open(os.path.join(resultsFolder, "outputs.csv"), 'w') as f:
            f.write("sim,FD,IVD,VAD,VPI,VDI,VCI,VSD,VAD_ICP_DCP,VAD_ICP,VAD_DCP,FD_DVC,TRBF,MaculaFlow\n")
            for results in tqdm(pool.imap_unordered(Func, list(fullListOfParams.values())),
                                desc='Running simulations',
                                total=len(fullListOfParams)):
                f.write(','.join(str(x) for x in results) + '\n')

    (Path(resultsFolder) / "OCTAMetrics").mkdir(exist_ok=True, parents=True) # Create directory 
    os.system(f"cp {os.path.join(resultsFolder, 'outputs.csv')} {os.path.join(resultsFolder, OCTAMetrics, 'OCTAMetrics.csv')}")
