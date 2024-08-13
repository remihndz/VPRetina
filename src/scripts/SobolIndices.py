from scipy.stats import sobol_indices, uniform, randint
from scipy.stats._sensitivity_analysis import sample_A_B, sample_AB
import vp_utils as vp
import metrics 
from SVC.CCO.cpp import svc
from SVC.CCO import cco_utils
from SVC.SSM import ssm_utils
import numpy as np

import time
import os
import pandas as pd


resultsFolder = "../Results/DataForPaper/SobolIndices" # '../Results/SobolIndices'
resultsFolder = os.path.abspath(resultsFolder)
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")
os.system("mkdir -p " + resultsFolder + "/Macula")
simsFolder = os.path.join(resultsFolder, "AVTrees")
os.system("mkdir -p " + simsFolder)

print("Results saved in:", resultsFolder)

verbose = False

d,s,n = 10,7,2**7
print(f"Running {n*(d+2)} simulations.")
# Non hyper-parameters
BaselineParameters = {'rCRA':163/2.0, 'vCRA':6.3, 'MAP':84, 'IOP':14.7, 'CRVP':15, 'capPressure':23,
                      'nTerm':200, # Unused
                      'i':0}

# BaselineParameters = {'nTerm0':600, 'nTerm1':500, 'i':0}

# backbones = np.array([os.path.abspath(f'./Results/DataForPaper/Baseline-Sims-nTerm-Fixed/Coarse/sim_{i}_artery.cco') for i in range(22)])
# for backbone in backbones:
#     os.system(f"cp {backbone} {os.path.join(resultsFolder, 'Backbones')}")    

df = {}
FOV = 0.4

def func(x:list):

    res = []
    for k in range(x.shape[1]):
        res.append(_func(x[:,k]))

    return np.array(res).T

def _func(xi:list):

    lLimFr, delta, eta, gamma, perfAreaFr, thetaMin, closeNeighFr, nTerm0, nTerm1, nPointsSVC = xi
    # backbone = int(backbone)
    
    params = {**BaselineParameters}
    # params['confFileName'] = simsFolder + f"/backbone_{backbone}_{BaselineParameters['i']}.conf"
    params['confFileNameMacula'] = f"{resultsFolder}/Macula/tmp_macula"

    params["backboneFileName"] = "../Results/DataForPaper/Baseline-Sims-nTerm-Fixed/Coarse/sim_0" # backbones[backbone]
    #params['outputFileName'] = simsFolder + f"/backbone_{backbone}_{BaselineParameters['i']}"
    params['lLimFr'] = lLimFr
    params['delta'] = delta
    params['eta'] = eta
    params['gamma'] = gamma
    params['perfAreaFr'] = perfAreaFr
    params['thetaMin'] = thetaMin
    params['closeNeighFr'] = closeNeighFr
    params['nTerm0'] = int(nTerm0)
    params['nTerm1'] = int(nTerm1)
    params['nPointsSVC'] = int(nPointsSVC)
    params['SVCFileName'] = resultsFolder + f"/Macula/sim_{params['i']}"
    params['AVTreeFileName'] = f"{resultsFolder}/AVTrees/sim_{params['i']}_AV.cco"

    params['confFileNameMacula_artery'] = params['confFileNameMacula'] + '_artery.conf'
    params['confFileNameMacula_vein'] = params['confFileNameMacula'] + '_vein.conf'
        
    cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula_artery'],
                             backboneFileName=params['backboneFileName']+'_artery.cco',
                             outputFileName=params['SVCFileName']+'_artery',
                             **{key:param for key,param in params.items() if 'filename' not in key.lower()})

    cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula_vein'],
                             backboneFileName=params['backboneFileName']+'_vein.cco',
                             outputFileName=params['SVCFileName']+'_vein',
                             **{key:param for key,param in params.items() if 'filename' not in key.lower()})
    
    svc.SVC_macula(params['confFileNameMacula_artery'], verbose)
    svc.SVC_macula(params['confFileNameMacula_vein'], verbose)    
    mergedTree = cco_utils.MergeArteryVeinCCO(arteryFile=params['SVCFileName']+'_artery.cco',
                                              veinFile=params['SVCFileName']+'_vein.cco',
                                              outputFile=params['AVTreeFileName'])
    
    with metrics.HiddenPrints():

        G = vp.VirtualRetinalVasculature(params['AVTreeFileName'],
                                         nPointsSVC=params['nPointsSVC'],
                                         nPointsICP=0, nPointsDCP=0)

        FD = metrics.FractalDimension(G,FOV)
        IVD = metrics.InterVesselDistance(G,FOV=FOV)
        vessels = metrics.SVPVessels(G, FOV=FOV)
        vesselDensity = metrics.VesselAreaDensity(vessels, FOV)
        VAD = metrics.VesselAreaDensity(vessels, FOV)
        VPI = metrics.VesselPerimeterIndex(vessels, FOV)
        VDI = metrics.VesselDiameterIndex(vessels, FOV)
        VCI = metrics.VesselComplexityIndex(vessels, FOV)
        VSD = metrics.VesselSkeletonDensity(vessels, FOV)

        df[BaselineParameters['i']] = {**params}
        df[BaselineParameters['i']]['FD'] = FD
        df[BaselineParameters['i']]['IVD'] = IVD
        df[BaselineParameters['i']]['VAD'] = VAD
        df[BaselineParameters['i']]['VPI'] = VPI
        df[BaselineParameters['i']]['VDI'] = VDI
        df[BaselineParameters['i']]['VCI'] = VCI
        df[BaselineParameters['i']]['VSD'] = VSD

        pd.DataFrame.from_dict(df, orient='index').to_csv(resultsFolder+'/MetricsAndParametersData.csv')
        
    BaselineParameters['i'] = BaselineParameters['i'] + 1 # Just making sure we are not overwriting sims.

    return [FD, IVD, VAD, VPI, VDI, VCI, VSD]

rng = np.random.default_rng()

dists = [
    uniform(loc=0.1, scale = 0.8), # lLimFr
    uniform(loc=0.1, scale = 0.8), # delta
    uniform(loc=0.1, scale = 0.5), # eta
    uniform(loc=2.5, scale = 0.5), # gamma
    uniform(loc=0.1, scale = 0.8), # perfAreaFr
    uniform(loc=0., scale = 0.4), # thetaMin
    uniform(loc=0., scale = 5), # closeNeighFr
#    randint(0, len(backbones)), # Backbones
    randint(300, 500),          # nTerm0
    randint(200,400),           # nTerm1
    randint(300,700)            # nPointsSVC
]

indices = sobol_indices(
    func=func,
    n = n,
    dists=dists,
    random_state=rng
)

print("First order:\n\t", indices.first_order)
print("Total order:\n\t", indices.total_order)

metrics = ['FD','IVD','VAD','VPI','VDI','VCI','VSD']
parameters = ["lLimFr", "delta", "eta", "gamma", "perfAreaFr", "thetaMin", "closeNeighFr", "backbone", "nTerm0", "nTerm1", "nPointsSVC"]

pd.DataFrame.from_dict({metric:{param:value for param, value in zip(parameters, values)}
                        for metric, values in zip(metrics, indices.first_order)}).to_csv(resultsFolder+'firstOrder.csv')

pd.DataFrame.from_dict({metric:{param:value for param, value in zip(parameters, values)}
                        for metric, values in zip(metrics, indices.total_order)}).to_csv(resultsFolder+'totalOrder.csv')


pd.DataFrame.from_dict(df, orient='index').to_csv(resultsFolder+'/MetricsAndParametersData.csv')

