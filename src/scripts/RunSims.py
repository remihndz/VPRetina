import vp_utils as vp
from SVC.CCO.cpp import svc
from SVC.CCO import cco_utils
from SVC.SSM import ssm_utils
import numpy as no
import matplotlib.pyplot as plt
import time
from tabulate import tabulate
import os
import pandas as pd
import numpy as np

verbose = False

n = 2
baseTerm = [200, 400, 300]
capPressure=35

times = []
nTerms = []


resultsFolder=os.path.abspath(f"../Results/TestingSpeedWith1Core")

os.system("mkdir -p " + os.path.join(resultsFolder, "SSM"))
os.system("mkdir -p " + os.path.join(resultsFolder, "Coarse"))
os.system("mkdir -p " + os.path.join(resultsFolder, "Macula"))
os.system("mkdir -p " + os.path.join(resultsFolder, "AVTrees"))
os.system("mkdir -p " + os.path.join(resultsFolder, "Images"))

ssm = ssm_utils.CreateShapeModels('../../data/Shape_Dataset/')
parameters = []
fullListOfParams = {}

# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

for i in range(n):

    print(f"Simulation {i}...")
    t0 = time.time()
    params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=baseTerm, capPressure=capPressure)
    params['CRVP'] = max(params['CRVP'], params['IOP'])
    params['pCRA'] = (2./3.)*params.get('MAP',84) - params.get('IOP', 15) # Guidoboni 2014
    # params['pCRA'] = np.random.normal(52.8, 9.3)
    params['nTerm'], params['nTerm0'], params['nTerm1'] = baseTerm

    headers = [k for k in params.keys() if "FileName" not in k]
    params['lengthCRA'] = 0.1
    parameters.append([params[k] for k in headers])
    fullListOfParams[i] = params
    # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))
    ssm.Generate(c=0.8, saveIn=params['ssmFileName'], **params)

    params['delta'] = 0.0

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

    params['confFileNameMacula_artery'] = params['confFileNameMacula'] + '_artery.conf'
    params['confFileNameMacula_vein'] = params['confFileNameMacula'] + '_vein.conf'

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
                                              outputFile=params['AVTreeFileName'],)


    times.append(time.time()-t0)
    print("Sim time:", times[-1])
    nTerms.append(2*(params['nTerm0']+params['nTerm1']+params['nTerm']))

    df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
    df.to_csv(resultsFolder + '/PopulationParameters.csv')

print(tabulate(parameters, headers=headers, tablefmt='fancy_grid'))
print("Simulation saved in " + resultsFolder)

print("Average simulation time: ", np.mean(times))
# plt.bar(nTerms, times);
# plt.show()
