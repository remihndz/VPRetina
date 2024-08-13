import vp_utils as vp
from SVC.CCO.cpp import svc
from SVC.CCO import cco_utils
from SVC.SSM import ssm_utils
import numpy as np
import matplotlib.pyplot as plt
import time
from tabulate import tabulate
import os
import pandas as pd

verbose = False

m,n = 3,10

resultsFolder = '../Results/SensitivityAnalysis/lLimFr'
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")

values = np.linspace(0.1, 1., n, endpoint=False)
print(values)

for i in range(m):
    os.system(f"cp ../Results/DataForConference/Macula/sim_{i}.cco " + resultsFolder+f"/baseline_{i}.cco")
    os.system(f"cp ../Results/DataForConference/Coarse/sim_{i}_AV.cco " + resultsFolder+f"/Backbones/backbone_{i}.cco")


times = []
nTerms = []

parameters = []
fullListOfParams = {}
# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

for i in range(n):
    print(f"Simulations {i}...")
    for j in range(m):
    
        t0 = time.time(); print(t0)
        params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=[100, 500, 600])
        
        params['nTerm0'] = 500
        params['nTerm1'] = 600
        
        headers = [k for k in params.keys() if "FileName" not in k]
        
        params['lengthCRA'] = 0.1
        
        parameters.append([params[k] for k in headers])

        # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))

        params['confFileNameMacula']=resultsFolder + f"/backbone_{j}_{values[i]}.conf"
        params['backboneFileName']=resultsFolder+f"/Backbones/backbone_{j}.cco"
        params['SVCFileName']=resultsFolder+f"/backbone_{j}_{values[i]}"
        params['lLimFr']=values[i]
        fullListOfParams[i*m + j] = params

        cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula'],
                                 backboneFileName=params['backboneFileName'],
                                 outputFileName=params['SVCFileName'],
                                 nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                                 lLimFr=values[i],        
                                 thetaMin=0.224) # Capillary bifurcation angles

        # svc.SVC_macula(params['confFileNameMacula'], verbose)
        times.append(time.time()-t0)
        nTerms.append(params['nTerm1']+params['nTerm'])

        # Save current population in case of a crash
        df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
        df.to_csv(resultsFolder + '/PopulationParameters.csv')
        

resultsFolder = '../Results/SensitivityAnalysis/delta'
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")

values = np.linspace(0.1, 1., n, endpoint=False)
print(values)

for i in range(m):
    os.system(f"cp ../Results/DataForConference/Macula/sim_{i}.cco " + resultsFolder+f"/baseline_{i}.cco")
    os.system(f"cp ../Results/DataForConference/Coarse/sim_{i}_AV.cco " + resultsFolder+f"/Backbones/backbone_{i}.cco")


times = []
nTerms = []

parameters = []
fullListOfParams = {}
# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

for i in range(n):
    print(f"Simulations {i}...")
    for j in range(m):
    
        t0 = time.time(); print(t0)
        params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=[100, 500, 600])
        
        params['nTerm0'] = 500
        params['nTerm1'] = 600
        
        headers = [k for k in params.keys() if "FileName" not in k]
        
        params['lengthCRA'] = 0.1
        
        parameters.append([params[k] for k in headers])

        # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))

        params['confFileNameMacula']=resultsFolder + f"/backbone_{j}_{values[i]}.conf"
        params['backboneFileName']=resultsFolder+f"/Backbones/backbone_{j}.cco"
        params['SVCFileName']=resultsFolder+f"/backbone_{j}_{values[i]}"
        params['lLimFr']=values[i]
        fullListOfParams[i*m + j] = params

        cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula'],
                                 backboneFileName=params['backboneFileName'],
                                 outputFileName=params['SVCFileName'],
                                 nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                                 delta=values[i],        
                                 thetaMin=0.224) # Capillary bifurcation angles

        svc.SVC_macula(params['confFileNameMacula'], verbose)
        times.append(time.time()-t0)
        nTerms.append(params['nTerm1']+params['nTerm'])

        # Save current population in case of a crash
        df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
        df.to_csv(resultsFolder + '/PopulationParameters.csv')
        

    
resultsFolder = '../Results/SensitivityAnalysis/gamma'
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")

values = np.linspace(2.5,3, n, endpoint=True)
print(values)

for i in range(m):
    os.system(f"cp ../Results/DataForConference/Macula/sim_{i}.cco " + resultsFolder+f"/baseline_{i}.cco")
    os.system(f"cp ../Results/DataForConference/Coarse/sim_{i}_AV.cco " + resultsFolder+f"/Backbones/backbone_{i}.cco")


times = []
nTerms = []

parameters = []
fullListOfParams = {}
# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

for i in range(n):
    print(f"Simulations {i}...")
    for j in range(m):
    
        t0 = time.time()
        print(t0)
        params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=[100, 500, 600])
        
        params['nTerm0'] = 500
        params['nTerm1'] = 600
        
        headers = [k for k in params.keys() if "FileName" not in k]
        
        params['lengthCRA'] = 0.1
        
        parameters.append([params[k] for k in headers])

        # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))

        params['confFileNameMacula']=resultsFolder + f"/backbone_{j}_{values[i]}.conf"
        params['backboneFileName']=resultsFolder+f"/Backbones/backbone_{j}.cco"
        params['SVCFileName']=resultsFolder+f"/backbone_{j}_{values[i]}"
        params['lLimFr']=values[i]
        fullListOfParams[i*m + j] = params

        cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula'],
                                 backboneFileName=params['backboneFileName'],
                                 outputFileName=params['SVCFileName'],
                                 nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                                 gamma=values[i],        
                                 thetaMin=0.224) # Capillary bifurcation angles

        svc.SVC_macula(params['confFileNameMacula'], verbose)
        times.append(time.time()-t0)
        nTerms.append(params['nTerm1']+params['nTerm'])

        # Save current population in case of a crash
        df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
        df.to_csv(resultsFolder + '/PopulationParameters.csv')
        

resultsFolder = '../Results/SensitivityAnalysis/eta'
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")

values = np.linspace(0.1, 0.6, n, endpoint=True)
print(values)

for i in range(m):
    os.system(f"cp ../Results/DataForConference/Macula/sim_{i}.cco " + resultsFolder+f"/baseline_{i}.cco")
    os.system(f"cp ../Results/DataForConference/Coarse/sim_{i}_AV.cco " + resultsFolder+f"/Backbones/backbone_{i}.cco")


times = []
nTerms = []

parameters = []
fullListOfParams = {}
# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

for i in range(n):
    print(f"Simulations {i}...")
    for j in range(m):
    
        t0 = time.time(); print(t0)
        params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=[100, 500, 600])
        
        params['nTerm0'] = 500
        params['nTerm1'] = 600
        
        headers = [k for k in params.keys() if "FileName" not in k]
        
        params['lengthCRA'] = 0.1
        
        parameters.append([params[k] for k in headers])

        # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))

        params['confFileNameMacula']=resultsFolder + f"/backbone_{j}_{values[i]}.conf"
        params['backboneFileName']=resultsFolder+f"/Backbones/backbone_{j}.cco"
        params['SVCFileName']=resultsFolder+f"/backbone_{j}_{values[i]}"
        params['lLimFr']=values[i]
        fullListOfParams[i*m + j] = params

        cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula'],
                                 backboneFileName=params['backboneFileName'],
                                 outputFileName=params['SVCFileName'],
                                 nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                                 eta=values[i],        
                                 thetaMin=0.224) # Capillary bifurcation angles
        
        svc.SVC_macula(params['confFileNameMacula'], verbose)
        times.append(time.time()-t0)
        nTerms.append(params['nTerm1']+params['nTerm'])

        # Save current population in case of a crash
        df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
        df.to_csv(resultsFolder + '/PopulationParameters.csv')
        

resultsFolder = '../Results/SensitivityAnalysis/thetaMin'
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")

values = np.linspace(0.0, 0.4, n, endpoint=True)
print(values)

for i in range(m):
    os.system(f"cp ../Results/DataForConference/Macula/sim_{i}.cco " + resultsFolder+f"/baseline_{i}.cco")
    os.system(f"cp ../Results/DataForConference/Coarse/sim_{i}_AV.cco " + resultsFolder+f"/Backbones/backbone_{i}.cco")


times = []
nTerms = []

parameters = []
fullListOfParams = {}
# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

for i in range(n):
    print(f"Simulations {i}...")
    for j in range(m):
    
        t0 = time.time(); print(t0)
        params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=[100, 500, 600])
        
        params['nTerm0'] = 500
        params['nTerm1'] = 600
        
        headers = [k for k in params.keys() if "FileName" not in k]
        
        params['lengthCRA'] = 0.1
        
        parameters.append([params[k] for k in headers])

        # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))

        params['confFileNameMacula']=resultsFolder + f"/backbone_{j}_{values[i]}.conf"
        params['backboneFileName']=resultsFolder+f"/Backbones/backbone_{j}.cco"
        params['SVCFileName']=resultsFolder+f"/backbone_{j}_{values[i]}"
        params['lLimFr']=values[i]
        fullListOfParams[i*m + j] = params

        cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula'],
                                 backboneFileName=params['backboneFileName'],
                                 outputFileName=params['SVCFileName'],
                                 nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                                 thetaMin=values[i])        
        #thetaMin=0.224) # Capillary bifurcation angles
        
        svc.SVC_macula(params['confFileNameMacula'], verbose)
        times.append(time.time()-t0)
        nTerms.append(params['nTerm1']+params['nTerm'])

        # Save current population in case of a crash
        df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
        df.to_csv(resultsFolder + '/PopulationParameters.csv')


resultsFolder = '../Results/SensitivityAnalysis/perfAreaFr'
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")

values = np.linspace(0.1, 1., n, endpoint=False)
print(values)

for i in range(m):
    os.system(f"cp ../Results/DataForConference/Macula/sim_{i}.cco " + resultsFolder+f"/baseline_{i}.cco")
    os.system(f"cp ../Results/DataForConference/Coarse/sim_{i}_AV.cco " + resultsFolder+f"/Backbones/backbone_{i}.cco")


times = []
nTerms = []

parameters = []
fullListOfParams = {}
# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]

for i in range(n):
    print(f"Simulations {i}...")
    for j in range(m):
    
        t0 = time.time(); print(t0)
        params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=[100, 500, 600])
        
        params['nTerm0'] = 500
        params['nTerm1'] = 600
        
        headers = [k for k in params.keys() if "FileName" not in k]
        
        params['lengthCRA'] = 0.1
        
        parameters.append([params[k] for k in headers])

        # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))

        params['confFileNameMacula']=resultsFolder + f"/backbone_{j}_{values[i]}.conf"
        params['backboneFileName']=resultsFolder+f"/Backbones/backbone_{j}.cco"
        params['SVCFileName']=resultsFolder+f"/backbone_{j}_{values[i]}"
        params['lLimFr']=values[i]
        fullListOfParams[i*m + j] = params

        cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula'],
                                 backboneFileName=params['backboneFileName'],
                                 outputFileName=params['SVCFileName'],
                                 nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                                 perfAreaFr=values[i],        
                                 thetaMin=0.224) # Capillary bifurcation angles
        
        svc.SVC_macula(params['confFileNameMacula'], verbose)
        times.append(time.time()-t0)
        nTerms.append(params['nTerm1']+params['nTerm'])

        # Save current population in case of a crash
        df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
        df.to_csv(resultsFolder + '/PopulationParameters.csv')
        

resultsFolder = '../Results/SensitivityAnalysis/closeNeighFr'
os.system("mkdir -p " + resultsFolder + "/Images")
os.system("mkdir -p " + resultsFolder + "/Backbones")

values = np.linspace(0.0, 5.0, n, endpoint=True)
print(values)

for i in range(m):
    os.system(f"cp ../Results/DataForConference/Macula/sim_{i}.cco " + resultsFolder+f"/baseline_{i}.cco")
    os.system(f"cp ../Results/DataForConference/Coarse/sim_{i}_AV.cco " + resultsFolder+f"/Backbones/backbone_{i}.cco")


times = []
nTerms = []

parameters = []
fullListOfParams = {}
# headers = ["rCRA","vCRA", "nTerm", "nTerm0", "nTerm1", "SVCFileName"]


for i in range(n):
    print(f"Simulations {i}...")
    for j in range(m):
    
        t0 = time.time(); print(t0)
        params = vp.GenerateParameters(resultsFolder, n=i, baseTerm=[100, 500, 600])
        
        params['nTerm0'] = 500
        params['nTerm1'] = 600
        
        headers = [k for k in params.keys() if "FileName" not in k]
        
        params['lengthCRA'] = 0.1
        
        parameters.append([params[k] for k in headers])

        # print(tabulate([[params[k] for k in headers]], headers=headers, tablefmt="fancy_grid"))

        params['confFileNameMacula']=resultsFolder + f"/backbone_{j}_{values[i]}.conf"
        params['backboneFileName']=resultsFolder+f"/Backbones/backbone_{j}.cco"
        params['SVCFileName']=resultsFolder+f"/backbone_{j}_{values[i]}"
        params['lLimFr']=values[i]
        fullListOfParams[i*m + j] = params

        cco_utils.ConfFileStage1(confFileName=params['confFileNameMacula'],
                                 backboneFileName=params['backboneFileName'],
                                 outputFileName=params['SVCFileName'],
                                 nTerm0=params['nTerm0'],nTerm1=params['nTerm1'],
                                 closeNeighFr=values[i],        
                                 thetaMin=0.224) # Capillary bifurcation angles

        svc.SVC_macula(params['confFileNameMacula'], verbose)
        times.append(time.time()-t0)
        nTerms.append(params['nTerm1']+params['nTerm'])

        # Save current population in case of a crash
        df = pd.DataFrame.from_dict({i:fullListOfParams[i] for i in range(len(fullListOfParams))}, orient='index')
        df.to_csv(resultsFolder + '/PopulationParameters.csv')
        


#print(tabulate(parameters, headers=headers, tablefmt='fancy_grid'))
plt.bar(nTerms, times);
plt.xlabel('# of terminal vessels')
plt.ylabel('Time (s)')
plt.show()
