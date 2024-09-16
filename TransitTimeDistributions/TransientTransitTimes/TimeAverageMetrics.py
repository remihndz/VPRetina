import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

df = pd.read_csv("./ExampleWaveForm_m.csv")
# df['area'] *= 1e2
# df['velocity'] *= 1e1
# df['flow'] = df['area']*df['velocity']
# df = df[df['time']>3.0]
# df['time'] = df['time']-3.0
df[['time', 'flow']] = df[df.columns]
df['flow'] *= 1e6

RVOFolder = Path().cwd()/"20RVOSimsWithNormalizedFlow"
HealthyFolder = Path().cwd()/"20NormalSimsWithNormalizedFlow"

def ExtractTimeAverage(pathdataFile):
    a = np.fromfile(pathdataFile, dtype=np.float32).reshape(-1,2)
    mask = ~np.logical_or(np.isnan(a), np.isinf(a)).any(1)
    tt, pr = a[mask].T
    assert (tt>0).all(), "Negative times found"
    assert np.logical_and(pr>=0, pr<=1).all(), "Probabilities <0 or >1 found"
    oef = np.exp(-0.004*tt) 
    pr /= pr.sum()
    del a
    
    MTT = (tt*pr).sum() # Expectancy on normalized flow
    df['MTT'] = MTT/df['flow']
    CHT = (((tt-MTT)**2)*pr).sum()**.5 # Standard deviation on normalized flow
    df['CHT'] = CHT/df['flow']
    df['OEF'] = df['flow'].apply(lambda f: 1.0 - (pr*(oef**(1/f))).sum()) # Expectancy = pr.sum()-(pr*exp(-0.004*tt/f))

    # Compute time averages
    MTT  = (df['time'].diff()*df['MTT']).sum()/(df['time'].max()-df['time'].min())
    CHT  = (df['time'].diff()*df['CHT']).sum()/(df['time'].max()-df['time'].min())
    MOEF = (df['time'].diff()*df['OEF']).sum()/(df['time'].max()-df['time'].min())
    
    return MTT, CHT, MOEF

healthyResults = {}
RVOResults = {}
for i, healthyEye in enumerate(HealthyFolder.glob('*.pathdata.bin')):
    print(i, healthyEye)
    MTT, CHT, MOEF = ExtractTimeAverage(healthyEye)
    healthyResults[i] = {'name':healthyEye.name.split('.')[0], 'MTT':MTT, 'CHT':CHT,'MOEF':MOEF}

    print(i, RVOFolder / healthyEye.name)
    MTT, CHT, MOEF = ExtractTimeAverage(RVOFolder / healthyEye.name)
    RVOResults[i] = {'name':healthyEye.name.split('.')[0], 'MTT':MTT, 'CHT':CHT,'MOEF':MOEF}

    pd.DataFrame.from_dict(healthyResults, orient='index').to_csv('TimeAverage_HealthyResults.csv')
    pd.DataFrame.from_dict(RVOResults, orient='index').to_csv('TimeAverage_RVOResults.csv')


