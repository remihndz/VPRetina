import sys
sys.path.append("../")
from metrics import *
from math import pi
from tabulate import tabulate
import pandas as pd
import seaborn as sns

FOV = 0.3

assert(len(sys.argv) > 1)
listOfFiles = sys.argv[1:]

# Give a value to use for all networks if wanted. Otherwise, the FAZ's area is substracted for each network.
domainAreaInCm = None # FOV**2 - (0.02**2)*pi

csvOutFile = os.path.relpath('/'.join(listOfFiles[0].split('/')[:-1]) + '/../OCTAMetrics.csv')
print(csvOutFile)

metrics = StatisticsMultipleTrees(listOfFiles, FOV=FOV,
                                  domainAreaInCm=domainAreaInCm,
                                  plot=False,
                                  nPointsSVC=6000,
                                  nPointsICP=5000,
                                  nPointsDCP=5000)

VAD, VPI, VDI, VCI, VSD, FD, IVD, VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DVC = metrics


table = tabulate([['VAD', f"{VAD.mean():.4f} ({VAD.std():.4f})", 0.505],
                ['VPI', f"{VPI.mean():.4f} ({VPI.std():.4f})", 0.05],
                ['VDI', f"{VDI.mean():.4f} ({VDI.std():.4f})", 24.073],
                ['VCI', f"{VCI.mean():.4f} ({VCI.std():.4f})", 17962],
                ['VSD', f"{VSD.mean():.4f} ({VSD.std():.4f})", 0.021],
                ['FD', f"{FD.mean():.4f} ({FD.std():.4f})", 1.42],
                ['IVD', f"{IVD.mean():.4f} ({IVD.std():.4f})", 22*1e-4],
                ['VAD_ICP', f"{VAD_ICP.mean():.4f} ({VAD_ICP.std():.4f})", 0.21],
                ['VAD_DCP', f"{VAD_DCP.mean():.4f} ({VAD_DCP.std():.4f})", 0.17],
                ['VAD_ICP_DCP', f"{VAD_ICP_DCP.mean():.4f} ({VAD_ICP_DCP.std():.4f})", 0.29],
                ['FD_DVC', f"{FD_DVC.mean():.4f} ({FD_DVC.std():.4f})", 1.48]],
               headers=['Metric', 'Mean (std)', 'Literature'],
               tablefmt='simple_grid')
print(table)
print(f"Analyzed {len(VAD):d} trees.")

df = pd.DataFrame.from_dict({i:{'VAD':VAD[i], 'VPI':VPI[i], 'VDI':VDI[i],
                                'VPI':VPI[i], 'VCI':VCI[i], 'VSD':VSD[i],
                                'FD':FD[i], 'IVD':IVD[i], 'VAD_ICP':VAD_ICP[i],
                                'VAD_DCP':VAD_DCP[i], 'VAD_ICP_DCP':VAD_ICP_DCP[i],
                                'FD_DVC':FD_DVC[i],
                                'sim':i,
                                'Filename':listOfFiles[i]}
                             for i in range(len(listOfFiles))}, orient = 'index')

#df.to_csv(csvOutFile)
ys = []
for x in [['VAD', 0.505],
                ['VPI', 0.05],
                ['VDI', 24.073],
                ['VCI', 17962],
                ['VSD', 0.021],
                ['FD', 1.42],
                ['IVD', 22*1e-4],
                ['VAD_ICP', 0.21],
                ['VAD_DCP', 0.17],
                ['VAD_ICP_DCP', 0.29],
                ['FD_DVC', 1.48]]:
    
    try:
        df[x[0]] /= x[-1]
        ys.append(x[0])
    except KeyError:
        pass
        
fig = plt.figure(figsize=(6,3))
rename = {'VAD_ICP':'VAD (ICP)', 'VAD_DCP':'VAD (DCP)', 'FD_ICP_DCP':'FD (ICP+DCP)', 'VAD_ICP_DCP':'VAD (ICP+DCP)'}
df = df[ys].rename(columns=rename)
ys = [rename[y] if y in rename.keys() else y for y in ys]
plot_order = list(df[ys].mean().sort_values().index)
g = sns.boxplot(data=(df[ys]-1)*100, color='b', orient='horizontal', order=plot_order, boxprops={'facecolor':'white'})
#g = sns.barplot(data=(df[ys]-1)*100, color='b', orient='horizontal', order=plot_order)
#g.plot([0,0], g.get_ylim(), 'k')
#g.set_yticklabels(['VAD', 'VDI', 'VCI', 'VSD', 'FD', 'IVD', 'VAD (ICP)', 'VAD (DCP)','VAD (ICP+DCP)', 'FD (ICP+DCP)'])
g.yaxis.grid(True)
plt.tight_layout()
plt.show()
