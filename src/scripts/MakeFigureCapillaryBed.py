# coding: utf-8
from CCO2Graph import Smoothing
from vp_utils import VirtualRetinalVasculature, intersect
from DVC.dvc_utils import CapillaryBed
from SVC.CCO.cco_utils import AVGraph
import networkx as nx, matplotlib.pyplot as plt, matplotlib, numpy as np
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


xlims = (.1058, .17)
ylims = (0.120, 0.160)

G = AVGraph("../Results/DataForPaper/Baseline-Sims-nTerm-Fixed/AVTrees/sim_14_AV.cco") 
Smoothing(G, alpha=.2, nPoints=2)

G.remove_nodes_from([n for n,p in G.nodes.data('position') if np.linalg.norm(p[:2])>.2])

vessels = np.array([[G.nodes[n]['position'][:2] for n in e] for e in G.edges])
Gcap = CapillaryBed(nPoints=4000, radiusROI=0.25, homogeneous=False,
      radiusFAZ=np.linalg.norm(vessels, axis=2).min()/2)
_l = np.array([[Gcap.nodes[n]['pos'][:2] for n in e] for e in Gcap.edges])
linestyles=['dashed'if np.any(intersect(line[0, np.newaxis], line[1,np.newaxis], vessels[:,0,:], vessels[:,1,:])) else 'solid'for line in _l]

lines_AV = LineCollection(
    vessels,
    colors=['b' if G.nodes[n]['nodeType']=='cap' else 'r' if G.nodes[n]['nodeType']=='art' else 'g' for n,_ in G.edges],
    linewidths=np.array([1e3/2.54*G[u][v]['radius'] if G.nodes[u]['nodeType']!='art' else 1e3/2.54*G[u][v]['radius']*2  for u,v in G.edges]),
)
lines_cap = LineCollection(_l,
    linewidths=[5*1e3/2.54*Gcap[u][v]['radius'] for u,v in Gcap.edges],
    linestyles=linestyles
    )


fig = plt.figure(layout='constrained')
gs = GridSpec(3,2, figure=fig) 

ax.add_collection(lines_AV)
ax.add_collection(lines_cap)
ax.set_xlim(-.12, .12)
ax.set_ylim(-.12, .12)

ax.add_patch(Rectangle((xlims[0],ylims[0]), xlims[1]-xlims[0], ylims[1]-ylims[0], fill=False, linewidth=2, color='k', linestyle='--'))
plt.axis('off')
plt.savefig("../../../Writing/VP retinal vasculature/imgs/CapillaryBed_with_AV.svg", format='svg', bbox_inches='tight', pad_inches=0)
plt.savefig("../../../Writing/VP retinal vasculature/imgs/CapillaryBed_with_AV.tiff", format='tiff', bbox_inches='tight', pad_inches=0, dpi=400)
plt.close('all')

print("First plot, done.")

vessels = np.array(
    [[G.nodes[n]['position'][:2] for n in e] for e in G.edges
#        if G.nodes[e[0]]['position'][1]<0 and G.nodes[e[0]]['position'][0]>0
#        if xlims[0]<G.nodes[e[0]]['position'][0]<xlims[1] and ylims[0]<G.nodes[e[0]]['position'][1]<ylims[1]
        ])
   
lines_AV = LineCollection(
    vessels,
    colors=['b' if G.nodes[n]['nodeType']=='cap' else 'r' if G.nodes[n]['nodeType']=='art' else 'g' for n,_ in G.edges],
    linewidths=np.array([10*1e3/2.54*G[u][v]['radius'] if G.nodes[u]['nodeType']!='art' else 1e3/2.54*G[u][v]['radius']*2  for u,v in G.edges])
)

lines_cap = LineCollection(
    _l,
    linewidths=[30*1e3/2.54*Gcap[u][v]['radius'] for u,v in Gcap.edges],
    linestyles=linestyles
    )

fig, ax = plt.subplots()
#ax = fig.add_subplot(gs[0,1])
ax.add_collection(lines_AV)
ax.add_collection(lines_cap)
ax.set_xlim(*xlims)
ax.set_ylim(*ylims)

ax.margins(0.3)

plt.axis('off')
plt.savefig("../../../Writing/VP retinal vasculature/imgs/CapillaryBed_with_AV_ZoomedIn.svg", format='svg', bbox_inches='tight', pad_inches=0)
plt.savefig("../../../Writing/VP retinal vasculature/imgs/CapillaryBed_with_AV_ZoomedIn.tiff", format='tiff', bbox_inches='tight', pad_inches=0, dpi=400)
#plt.show();
plt.close('all')

print("Second plot, done.")

lines_AV = LineCollection(
    vessels,
    colors=['b' if G.nodes[n]['nodeType']=='cap' else 'r' if G.nodes[n]['nodeType']=='art' else 'g' for n,_ in G.edges],
    linewidths=np.array([1e3/2.54*G[u][v]['radius'] if G.nodes[u]['nodeType']!='art' else 1e3/2.54*G[u][v]['radius']*2  for u,v in G.edges])
)

Gcap.remove_edges_from((e for e,sty in zip(Gcap.edges, linestyles) if sty!='solid'))
_l = np.array([[Gcap.nodes[n]['pos'][:2] for n in e] for e in Gcap.edges])
lines_cap = LineCollection(_l,
    linewidths=[5*1e3/2.54*Gcap[u][v]['radius'] for u,v in Gcap.edges] 
    )


fig, ax = plt.subplots()
#ax = fig.add_subplot(gs[1:,:])
ax.add_collection(lines_AV)
ax.add_collection(lines_cap)
ax.set_xlim(-.12, .12)
ax.set_ylim(-.12, .12)

ax.margins(0.25)
plt.axis('off')
plt.savefig("../../../Writing/VP retinal vasculature/imgs/CapillaryBed_with_AV_IntersectionRemoved.svg", format='svg', bbox_inches='tight', pad_inches=0)
plt.savefig("../../../Writing/VP retinal vasculature/imgs/CapillaryBed_with_AV_IntersectionRemoved.tiff", format='tiff', bbox_inches='tight', pad_inches=0, dpi=400)
#plt.show()
plt.close('all')

print("Last plot, done.")

#fig.savefig("../../../Writing/VP retinal vasculature/imgs/Figure2.svg", format='svg')
#fig.savefig("../../../Writing/VP retinal vasculature/imgs/Figure2.tiff", format='tiff', dpi=400)
#plt.show()
plt.close('all')
