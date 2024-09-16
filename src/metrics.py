'''
Tools for analysis of a .cco file
'''
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import skimage.io as io
import numpy as np
import networkx as nx
import pandas as pd
from vp_utils import VirtualRetinalVasculature
import sys, os
import cv2
from scipy.ndimage import distance_transform_edt
from multiprocessing import Process, Manager
from tqdm import tqdm

class HiddenPrints:
    # Remove prints to terminal. 
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def StatisticsMultipleTrees(ccoFiles:list, FOV:float, domainAreaInCm:float|None=None, upToStage:int=100,
                            **kwargs):
    '''
    Compute vessel density mean and std 
    from all the trees in ccoFiles.
    '''
    VAD, VPI, VDI, VCI, VSD, FD, IVD = [], [], [], [], [], [], []
    VAD_ICP_DCP, VAD_ICP, VAD_DCP, FD_DCP = [], [], [], []

    for ccoFile in tqdm(ccoFiles):
        #print("Reading " + ccoFile, end='...\t')
        if ccoFile.endswith('.graph'):
            G = VirtualRetinalVasculature()
            G.Read(ccoFile)
        else:
            with HiddenPrints():
                try:
                    G = VirtualRetinalVasculature(ccoFile, **kwargs)
                except IndexError:
                    print("Failed to create the graph for " + ccoFile + ". Skipping it.")
                    continue

        FD.append(FractalDimension(G,FOV, **kwargs, imageFileName=G._ccoFile[:-4]+"_SVC_FD.png"))
        IVD.append(InterVesselDistance(G,FOV=FOV, **kwargs, imageFileName=G._ccoFile[:-4]+"_SVC_IVD.png"))

        if not domainAreaInCm:
            vessels, domainAreaInCm = SVPVessels(G, FOV=FOV)
        else:
            vessels, _ = SVPVessels(G, FOV=FOV)
            
        vesselDensity = VesselAreaDensity(vessels, domainAreaInCm)
        VAD.append(VesselAreaDensity(vessels, domainAreaInCm))            
        VPI.append(VesselPerimeterIndex(vessels, domainAreaInCm))
        VDI.append(VesselDiameterIndex(vessels, domainAreaInCm))
        VCI.append(VesselComplexityIndex(vessels, domainAreaInCm))
        VSD.append(VesselSkeletonDensity(vessels, domainAreaInCm))
        del vessels
	
        if kwargs.get('nPointsICP', 0):

            if not domainAreaInCm:
                vessels, domainAreaInCm = DVCVessels(G, FOV=FOV)
            else:
                vessels, _ = DVCVessels(G, FOV=FOV)
            VAD_ICP_DCP.append(VesselAreaDensity(vessels, domainAreaInCm))

            if not domainAreaInCm:
                vessels, domainAreaInCm = DCPVessels(G, FOV=FOV)
            else:
                vessels, _ = DCPVessels(G, FOV=FOV)
            VAD_DCP.append(VesselAreaDensity(vessels, domainAreaInCm))

            if not domainAreaInCm:
                vessels, domainAreaInCm = ICPVessels(G, FOV=FOV)
            else:
                vessels, _ = ICPVessels(G, FOV=FOV)
            VAD_ICP.append(VesselAreaDensity(vessels, domainAreaInCm))
            del vessels
            
            FD_DCP.append(FractalDimension(G, FOV, plexus=2, imageFileName=G._ccoFile[:-4]+"_DCP_FD.png"))
        else:
            VAD_ICP_DCP.append(0)
            VAD_DCP.append(0)
            VAD_ICP.append(0)
            FD_DCP.append(0)

    directory = os.path.dirname(G._ccoFile)
    os.system(f"mv {directory}/*.png {os.path.abspath(os.path.join(directory, os.pardir, 'Images'))}")

    return (np.array(VAD), np.array(VPI), np.array(VDI), np.array(VCI),
            np.array(VSD), np.array(FD), np.array(IVD),
            np.array(VAD_ICP_DCP), np.array(VAD_ICP), np.array(VAD_DCP),
            np.array(FD_DCP))

def SVPVessels(G:nx.DiGraph, FOV:float=0.6):
    '''
    Takes a graph of the vasculature (should be a subclass of nx.DiGraph)
    and return a DataFrame of its SVP vessels within the FOV.

    ### Parameters
    G: networkx.DiGraph like
        Graph of the vasculature.
    FOV: float
        The diameter of the region of interest (roi; a square), in cm.
    '''
    fov2 = FOV/2.
    subG = G.subgraph(n for n,p in G.nodes.data('position') if np.linalg.norm(p[:2], np.inf)<fov2 and G.nodes[n]['plexus']==0)
    vessels = [d for _,_,d in subG.edges.data()]
    domainAreaInCm = FOV**2 - np.pi*min(np.linalg.norm(p[:2]) for n,p in subG.nodes.data('position'))**2
    # vessels = [
    #     d for n1,n2,d in G.edges.data() if
    #     G.nodes[n1]['plexus']==G.nodes[n2]['plexus']==0  # Edge in SVP
    #     and np.linalg.norm(G.nodes[n1]['position'], np.inf)<fov2 # Node in roi
    #     and np.linalg.norm(G.nodes[n2]['position'], np.inf)<fov2 # Node in roi
    # ]
    # domainAreaInCm = FOV**2 - np.pi*(min(np.linalg.norm([p[:2] for n,p in G.nodes.data('position')])))**2
    # #print(f"Kept {len(vessels)}/{G.number_of_edges()} SVP vessels in the FOV.")
    return pd.DataFrame(vessels), domainAreaInCm

def ICPVessels(G:nx.DiGraph, FOV:float=0.6):
    '''
    Takes a graph of the vasculature (should be a subclass of nx.DiGraph)
    and return a DataFrame of its DCP vessels within the FOV.

    ### Parameters
    G: networkx.DiGraph like
        Graph of the vasculature.
    FOV: float
        The diameter of the region of interest (roi; a square), in cm.
    '''

    fov2 = FOV/2. 
    subG = G.subgraph(n for n,p in G.nodes.data('position') if np.linalg.norm(p[:2], np.inf)<fov2 and G.nodes[n]['plexus']==1)
    vessels = [d for _,_,d in subG.edges.data()]
    domainAreaInCm = FOV**2 - np.pi*min(np.linalg.norm(p[:2]) for n,p in subG.nodes.data('position'))**2
    # vessels = [
    #     d for n1,n2,d in G.edges.data() if
    #     G.nodes[n1]['plexus']==G.nodes[n2]['plexus']==1  # Edge in ICP
    #     and np.linalg.norm(G.nodes[n1]['position'], np.inf)<fov2 # Node in roi
    #     and np.linalg.norm(G.nodes[n2]['position'], np.inf)<fov2 # Node in roi
    # ]
    # domainAreaInCm = FOV**2 - np.pi*(min(np.linalg.norm([p[:2] for n,p in G.nodes.data('position')])))**2
    return pd.DataFrame(vessels), domainAreaInCm

def DCPVessels(G:nx.DiGraph, FOV:float=0.6):
    '''
    Takes a graph of the vasculature (should be a subclass of nx.DiGraph)
    and return a DataFrame of its DCP vessels within the FOV.

    ### Parameters
    G: networkx.DiGraph like
        Graph of the vasculature.
    FOV: float
        The diameter of the region of interest (roi; a square), in cm.
    '''

    fov2 = FOV/2.
    subG = G.subgraph(n for n,p in G.nodes.data('position') if np.linalg.norm(p[:2], np.inf)<fov2 and G.nodes[n]['plexus']==2)
    vessels = [d for _,_,d in subG.edges.data()]
    domainAreaInCm = FOV**2 - np.pi*min(np.linalg.norm(p[:2]) for n,p in subG.nodes.data('position'))**2
    # vessels = [
    #     d for n1,n2,d in G.edges.data() if
    #     G.nodes[n1]['plexus']==G.nodes[n2]['plexus']==2  # Edge in DCP
    #     and np.linalg.norm(G.nodes[n1]['position'], np.inf)<fov2 # Node in roi
    #     and np.linalg.norm(G.nodes[n2]['position'], np.inf)<fov2 # Node in roi
    # ]
    # print(f"Kept {len(vessels)}/{G.number_of_edges()} DCP vessels in the FOV.")
    return pd.DataFrame(vessels), domainAreaInCm

def DVCVessels(G:nx.DiGraph, FOV:float=0.6):
    '''
    Takes a graph of the vasculature (should be a subclass of nx.DiGraph)
    and return a DataFrame of its DCV vessels within the FOV.

    ### Parameters
    G: networkx.DiGraph like
        Graph of the vasculature.
    FOV: float
        The diameter of the region of interest (roi; a square), in cm.
    '''
    fov2 = FOV/2.
    subG = G.subgraph(n for n,p in G.nodes.data('position') if np.linalg.norm(p[:2], np.inf)<fov2 and (G.nodes[n]['plexus']==2 or G.nodes[n]['plexus']==1))
    vessels = [d for _,_,d in subG.edges.data()]
    domainAreaInCm = FOV**2 - np.pi*min(np.linalg.norm(p[:2]) for n,p in subG.nodes.data('position'))**2

    # vessels = [
    #     d for n1,n2,d in G.edges.data() if
    #     (G.nodes[n1]['plexus']==G.nodes[n2]['plexus']==1  # Edge in ICP
    #      or G.nodes[n1]['plexus']==G.nodes[n2]['plexus']==2) # edge in DCP
    #     and np.linalg.norm(G.nodes[n1]['position'], np.inf)<fov2 # Node in roi
    #     and np.linalg.norm(G.nodes[n2]['position'], np.inf)<fov2 # Node in roi
    # ]
    # print(f"Kept {len(vessels)}/{G.number_of_edges()} DVC vessels in the FOV.")
    return pd.DataFrame(vessels), domainAreaInCm

'''
Five-index quantitative metrics of foveal microvasculature, proposed by Chu et al. 2016.
'''
def VesselAreaDensity(vessels:pd.DataFrame, FOVArea:float=0.036):
    VAD = np.multiply(vessels.radius, vessels.length).sum()*2
    return VAD/(FOVArea)

def VesselPerimeterIndex(vessels:pd.DataFrame, FOVArea:float=0.6):
    VPI = vessels.length.sum()*2
    return (VPI/(FOVArea**2))*1e-4 # Converts from cm^-1 to micron^-1

def VesselDiameterIndex(vessels:pd.DataFrame, FOVArea:float=0.6):
    totalLength = vessels.length.sum()
    totalArea   = 2.0*np.multiply(vessels.radius, vessels.length).sum()
    return totalArea/totalLength * 1e4  # Convert from cm to micron

def VesselComplexityIndex(vessels:pd.DataFrame, FOVArea:float=0.6):
    totalPerimeter = vessels.length.sum()*2
    totalArea      = np.multiply(vessels.radius, vessels.length).sum()*2
    return (totalPerimeter**2/(4*np.pi*totalArea))/1.5711 # Normalized as in Chu 2016

def VesselSkeletonDensity(vessels:pd.DataFrame, FOVArea:float=0.6):
    return (vessels.length.sum()/(FOVArea**2))*1e-4

def ReadCCO(ccoFileName, FOV=0.6):
    # FOV is the side length of the
    # square centered at (0,0).
    # Vessels outside the square
    # are not loaded
    with open(ccoFileName, 'r') as f:
        # Unused *Tree information
        row = f.readline()
        row = f.readline()
        row = f.readline()

        # Vessels information
        row = f.readline()
        # print('Reading', row.strip(), '...')
        nVessels = int(f.readline())
        print("Reading", ccoFileName, "with", nVessels, "vessels in the tree (including root).")
        # print(nVessels, "vessels in the tree.")

        vessels = {}
        maxDist = FOV/2.0
        for i in range(nVessels):
            row = (f.readline()).split() # Split all columns in a list
            Id, xProx, xDist, r, q, stage = int(row[0]), [float(xi) for xi in row[1:4]], [float(xi) for xi in row[4:7]], float(row[12]), float(row[10]), int(row[-1])
            l = sum([(a-b)**2 for a,b in zip(xProx, xDist)])**.5

            if all([x<=maxDist for x in xDist]) or all([x<=maxDist for x in xProx]):
                vessels[Id] = [r, l, q, stage]

        row = f.readline()
        row = f.readline()
        # print('Reading', row.strip())
        connectivity = {}
        for i in range(nVessels):
            row = (f.readline().split())
            if (len(row)==2):
                Id, parentId = int(row[0]), int(row[1])
                connectivity[Id] = [parentId, []]
            else:
                Id, parentId, children = int(row[0]), int(row[1]), [int(i) for i in row[2:]]
                connectivity[Id] = [parentId, children]

        return vessels, connectivity

''' Image analysis.'''

        
def CreatePNG(G:VirtualRetinalVasculature, plexus:int=0, FOV:float=0.3,
              **kwargs) -> np.ndarray:
    '''
    Takes a graph of the vasculature (should be a subclass of nx.DiGraph)
    and return a DataFrame of its SVP vessels within the FOV.

    ### Parameters
    G: networkx.DiGraph like
        Graph of the vasculature.
    plexus: int
        Which plexus to include. 0,1,2 for SVP, ICP, DCP
    FOV: float
        The diameter of the region of interest (roi; a square), in cm.

    ### Returns
    im: np.ndarray (Lx,Ly,1)
        An 'artificial OCTA' where vessels are white pixels on a black background. 
    '''
    fov2 = FOV/2.
    if isinstance(plexus, int):
        vessels = [(n1,n2,r)
                   for n1,n2,r in G.edges.data('radius') if
                   G.nodes[n1]['plexus']==G.nodes[n2]['plexus']==plexus  # Edge in plexus
                   and np.linalg.norm(G.nodes[n1]['position'],np.inf)<fov2 # Node in roi
                   and np.linalg.norm(G.nodes[n2]['position'],np.inf)<fov2 # Node in roi
                   ]
    else:
        vessels = [(n1,n2,r)
                   for n1,n2,r in G.edges.data('radius') if
                   G.nodes[n1]['plexus'] in plexus
                   and G.nodes[n2]['plexus'] in plexus
                   and np.linalg.norm(G.nodes[n1]['position'],np.inf)<fov2 # Node in roi
                   and np.linalg.norm(G.nodes[n2]['position'],np.inf)<fov2 # Node in roi
                   ]
        
    #print(f"Kept {len(vessels)}/{G.number_of_edges()} SVP vessels in the FOV.")

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (8,8))) # Fig size in inches

    s = fig.get_size_inches()[0]/(0.394*FOV) # The scaling of FOV into the figure
    radiusInCmToPoints = kwargs.get('radiusInCmToPoints',
                                    lambda r: kwargs.get('radiusScaling',0.05)*r*2.0*s/0.0139) # Return the line width adapted for the figure size

    lines = LineCollection([[G.nodes[n1]['position'][:2], G.nodes[n2]['position'][:2]] for n1,n2,d in vessels],
                           colors='white',
                           linewidths=[radiusInCmToPoints(r) for n1,n2,r in vessels])

    ax.add_collection(lines)
    ax.set_xlim(-fov2, fov2); ax.set_ylim(-fov2,fov2);
    plt.axis('off')
    imgFile = kwargs.get('imageFileName', G._ccoFile[:-4]+".png")
    fig.savefig(imgFile, bbox_inches='tight', pad_inches=0, dpi=kwargs.get('dpi',500), transparent=False)
    plt.close(fig)
    plt.style.use('default') # Restaure default style (white background)
    im = io.imread(imgFile, as_gray=True)
    return im

def boxCount(img:np.ndarray, k:int):

     S = np.add.reduceat(
         np.add.reduceat(img, np.arange(0, img.shape[0], k), axis=0),
         np.arange(0, img.shape[1], k), axis=1)

     # We count non-empty (0) and non-full boxes (k*k)
     return len(np.where((S > 0) & (S < k*k))[0])


def FractalDimension(G:VirtualRetinalVasculature, FOV:float=0.3,
                     method='histogramdd', **kwargs):

    # try:
    #     im = io.imread(kwargs.get('imageFile', None), as_gray=True)
    # except FileNotFoundError:
    im = CreatePNG(G, FOV=FOV, **kwargs)            

    Ly,Lx = im.shape[:2]
    Ns = []

    if method=='histogramdd':
        scales = np.logspace(0.01,5, num=20, endpoint=False, base=2)
        pixels = np.transpose((im>0).nonzero())
        for scale in scales:
            H, edges = np.histogramdd(pixels, bins=(np.arange(0,Lx, scale), np.arange(0,Ly,scale)))
            Ns.append(np.sum(H>0))

    else:
        p = min(Ly,Lx)
        n = 2**np.floor(np.log(p)/np.log(2))
        n = int(np.log(n)/np.log(2))
        scales = 2**np.arange(n, 1, -1)
        for scale in scales:
            Ns.append(boxCount(im, scale))

    coeffs = np.polyfit(np.log(scales), np.log(Ns), 1)
    return -coeffs[0]


def IntercapillaryDistanceMap(G:VirtualRetinalVasculature,
                              FOV:float=0.4, **kwargs):
    # try:
    #     im = io.imread(kwargs.get('imageFile', None), as_gray=True)
    # except FileNotFoundError:
    im = CreatePNG(G, FOV=FOV, **kwargs)

    # # Erode
    kernel = np.ones((2,2), np.uint8)
    im = cv2.erode(im, kernel, iterations=2)

    # Compute distance map
    dHeight, dWidth = FOV/im.shape[0], FOV/im.shape[1] # Size of a pixel in cm
    D = distance_transform_edt(im==0, return_distances=True, sampling=[dHeight,dWidth])
    
    if kwargs.get('plot', False):
        plt.subplot(121)
        plt.imshow(im)
        
        plt.subplot(122)
        plt.imshow(D, cmap='hot')
        plt.colorbar()
        plt.show()

    return D


def InterVesselDistance(G:VirtualRetinalVasculature,
                        FOV:float=0.3, **kwargs):
    
    D = IntercapillaryDistanceMap(G, FOV, **kwargs)
    h,w = D.shape

    center = (int(w/2), int(h/2))
    rFAZ = kwargs.get('rFAZ', 0.016) # In cm
    radius = rFAZ*min(h,w)/FOV
    
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    ## Remove FAZ and limit to a circular FOV
    #mask = (dist_from_center >= radius) & (dist_from_center <= min(center[0], center[1], w-center[0], h-center[1]))
    ## Limit to a circular FOV
    mask = (dist_from_center <= min(center[0], center[1], w-center[0], h-center[1]))
    D[~mask] = 0
    # plt.imshow(D)
    # plt.show()
    
    return D.sum()/mask.sum()
