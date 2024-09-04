import os
import numpy as np
import matplotlib.pyplot as plt


from scipy.linalg import orthogonal_procrustes
from scipy.special import comb
from scipy.spatial.transform import Rotation as rotation
import scipy.interpolate as spi

class SSM_RetinalArcades(object):
    def __init__(self, trainingShapes:dict):
        self._subModels = {key:SSM_oneArcade(arcade) for key,arcade in trainingShapes.items()}

    def PlotVariances(self):
        fig, axs = plt.subplots(2,2, figsize=(8,6),sharex=True, sharey=True); axs=axs.ravel()
        for model, ax in zip(self._subModels.items(), axs):
            model[1].PlotVariance(ax=ax)
            ax.set_title(model[0])
            ax.set(xlabel=None, ylabel=None)
        fig.text(0.5, 0.04, 'Variance explained [%]', ha='center')
        fig.text(0.04, 0.5, 'Number of modes of variation', va='center', rotation='vertical')
        plt.show()

    def PointClouds(self, n=200):
        fig, axs = plt.subplots(2,2, figsize=(8,6),sharex=True, sharey=True); axs=axs.ravel()
        for model, ax in zip(self._subModels.items(), axs):
            model[1].PointClouds(ax=ax)
            ax.set_title(model[0])
            ax.set(xlabel=None, ylabel=None)
            ax.axis('off')
        plt.show()

    def Generate(self, m=5, c=1.5,
                 saveIn:str=None,
                 rCRA:float=0.008, # In cm, from Dorner et al. (2009)
                 k:float=1.11, # From Goldenberg et al. (2013)
                 span:float=2, # In cm, the length (straight line from OD to end of longest arcade)
                               # to rescale the shapes. See with Savita what should be a good length.
                 vCRA:float=6.3, # In cm/s, the flow velocity through the CRA
                 **kwargs # Additional parameters for the .cco file if saveIn is given
                ):
        """
        Generate a new shape using the statistical shape model.
        @param m: number of modes to use. Can be int or list[int]
                  with a number for each sub-models.
        @param c: Scaling of variance in the generated shapes. 
        @param saveIn: file to save the generated arcades.
        @param rCRA: radius of the CRA. Used if saveIn.
        @param k: ratio rCRV/rCRA. Used if saveIn.
        """
        if isinstance(m,int):
            arcades = {key:model.Generate(m=m,c=c) for key,model in self._subModels.items()}
        elif isinstance(m,list):
            arcades = {key:model.Generate(m=modes,c=c) for (key,model),modes
                       in zip(self._subModels.items(), m)}
        else:
            raise ValueError(f"m must be an int or a list of 4 ints. Instead got {type(m)}.")

        # Find a suitable point for the OD
        OD = np.mean([arcade[0] for arcade in arcades.values()],axis=0)
        # Smooth the arcades (Catmull-Rom should avoid self loops)
        smoothArcades = SmoothShapes([shape for shape in arcades.values()], method='Catmull-Rom', nPts=10)
        arcades = {key:smoothArcade for key, smoothArcade in zip(arcades.keys(), smoothArcades)}
        # Translate all the arcades to have the same starting point (OD)
        arcades = {key:arcade + (OD-arcade[0,:])[np.newaxis,:] for key, arcade in arcades.items()}

        if saveIn:
            self.SaveArcades(arcades, saveIn, rCRA, k, span, vCRA, **kwargs)
            return
        else:
            return arcades

    def SaveArcades(self, arcades,
                    fileName, # Base file name. _artery/_vein will be appended to differentiate arcades
                    radiusCRA,
                    k,
                    span,
                    vCRA,
                    **kwargs):

        # dp = kwargs.get('dp', 16.6*133.3)
        # refPressure = kwargs.get('refPressure', 60*133.3)
        
        def WriteTwoArcades(IA, SA, vesselType:str):
            if vesselType == 'artery':
                filename = fileName+"_artery.cco"
                radius, q = radiusCRA, vCRA*(radiusCRA**2)*np.pi
                refPressure = kwargs.get('pCRA', 52.8)
                #refPressure = (2./3.)*kwargs.get('MAP',62.62) - kwargs.get('IOP', 15) # Guidoboni 2014
                dp = refPressure - kwargs.get('capPressure', 23) # Takahashi 2009
            else:
                filename = fileName+"_vein.cco"
                radius, q = radiusCRA*k, vCRA*(radiusCRA**2)*np.pi  # vCRA*(radiusCRA**2)*np.pi*k*k
                refPressure = kwargs.get('pCRV', kwargs.get('CRVP', 14)) # Guidoboni 2014
                dp = abs(kwargs.get('capPressure',23)-refPressure) # Takahashi 2009                                


            refPressure *= 133.332*10
            dp *= 133.332 * 10
                
            with open(filename, 'w') as f:

                nVessels = IA.shape[0] + SA.shape[0]-1 # Each shapes nPts-1 vessels. Add one for the CRA  
                f.write("*Tree\n")
                treeInfo = f"{IA[0,0]} {IA[0,1]} -{kwargs.get('lengthCRA',kwargs.get('lengthCRV',0.1))} {q} {kwargs.get('psiFactor',9e-6)}"
                treeInfo+= f" {dp} 2 {refPressure} {IA.shape[0]+SA.shape[0]} {radius} {kwargs.get('tol', 1e-6)}\n\n"
                f.write(treeInfo)

                f.write("*Vessels\n")
                f.write(f"{nVessels}\n")
                branchingMode = 2
                vesselsInfo = [f"0 {IA[0,0]} {IA[0,1]} -{kwargs.get('lengthCRA',kwargs.get('lengthCRV',0.1))} {IA[0,0]} {IA[0,1]} 0.0 0.0 0.0 0.0 0.0 {branchingMode} {radius} {q} 0.0 0.0 0.0 0 0.0 0.0 -100"]
                branchingMode = kwargs.get('branchingMode', 1) # enum  BRANCHING_MODE {
                # NO_BRANCHING, RIGID_PARENT, DEFORMABLE_PARENT,
                # DISTAL_BRANCHING, ONLY_AT_PARENT_HOTSPOTS}
                vesselsInfo += [f"{i} {xprox[0]} {xprox[1]} 0.0 {xdist[0]} {xdist[1]} 0.0 0.0 0.0 0.0 0.0 {branchingMode} {radius} {q/2.0} 0.0 0.0 0.0 0 0.0 0.0 0"
                                for i,(xprox,xdist) in enumerate(zip(IA[:-1], IA[1:]), start=1)]
                IA_start, IA_end = int(vesselsInfo[1].split(' ')[0]), int(vesselsInfo[-1].split(' ')[0])
                vesselsInfo += [f"{i} {xprox[0]} {xprox[1]} 0.0 {xdist[0]} {xdist[1]} 0.0 0.0 0.0 0.0 0.0 {branchingMode} {radius} {q/2.0} 0.0 0.0 0.0 0 0.0 0.0 0"
                                for i,(xprox,xdist) in enumerate(zip(SA[:-1], SA[1:]), start=IA_end+1)]
                SA_start, SA_end = IA_end+1, int(vesselsInfo[-1].split(' ')[0])
                f.write('\n'.join(vesselsInfo))

                f.write('\n\n*Connectivity\n')
                vesselsConn = [f"0 -1 {IA_start} {IA_end+1}", f"{IA_start} 0 {IA_start+1}", f"{IA_end+1} 0 {IA_end+2}"]
                vesselsConn += [f"{i} {i-1} {i+1}" for i in range(IA_start+1, IA_end)] + [f"{IA_end} {IA_end-1}"]
                vesselsConn += [f"{i} {i-1} {i+1}" for i in range(SA_start+1, SA_end)] + [f"{SA_end} {SA_end-1}"]
                f.write('\n'.join(vesselsConn))
            return

        WriteTwoArcades(arcades['Inf. Art. arcade'], arcades['Sup. Art. arcade'], 'artery')
        WriteTwoArcades(arcades['Inf. Vei. arcade'], arcades['Sup. Vei. arcade'], 'vein')
        
        return 
            
            
    def Save(self, fileName:str):
        """
        Saves a model to a text file.
        """
        with open(fileName, 'w') as f:
            print("Writing arcades model to", fileName)
            f.write("# Four arcade models.")
            for key, model in self._subModels.items():
                f.write(f"\n# {key}\n")
                f.write(f"## eigenvectors\n")
                f.write('\n'.join([' '.join([x for x in e]) for e in model._evect.astype(str)]))
                f.write('\n##Eigenvalues\n')
                f.write(' '.join([x for x in model._eval.astype(str)]))

    ### TODO loading methods.
    # def Load(self, fileName:str):
    #     """
    #     Load a model from a text file.
    #     """
    #     self._subModels = {}
    #     with open(fileName, 'w') as f:
    #         print("Reading arcades model from", fileName)
    #         line = f.readline() # Header
    #         key = f.readline() # Model name
    #         self._subModels[key] = SSM_oneArcadel()


class SSM_oneArcade(object):
    def __init__(self, trainingPoints:list):
        '''
        Train the statistical shape model from the training points.

        ### Parameters
        trainingPoints: list
            List of numpy arrays giving the location of the points. 
        '''
        self._X = np.array([shape.ravel() for shape in trainingPoints]) # The data matrix
        self.k, self.n = self._X.shape # # of samples, # of points per shape
        self._cov = np.corrcoef(self._X, rowvar=False) # Correlation matrix
        self._eval, self._evect = np.linalg.eigh(self._cov) # Eigen pairs
        self._xbar = np.mean(self._X, axis=0).reshape((-1,2))
        
    def PlotVariance(self, ax=None):
        # if not np.all(self._eval>0):
        #     print("Not all eigenvalues are real positive.")
        if ax:
            plot=ax.plot
        else:
            plot=plt.plot
        variance = 100*np.cumsum(np.flip(abs(self._eval), axis=0))/np.sum(abs(self._eval))
        mask = [True] + [True if variance[i-1]<95 else False for i in range(1,variance.size)]
        plot(np.arange(len(self._eval))[mask]+1, variance[mask])
        plt.xlabel("Number of modes of variation.")
        plt.ylabel("Explained variance (%).")
        if ax:
            return
        plt.show()
        
    def Generate(self, b=None, m=5, c=3):
        """
        Generate a new shape given a coefficient vector b (b.shape=(m,)).
        m is the number of variation modes to use.
        c dictates the variation allowed. Larger c allows for more variation in
        the modes of variations. c can be a vector of size m to control
        variations in all modes.
        """
        P = np.real(self._evect[-m:]).T
        scaleFactor = 2.5/540.0 # cm/pixel, the rescaling factor estimated with the rule
                                # of thumb: 10 degree FOV = 5mm. DRIVE images are 50 degree
        if b:
            return ( self._xbar + (P.dot(b)).reshape((-1,2)) )*scaleFactor
        else:            
           # return self._xbar + (P.dot(np.random.normal(0.0, c*np.sqrt(np.real(self._eval[-m:]))))).reshape((-1,2))
            return ( self._xbar + (P.dot(np.random.normal(0.0, c*np.real(self._eval[-m:])))).reshape((-1,2)) )*scaleFactor 

    def PointClouds(self, ax=None, n = 200):
        colors = plt.cm.rainbow(np.linspace(0, 1, self._xbar.shape[0]))
        if not ax:
            fig, ax = plt.subplots()
        for shape in [self.Generate(c=3, m=2) for _ in range(n)]:
            ax.scatter(*shape.T, c=colors, alpha=0.2, s=1)
        ax.set_title("Probability distribution of generate shapes.")
        if not ax:
            plt.tight_layout()
            plt.show()

def bernstein_poly(i, n, t):
    """
     The Bernstein polynomial of n, i as a function of t
    """
    return comb(n, i) * ( t**(n-i) ) * (1 - t)**i

def bezier_curve(points, nTimes=1000):
    """
       Given a set of control points, return the
       bezier curve defined by the control points.
       points should be a list of lists, or list of tuples
       such as [ [1,1], 
                 [2,3], 
                 [4,5], ..[Xn, Yn] ]
        nTimes is the number of time steps, defaults to 1000
        See http://processingjs.nihongoresources.com/bezierinfo/
    """
    nPoints = len(points)
    xPoints = np.array([p[0] for p in points])
    yPoints = np.array([p[1] for p in points])
    t = np.linspace(0.0, 1.0, nTimes)
    polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])
    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)
    return xvals, yvals

def f(shape):
    """
    Returns the (unnormalized) curvature of the shape,
    i.e., the function f(t)=X''(t)Y'(t)-Y''(t)X'(t).
    @param shape: a nx2 numpy array. 
    """
    s = np.sum(np.diff(shape, axis=0)**2, axis=1)**.5 # Arc length
    J = np.divide(np.diff(shape, n=1, axis=0)[1:],s[1:,np.newaxis]) # The first derivative
    H = np.diff(shape, n=2, axis=0)/(s[1:, np.newaxis]**2) # The second derivative
    f = J[:,1]*H[:,0]-J[:,0]*H[:,1] # The curvature
    return f

def SplineInflectionPoints(spline, nPoints=100, tol=1e-6):
    xp,yp = spline[0].derivative(nu=1), spline[1].derivative(nu=1)
    xpp,ypp = spline[0].derivative(nu=2), spline[1].derivative(nu=2)
    print(xpp.k, ypp.k)
    t = np.linspace(0,1,nPoints,endpoint=True)[1:-1] # End points will be added regardless of inflection points
    # curv = yp(t)*xpp(t) - ypp(t)*xp(t)
    roots = [spi.sproot(xpp.tck), spi.sproot(ypp.tck)]
    # Find t_i such that x''=y''=0, i.e., inflection points
    ts = roots[0][np.linalg.norm(np.array(roots[0])[np.newaxis,:]-np.array(roots[1]), axis=0)<tol]
    # pts = np.where(np.diff(np.sign(curv)))[0] + 1
    pts = np.zeros((ts.size+2,2))
    pts[1:-1,:] = np.array([sp(ts) for sp in spline]).T
    pts[-1,:] = np.array([spline[0](1), spline[1](1)]).T
    pts[0,:] = np.array([spline[0](0), spline[1](0)]).T
    return pts

def InflectionPoints(shape):
    curv = f(shape) 
    pts = np.where(np.diff(np.sign(curv)))[0] + 1# Find where curv changes sign
    if pts[-1]!=shape.shape[0]-1:    
        pts = np.append(pts, shape.shape[0]-1)
    if pts[0]!=0:
        pts = np.insert(pts, 0,0)
    return pts

def UpSampleShape(shape, nmax):
    """
    Upsamples a shape by adding the mid-point of the 
    nmax-len(shape) longest segments to the point list.
    @param shape: a nx2 numpy array
    @param nmax: the desired number of landmarks (i.e.,
                 the maximum number of inflections points
                 in the training set).
    """
    upSampledShape = shape
    while len(upSampledShape)<nmax:
        # plt.scatter(*upSampledShape.T, marker='x')
        n = nmax-upSampledShape.shape[0] # Number of points needed
        # Find the n longest segments and add the midpoint to the training list
        newPts = [[i+1, (upSampledShape[i]+upSampledShape[i+1])/2.0] for _,i in sorted(zip(np.linalg.norm(upSampledShape[1:]-upSampledShape[:-1],axis=1), range(upSampledShape.shape[0]-1)))[:n]]
        upSampledShape = np.insert(upSampledShape, [i for i,_ in newPts], [p for _,p in newPts], axis=0)
    return upSampledShape

def GPA(Shapes:list, tol:float=1e-2, maxIter:int=3,
        centers:list=[], 
        with_scaling=True, with_reflection=False):
    '''
    Procrust alignment of a list of shapes.
    '''
    n = len(Shapes)
    if centers:
        # Translate all shapes to according to the given
        # center points.
        if len(centers)!=n:
            print(f"Expected {len(Shapes)} shape centers. Got {len(centers)}. Ignoring the centering.")
            centeredShapes = Shapes
        else:
            centeredShapes = [shape-center[np.newaxis,:] for shape,center in zip(shapes,centers)]
    else:
        centeredShapes = Shapes
        
    def OptimalScalingAndRotation(matricesTuple):
        A,B = matricesTuple
        A,B = A - A.mean(0), B - B.mean(0)
        normA, normB = np.linalg.norm(A), np.linalg.norm(B)
        A,B = A/normA, B/normB
        
        u,s,v = np.linalg.svd(A.T.dot(B))
        Q = v.dot(u.T) # Optimal roation
        if ~with_reflection and np.linalg.det(Q)<0: # Cancel reflections
            v[:,-1] *= -1
            s[-1] *= -1
            Q = v.dot(u.T)
        if with_scaling:
            traceQA_TB = s.sum()
            # Optimal scaling
            scale = traceQA_TB * normB/normA
            Q *= scale
    
        return Q

    refShape = sum(centeredShapes)/n
    i = 0
    
    while i<maxIter:
        dist=0
        optimalScalingAndRotation = list(map(OptimalScalingAndRotation, [(shape,refShape) for shape in Shapes]))
        centeredShapes = list(map(lambda x: x[0].dot(x[1]), # Apply the transformation (x[1]) to the array (x[0])
                              list(zip(centeredShapes,optimalScalingAndRotation))))
        # A similarity metric defined to be 1 when perfectly aligned and decrease when shapes spread apart
        similarity = 1./(1+sum([np.linalg.norm(refShape-shape) for shape in centeredShapes]))
        refShape = sum(centeredShapes)/n
        i+=1
        if similarity < tol:
            break
    print(f"Shapes aligned with similarity {1./(1+sum([np.linalg.norm(refShape-shape) for shape in centeredShapes]))}")
    return centeredShapes

def Kabsch_Umeyama(Shapes:list, maxIter:int=2, tol:float=1e-3):
    def kabsch_umeyama(A, B):
        """
        Align two shapes by rotating and scaling B.
        From https://zpl.fi/aligning-point-patterns-with-kabsch-umeyama-algorithm/
        """
        assert A.shape == B.shape
        n, m = A.shape

        EA = np.mean(A, axis=0)
        EB = np.mean(B, axis=0)
        VarA = np.mean(np.linalg.norm(A - EA, axis=1) ** 2)

        H = ((A - EA).T @ (B - EB)) / n
        U, D, VT = np.linalg.svd(H)
        d = np.sign(np.linalg.det(U) * np.linalg.det(VT))
        S = np.diag([1] * (m - 1) + [d])
        
        R = U @ S @ VT
        c = VarA / np.trace(np.diag(D) @ S)
        t = EA - c * R @ EB
        return R, c, t

    refShape = sum(Shapes)/len(Shapes) # Mean shape
    # refShape/= ((refShape-refShape.mean()).sum()**2).sum()
    copyShapes = Shapes
    i=0
    while i<maxIter:
        alignedShapes = [np.array([t+c*R@b for b in shape]) for shape, (R,c,t)
                         in zip(copyShapes, map(lambda x: kabsch_umeyama(refShape,x), copyShapes))]

        similarity = 1./(1+sum([np.linalg.norm(refShape-shape) for shape in alignedShapes]))

        i+=1
        if similarity < tol:
            break
        copyShapes = alignedShapes

    print(f"Shapes aligned with similarity {1./(1+sum([np.linalg.norm(refShape-shape) for shape in alignedShapes]))}")
    return alignedShapes        
    

def Kabsch(Shapes:list, scaling=True, error = 1e-1, maxIter = 2) -> list:
    """
    Align the shapes using rotations and scaling.
    """
    # Shapes is a list of k shapes made of (L, 2) L 2D landmarks
    centeredShapes = Shapes
    # centeredShapes = [shape[:-1]-shape[-1] for shape in Shapes]

    k = len(centeredShapes)
    newShapes = np.array(centeredShapes)
    newShapes = np.pad(newShapes, ((0,0), (0,0), (0,1)))

    i = 0
    print(f'Running generalized procruste analysis on the data.')
    while i < maxIter:
        # Is it ok to reassign refShape every time?    
        dist = 0
        refShape = np.mean(newShapes, axis=0)
        muB = refShape.mean(0)
        for i in range(k):
            # R, rssd = rotation.align_vectors(refShape, newShapes[i])
            muA = newShapes[i].mean(0)
            sigma2 = newShapes[i].var(axis=0)
            u,D,v = np.linalg.svd(((newShapes[i]-muA).T.dot(refShape-muB)).T)
            d = np.sign(np.linalg.det(u)*np.linalg.det(v.T))
            s = np.diag([1 for _ in range(len(D)-1)]+[d])
            c = sigma2/np.trace(np.diag(D).dot(s))
            newShapes[i] = muA + (newShapes[i]-muA).dot(c*v.dot(s).dot(u.T))
            dist += np.linalg.norm(newShapes[i]-refShape)
            # dist += rssd
            # newShapes[i] = R.apply(newShapes[i])
                

        print(f'\tDistance between mean shape and reference shape {dist=}.')
        i+=1
        if dist < error:
            print('\tError reached tolerance, terminating GPA.')
            print(f'\tTolerance achieved in {i} iterations. Stopping the alignement procedure.')
            break
    newShapes[0] = refShape
    return newShapes[...,:-1]

def CatmullRomSpline(P0, P1, P2, P3, a, nPoints):
      """
      P0, P1, P2, and P3 should be (x,y) point pairs that define the Catmull-Rom spline.
      nPoints is the number of points to include in this curve segment.
      """
      # Convert the points to numpy so that we can do array multiplication
      P0, P1, P2, P3 = map(np.array, [P0, P1, P2, P3])
      
      # Calculate t0 to t4
      alpha = a
      def tj(ti, Pi, Pj):
          xi, yi = Pi
          xj, yj = Pj
          return ( ( (xj-xi)**2 + (yj-yi)**2 )**0.5 )**alpha + ti
      
      t0 = 0
      t1 = tj(t0, P0, P1)
      t2 = tj(t1, P1, P2)
      t3 = tj(t2, P2, P3)
      
      # Only calculate points between P1 and P2
      t = np.linspace(t1,t2,nPoints)[1:]
      
      # Reshape so that we can multiply by the points P0 to P3
      # and get a point for each value of t.
      t = t.reshape(len(t),1)
      
      A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
      A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
      A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
      
      B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
      B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3
      
      C  = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
      return C
  
def CatmullRomChain(Points,alpha,nPoints):
    """
    Calculate Catmull Rom for a chain of points and return the combined curve.
    """
    P = Points
    x1=P[0][0]
    x2=P[1][0]
    y1=P[0][1]
    y2=P[1][1]
    x3=P[-2][0]
    x4=P[-1][0]
    y3=P[-2][1]
    y4=P[-1][1]
    dom=max(P[:,0])-min(P[:,0])
    rng=max(P[:,1])-min(P[:,1])
    pctdom=1
    pctdom=float(pctdom)/100
    prex=x1+np.sign(x1-x2)*dom*pctdom
    prey=(y1-y2)/(x1-x2)*(prex-x1)+y1
    endx=x4+np.sign(x4-x3)*dom*pctdom
    endy=(y4-y3)/(x4-x3)*(endx-x4)+y4
    P=list(P)
    P.insert(0,np.array([prex,prey]))
    P.append(np.array([endx,endy]))
    
    sz = len(P)
    
    # The curve C will contain an array of (x,y) points.
    C = [P[0]]
    for i in range(sz-3):
        c = CatmullRomSpline(P[i], P[i+1], P[i+2], P[i+3],alpha, nPoints)
        C.extend(c)
    return C
      
def SmoothShapes(Shapes:list, nPts:int=50, smoothness:float=0, k:int=3, method:str='BSpline'):
    if method=='BSpline':
        t = np.linspace(0,1,nPts, endpoint=True)[1:-1]
        taus = [np.cumsum(np.sum(np.diff(shape,axis=0)**2, axis=1)**.5) for shape in Shapes]
        taus = [np.insert(tau, 0,0)/tau[-1] for tau in taus]
        return [np.array([spi.BSpline(*spi.splrep(tau, shape[:,0], s=smoothness, k=k))(t),
                          spi.BSpline(*spi.splrep(tau, shape[:,1], s=smoothness, k=k))(t)]).T for shape, tau in zip(Shapes,taus)]
    elif method=='Catmull-Rom':
        points =  [np.array(CatmullRomChain(shape, alpha=1, nPoints=nPts))
                   for shape in Shapes]
        return points
    else:
        # Use Bezier curves
        return [np.array(bezier_curve(shape, nPts)).T for shape in Shapes]


def CreateShapeModels(dataFolder:str, saveAs:str=None) -> SSM_RetinalArcades:
    '''
    Assumes dataFolder is structured as:
    -dataFolder
        -Arteries
            -Superior_arcade
                -SegmentedArcadeFiles
            -Inferior_arcade
                -SegmentedArcadeFiles
        -Veins
            -Superior_arcade
                -SegmentedArcadeFiles
            -Inferior_arcade
                -SegmentedArcadeFiles
        -Fovea
            -FoveaLocationFiles
    
    ### Parameters
    dataFolder: str
        The path to the data.
    saveAs: str
        If provided, the trained model is saved in saveAs.

    ### Return
    The trained arcades model.
    '''

    ## Load the data
    ### Position of the fovea
    fovea_dir = dataFolder + "/Fovea/"
    foveas = {filename:np.loadtxt(fovea_dir+filename, skiprows=1)[1:] for filename in os.listdir(fovea_dir)}
    ### Superior arcades
    SAA_dir = dataFolder+"/Arteries/"+"Superior_arcade/"
    SAA_dataset = [np.loadtxt(SAA_dir+filename)-foveas.get(filename, 0.0) for filename in os.listdir(SAA_dir)]
    SVA_dir = dataFolder+"/Veins/"+"Superior_arcade/"
    SVA_dataset = [np.loadtxt(SVA_dir+filename)-foveas.get(filename, 0.0) for filename in os.listdir(SVA_dir)]
    
    ### Inferior arcades
    IAA_dir = dataFolder+"/Arteries/"+"Inferior_arcade/"
    IAA_dataset = [np.loadtxt(IAA_dir+filename)-foveas.get(filename, 0.0) for filename in os.listdir(IAA_dir)]
    IVA_dir = dataFolder+"/Veins/"+"Inferior_arcade/"
    IVA_dataset = [np.loadtxt(IVA_dir+filename)-foveas.get(filename, 0.0) for filename in os.listdir(IVA_dir)]

    ## Flip the shapes to correspond to left eyes
    ### TODO add the translation to the have the center of the fovea at (0,0)
    IAA_dataset = [shape if shape[0,0]<shape[-1,0] else np.array([-shape[:,0], shape[:,1]]).T
                                  for shape in IAA_dataset]
    IVA_dataset = [shape if shape[0,0]<shape[-1,0] else np.array([-shape[:,0], shape[:,1]]).T
                                  for shape in IVA_dataset]
    SAA_dataset = [shape if shape[0,0]<shape[-1,0] else np.array([-shape[:,0], shape[:,1]]).T
                                  for shape in SAA_dataset]
    SVA_dataset = [shape if shape[0,0]<shape[-1,0] else np.array([-shape[:,0], shape[:,1]]).T
                                  for shape in SVA_dataset]

    ## Smooth the shapes, find training shapes and align them
    smoothingMethod = 'BSpline'
    smoothness, k, nPts = 0, 3, 30 # Smoothness and degree of splines and number of sample points
    
    # IAA_smooth = SmoothShapes(IAA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)
    # IVA_smooth = SmoothShapes(IVA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)
    # SAA_smooth = SmoothShapes(SAA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)
    # SVA_smooth = SmoothShapes(SVA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)

    ### The shapes may not have the same number of points
    IAA_training = [shape[InflectionPoints(shape)] for shape in SmoothShapes(IAA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)]
    IVA_training = [shape[InflectionPoints(shape)] for shape in SmoothShapes(IVA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)]
    SAA_training = [shape[InflectionPoints(shape)] for shape in SmoothShapes(SAA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)]
    SVA_training = [shape[InflectionPoints(shape)] for shape in SmoothShapes(SVA_dataset, nPts=nPts, smoothness=smoothness, k=k, method=smoothingMethod)]

    ### Upsample shapes to have the same number of points
    Align = Kabsch_Umeyama # Choices are GPA, Kabsch, Kabsch_Umeyama
    n = len(IAA_training) # Should be the same, if not, have different n's
    IAA_training = Align([UpSampleShape(shape, nmax) for shape,nmax in zip(IAA_training, n*[max([len(shape) for shape in IAA_training])])])
    IVA_training = Align([UpSampleShape(shape, nmax) for shape,nmax in zip(IVA_training, n*[max([len(shape) for shape in IVA_training])])])
    SAA_training = Align([UpSampleShape(shape, nmax) for shape,nmax in zip(SAA_training, n*[max([len(shape) for shape in SAA_training])])])
    SVA_training = Align([UpSampleShape(shape, nmax) for shape,nmax in zip(SVA_training, n*[max([len(shape) for shape in SVA_training])])])

    ## Create the statistical shape model
    ssm = SSM_RetinalArcades({'Inf. Art. arcade':IAA_training, 'Inf. Vei. arcade':IVA_training,
                                                        'Sup. Art. arcade':SAA_training, 'Sup. Vei. arcade':SVA_training})
    #ssm.PlotVariances()
    #ssm.PointClouds(n=300)
    #arcades = ssm.Generate(m=2,c=4)
    #meanArcades = ssm.Generate(c=0, saveIn="test")
    #colors = plt.cm.rainbow(np.linspace(0,1,len(arcades)))
    # for (key, arcade),mean, color in zip(arcades.items(), meanArcades.values(), colors):
    #     plt.plot(*bezier_curve(arcade, 50), label=key, c=color)
    #     plt.plot(*mean.T, '--', c=color)
    #     plt.legend()
    #     plt.show()
        
    if saveAs:
        ssm.Save(saveAs)
    return ssm
