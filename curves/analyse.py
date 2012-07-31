# Curves force analysis software
# JWA (c)2009

# Analyse force curves
import numpy as np
import scipy.signal.signaltools as signaltools
import mlpy
import stats
import time

def find_peaks_wavelet(x,y,gap_thresh=5,window=10,ridge_thresh=10,snr_thresh=0.0,peak_thresh=0.0,scaleCutoff=(4.0,32.0),order=4,verbose=False):
    """ Find peaks by identifying ridges in CWT. """
    # Calculate CWT
    #dt=1
    dt=abs(x[1]-x[0])
    dj=0.1
    d,s=mlpy.cwt(y,dt,dj,wf='dog',p=order,extmethod='reflection') # p=4 good for peaks, p=1 good for steps
    assert d.shape[1]==x.size
    cutoff_i=(np.argmin(abs(s-scaleCutoff[0]))+1,np.argmin(abs(s-scaleCutoff[1]))+1)
    cwt=np.absolute(d[cutoff_i[0]:cutoff_i[1],:])
    scales=s[cutoff_i[0]:cutoff_i[1]]

    # Find peaks by finding ridges in CWT and subjecting to thresholds
    scaleRange=(max(scaleCutoff[0],scaleCutoff[1]*0.1),scaleCutoff[1])
    ridges,lMax=find_ridges(scales,cwt,gap_thresh,window)
    pk=[]
    for r in ridges:
        # Ridge length threshold
        if len(r)>=ridge_thresh:
            pkStrength,pkNoiseInd=peak_strength(r,scales,cwt,scaleRange)
            pkInd=r[-1][1]
            #pkInd=pkNoiseInd
            if pkStrength>peak_thresh:
                noise=local_noise(pkNoiseInd,y,window)
                snr=pkStrength/noise
                if verbose:
                    print "Ridge:",[(scales[i],x[j]) for i,j in r]
                    print "Peak strength:",pkStrength
                    print "SNR:",snr
                if snr>snr_thresh:
                    # Peak tuple: index, x, mean y over wind, max y over wind,
                    #             strength, snr
                    ywindow=y[stencil(pkInd,window,y.size)]
                    pk.append((pkInd,x[pkInd],np.mean(ywindow),np.max(ywindow),
                               pkStrength,snr))
    return pk,scales,cwt,ridges,lMax


def find_ridges(s,da,gap_threshold,window=10):
    """ Find ridges by linking local maxima at different scales. """
    # Find local maxima at largest scale
    minWindow=5
    dist_threshold=max(minWindow,2*window)
    #window=max(2*int(s[-1])+1,minWindow)
    lMax=local_maxima(da[-1,:],window)
    lMaxList=[lMax]
    ridges=[[(-1,x)] for x in lMax]
    ridgeLines=[] # completed ridges

    # Set gap numbers to zero
    gaps=[0]*lMax.size

    # Go backwards down the scales
    for i in range(-2,-s.size-1,-1):
        # Find local maxima at this scale
        #window=max(int(s[i])+1,minWindow)
        lMax=local_maxima(da[i,:],window)
        lMaxList.append(lMax)

        # Longest ridges get priority, so we sort lists (ridges and gaps together)
        zipped=zip(ridges,gaps)
        zipped.sort(key=lambda x: len(x[0]),reverse=True)
        unzipped=zip(*zipped)
        ridges=list(unzipped[0])
        gaps=list(unzipped[1])

        for j,r in enumerate(ridges):
            if gaps[j]<gap_threshold:
                # Find nearest local maxima of last maximum in ridge
                if lMax.size>0:
                    nearest=np.argmin(np.abs(lMax-r[-1][1]))
                    if np.abs(lMax[nearest]-r[-1][1])<dist_threshold:
                        gaps[j]=0
                        ridges[j].append((i,lMax[nearest]))
                        lMax=np.delete(lMax,nearest) # remove linked maximum from available list
                    else:
                        gaps[j]+=1
                else:
                    gaps[j]+=1
            else:
                # Save finished ridges, and remove from list
                if len(ridges[j])>1:
                    ridgeLines.append(ridges[j])
                del ridges[j]
                del gaps[j]
        # Initiate unlinked maxima as new ridges
        ridges.extend([[(i,x)] for x in lMax])
        gaps.extend([0]*lMax.size)

    # Add any remaining ridges
    ridgeLines.extend(filter(lambda r: len(r)>1, ridges))

    # Make local maximum list have a logical order
    lMaxList.reverse()
    return ridgeLines,lMaxList
    
def local_maxima(y,window=10):
    """ Find local maxima using a sliding window. """
    
    # Extend to multiple of window
    extendBy=window-np.mod(y.size,window)
    y2=np.append(y,np.repeat(y[-1],extendBy))
    z=np.reshape(y2,(-1,window))

    # Index of max of each window
    zArgMax=np.argmax(z,axis=1)

    # Ignore boundary maximums
    selInd=np.where((zArgMax!=0) & (zArgMax!=window-1))[0]
    localMaxima=selInd*window + zArgMax[selInd]

    # Shift by window/2 and repeat
    shift=window/2
    z=np.roll(z,shift)
    
    # Index of max of each window
    zArgMax=np.argmax(z,axis=1)
    # Ignore boundary maximums
    selInd=np.where((zArgMax!=0) & (zArgMax!=window-1))[0]
    localMaxima=np.append(localMaxima,selInd*window + zArgMax[selInd] - shift)
    
    # Sort the two additions to localMaxima
    localMaxima=np.sort(localMaxima)

    # If local maxima  are closer than window, kill smaller of two
    selInd=np.where(np.diff(localMaxima) < window)[0]
    if (selInd.size>0):
        selInd1=selInd
        selInd2=selInd+1

        yDiff=y[selInd1]-y[selInd2]

        localMaxima=np.delete(localMaxima,selInd1[np.where(yDiff<=0)[0]])
        localMaxima=np.delete(localMaxima,selInd2[np.where(yDiff>0)[0]])

    # Delete any maxima arising from the extension or edge
    localMaxima=np.delete(localMaxima,np.where((localMaxima>=y.size-1) | (localMaxima==0))[0])

    return localMaxima

def peak_strength(ridge,scale,cwt,scaleRange):
    """ Calculate peak strength for a given ridge. """
    assert scaleRange[1]>scaleRange[0]
    pkStr=0.0
    pkInd=0
    #pkInd=ridge[0][1]
    for s,i in ridge:
        if scale[s]>=scaleRange[0] and scale[s]<=scaleRange[1] and cwt[s,i]>pkStr:
            pkStr=cwt[s,i]
            pkInd=i
    return pkStr,pkInd

def local_noise(i,y,window=10):
    """ Calculate the signal-to-noise ratio. """
    #n=cwt.shape[1]
    n=y.size
    #noiseSamples=cwt[0,stencil(i,window/2,n)]
    noiseSamples=y[stencil(i,window/2,n)]
    #noise=np.mean(noiseSamples)
    noise=np.std(noiseSamples)
    #noiseSamples=np.sort(noiseSamples)
    # define local noise as 95th quantile
    #noise=noiseSamples[np.floor(noiseSamples.size*0.95)]
    return noise

def flatten(l):
    """ Flatten list of lists. """
    for el in l:
        if isinstance(el, list):
            for sub in el:
                yield sub
        else:
            yield el




# Input: baseline subtracted data
# Output: list of x-coord,amplitude,width,x-index,left_index,right_index
def find_peaks(x,y,width_thresh=50.0,width_height=.7,dy_thresh=0.1,amp_thresh=-1,invert=False,verbose=False,window=10):
    """Finds peaks for the x, y data given. Parameters are derivative threshold,
    width at % height threshold, amplitude threshold. Returns list of coordinate
    pairs."""
    # Invert?
    if invert:
        y = -y
        if verbose:
            print "Inverting"
    # Smooth first
    #y = smooth(x,y,window,method="median")
    
    # Find peaks as roots of 1st derivative
    sdy = derivative(x,y,x[1]-x[0])

    peaks = []
    derivLeft=0
    derivRight=0
    startPoint=0

    while startPoint<y.size-1:
        # Find root of derivative
        derivLeft = find_sign_change(sdy[startPoint:])
        if derivLeft == -1:
            break
        derivLeft += startPoint
        if derivLeft >= sdy.size-2:
            break
        # Find where crosses back
        derivRight = find_sign_change(sdy[derivLeft+1:])
        if derivRight == -1:
            break
        derivRight += derivLeft+1
        if derivRight >= sdy.size-1:
            break
        assert derivRight>derivLeft

        # Check if most negative derivative exceeds thresh
        #(-ve since looking for +ve peaks)
        max_deriv=max(sdy[derivLeft:derivRight])
        # Only look at intervals larger than 4 points and with large derivatives
        if (derivRight-derivLeft)>4 and max_deriv > dy_thresh:
            if verbose:
                print "Potential: {0:.3}[{1}]-{2:.3}[{3}] {4:.3}".format(x[derivLeft],derivLeft,x[derivRight],derivRight,max(y[derivLeft:derivRight]))
            # Centre of peak at max
            #pk_centre = derivLeft+argmax(y[derivLeft:derivRight])+1

            # Fit parabola over three points around peak
            # to find true height estimate
            #P = polyfit(x[pk_centre-2:pk_centre+2],y[pk_centre-2:pk_centre+2],2)
            
            # Subtract mean for conditioning
            meanx=np.mean(x[derivLeft:derivRight])
            xp=x[derivLeft:derivRight]-meanx
            # Subtract off linear trend
            #meany=mean(y[derivLeft:derivRight])
            L=np.polyfit(xp,y[derivLeft:derivRight],1)
            yp=y[derivLeft:derivRight]-np.polyval(L,xp)
            #print y[derivLeft:derivRight],yp
            P = np.polyfit(xp,yp,2)
            dProots=np.real(np.roots(polyder(P)))
            assert len(dProots)==1
            pk_x_uncorr = dProots[0]
            pk_height_uncorr = np.polyval(P,pk_x_uncorr)
            pk_centre = derivLeft+np.argmin(xp-pk_x_uncorr)
            pk_x = pk_x_uncorr+meanx # Add mean back in
            pk_height=pk_height_uncorr + np.polyval(L,pk_x_uncorr) # Add back correction for trend subtraction


            #pk_height=max(y[derivLeft:derivRight])
            #pk_x=x[pk_centre]
            if verbose:
                print "Centre:",pk_x,"Height:",pk_height

            # Test amplitude threshold if required and make sure peak inside dataset
            if pk_x>min(x) and pk_x<max(x):
                if (amp_thresh <= 0 or pk_height > amp_thresh):

                    # Height to test width at
                    test_height = width_height * pk_height_uncorr
            
                    # Find extent of peak at test height
                    #pk_left = find_last_less(y[:pk_centre],test_height)
                    #if pk_left == -1:
                    #    pk_left = 0
                    #pk_right = find_first_less(y[pk_centre:],test_height)
                    #if pk_right == -1:
                    #    pk_right = y.size-1
                    #pk_right += pk_centre
                    #if pk_right > y.size-1:
                    #    pk_right = y.size-1

                    # Test: Use parabolic fit around peak to estimate width
                    assert len(P)==3
                    P[2]=P[2]-test_height # adjust constant coeff
                    pk_bounds=np.real(np.roots(P))
                    assert len(pk_bounds)==2
                    pk_right=max(pk_bounds) + meanx
                    pk_left=min(pk_bounds) + meanx

                    #pk_width = x[pk_right]-x[pk_left]
                    pk_width=pk_right-pk_left
                    if verbose:
                        print "Width:",pk_width
                    #print derivLeft, derivRight, pk_centre, pk_left, pk_right, x[pk_left],x[pk_right],pk_width
                    # Check peak width
                    # Not sure which relation is most useful
                    #if pk_width < width_thresh:
                    if pk_width > width_thresh:
                        if verbose:
                            print "Peak found at",x[pk_centre],"height",y[pk_centre],"width",pk_width,"derivative",max_deriv
                        peaks.append([pk_x,-pk_height if invert else pk_height,pk_width,pk_centre,pk_left,pk_right,max_deriv])
                    elif verbose:
                        print "Rejected on width"
                elif verbose:
                    print "Rejected on amplitude"
            elif verbose:
                print "Rejected because peak centre outside data range"
        # Move on
        startPoint = derivRight+1

    return peaks # List of tuples

# Find steps in force curve
# see Canny, IEEE Trans. Patt. Anal. Mach. Intel. 8:679-698 (1986)
# w: convolution stencil width
# stddevs: number of Gaussian standard deviations to convolve in stencil
def find_steps(x,y,w=201,cw=4,stddevs=4,peak_width=1,peak_height=.7,peak_amp=-1,
           invert=False,step_thresh=None):
    #h = abs(x[1]-x[0])
    # make stencil odd
    if w % 2 == 0:
        w += 1
    ww = w
    # Build Gaussian derivative
    # Number of +/-standard deviations over window 
    gx = np.linspace(-stddevs,stddevs,w)
    v = 1
    gy = -np.sqrt(2)*gx*np.exp(-.5*gx**2/v)/v**3
    i = cw

    # Edge response
    ey = np.empty(x.size)
    while i<x.size:
        # Stencil for convolution
        ind = stencil(i,w//2,y.size)
        # Convolute. Converts steps to peaks
        con = signaltools.convolve(y[ind],gy,mode='same')
        # Trim down to w
        if i+cw+1>x.size:
          ey[i-cw:x.size] = con[w/2-cw:w/2+x.size-i]
        else:
          ey[i-cw:i+cw+1] = con[w/2-cw:w/2+cw+1]
        i+=cw

    # Peak find to get steps
    #dey = derivative(x,ey,x[1]-x[0])
    peaks = find_peaks(x,ey,width_height=peak_height,width_thresh=peak_width,amp_thresh=peak_amp,invert=invert)
    # peaks is a list of lists (Z, amp, width, peak_centre, peak_left, peak_right)
    # convert to list of steps (Z, halfamp, height, width, step_centre)
    steps = []
    for p in peaks:
        pk_x = p[0]
        pk_width = p[2]
        pk_centre = p[3]
        pk_left = p[4]
        pk_right = p[5]
        #step_halfamp = mean(y[pk_centre-1:pk_centre+1])
        # Take mean of 5 points at each side to determine step height
        #step_heightLeft = mean(y[max(0,pk_left-5):pk_left])
        #step_heightRight = mean(y[pk_right:min(pk_right+5,y.size-1)])
        #step_height = step_heightRight - step_heightLeft 
        step_height = y[pk_right]-y[pk_left]
        step_halfamp=0.5*(y[pk_right]+y[pk_left])
        if step_thresh == None or step_height > step_thresh:
            s = []
            s.append(pk_x)
            s.append(step_halfamp)
            s.append(step_height)
            s.append(pk_width)
            s.append(pk_centre)
            steps.append(s)

    return (steps, x, ey)


def find_sign_change(vec):
    """ Return index just before sign change"""
    for i in xrange(0,vec.size-1):
        if (vec[i]>=0 and vec[i+1]<=0) or (vec[i]<=0 and vec[i+1]>=0):
            return i
    return -1

def find_first_less(vec,value):
    for i in xrange(0,vec.size-1):
        if vec[i] < value:
            return i
    return -1

def find_last_less(vec,value):
    for i in xrange(vec.size-1,0,-1):
        if vec[i] < value:
            return i
    return -1


def find_contact_point(x,y,window=20,std_devs=3.0):
    """ Locates the contact point by fitting a line through the contact regime
    and the fully retracted regime, then finding their intersection. """

    retracted=fit_baseline(x,y)
    contact=fit_contact(x,y)
    cx=-(retracted[1]-contact[1])/(retracted[0]-contact[0])
    cy=np.polyval(retracted,cx)
    return (cx,cy)


def fit_baseline(x,y,window=40):
    """ Fit a line through fully-retracted regime to define a baseline.
    Supercedes fit_line for this usage. """

    # If fully retracted position is at end, flip
    if x[-1]>x[0]:
        y=y[::-1]
        x=x[::-1]
    assert x[0]>x[-1]
    
    # How many points to use
    fitlim=significant_diff(y,window,thresh=5.0)

    # Fit line
    P=np.polyfit(x[0:fitlim],y[0:fitlim],1)
    return P

def fit_contact(x,y,window=40):
    """ Fit a line through contact regime to define a contact-line.
    Supercedes fit_slope for this usage. """
    
    # If fully extended position is at end, flip
    if x[-1]<x[0]:
        y=y[::-1]
        x=x[::-1]
    assert x[0]<x[-1]

    # How many points to use
    fitlim=significant_diff(y,window,increasing=False)
    #fitlim=max_diff(y,window)
    # Fit line
    P=np.polyfit(x[0:fitlim],y[0:fitlim],1)
    return P

def significant_diff(y,window,increasing=True,thresh=2.0):
    """ Return last index of first window whose average deviates by more than
    a given number of standard deviations from the rest. """
    
    # Trim to multiple of window, divide into windows
    z=np.reshape(y[0:window*(y.size/window)],(-1,window))
    
    # Find average of each window, and take their differences
    ad=np.abs(np.diff(np.mean(z,axis=1)))

    # Go along list of differences until one differs by more x times mean
    # of the previous
    minWindows=4
    curr_mean=np.mean(ad[0:minWindows+1]) # start with at least four windows
    for i in xrange(minWindows+1,ad.size):
        #print i, curr_mean,ad[i]
        if increasing:
            if ad[i]>curr_mean*thresh:
                break
        else:
            if ad[i]<curr_mean/thresh:
                break
        curr_mean=np.mean(ad[0:i])
    fitlim=i*window
    #print '#',fitlim
    return fitlim

def max_diff(y,window):
    """ Return last index of first window whose average deviates by more than
    a given number of standard deviations from the rest. """
    
    # Trim to multiple of window, divide into windows
    z=np.reshape(y[0:window*(y.size/window)],(-1,window))
    
    # Find average of each window, and take their differences
    ad=np.diff(np.mean(z,axis=1))
    # maybe take std(ad) and use some fraction of that?

    # Go along list of differences until one differs by more than std_devs std devs
    # of the previous
    i=np.argmax(ad)
    fitlim=(i-1)*window
    #print fitlim
    return fitlim

# Fits a line up to a threshold as a percentage of max derivative
# Direction: 0 for downward, 1 for upward
# Defl_thresh=None for no thresholding, i.e. use full data
def fit_line(x,y,direction=0,defl_thresh=0.2,window=10,model='linear'):
    # Deflection threshold to detect contact
    # Smooth out curve by averaging points
    sy = smooth(x,y,window)
    # Calculate derivative
    dy = derivative(x,sy,x[1]-x[0])
    # Smooth derivative
    sdy = smooth(x,dy,window)
    max_dy=np.max(np.abs(sdy))
    if defl_thresh==None:
        return 0,stats.fit(x,sy,model)
    if direction==0:
        i=sdy.size-2
        while i>=0 and np.abs(sdy[i])/max_dy < defl_thresh:
            i = i-1
        # i to end is free approach region
        if len(x[i:-1]) == 0:
            i = 0 # if can't define free approach region, fit all
        return i,stats.fit(x[i:-1],sy[i:-1],model)
    else:
        i=1
        while i<sdy.size-1 and np.abs(sdy[i])/max_dy > defl_thresh:
            i += 1
        # 0 to i is linear regime
        return i,stats.fit(x[0:i],sy[0:i],model)

# Fits a line up to a change in derivative sign
def fit_slope(x,y,direction=0,window=10):
    # Smooth out curve by averaging points
    sy = smooth(x,y,window)
    # Calculate derivative
    dy = derivative(x,sy,x[1]-x[0])
    # Smooth derivative
    sdy = smooth(x,dy,window)
    # Sign
    if direction==0: # Right to left
        i=sdy.size-2
        sgn = sign(sdy[i])
        while i>=0 and sign(sdy[i])==sgn:
            i = i-1
        # i to end is free approach region
        if len(x[i:-1]) == 0:
            i = 0 # if can't define free approach region, fit all
        return i,stats.fit(x[i:-1],sy[i:-1],model='linear')
    else: # Left to right
        i=1
        sgn = np.sign(sdy[i])
        while i<sdy.size-1 and np.sign(sdy[i])==sgn:
            i += 1
        # 0 to i is linear regime
        #print i
        return i,stats.fit(x[0:i],sy[0:i],model='linear')

# Window smoothing  
def smooth(x,y,w,method='mean'):
    if w<1:
        print 'Smoothing window must be greater than zero.'
        return
    sy = np.empty(y.size) # Smoothed y
    for i in range(0,y.size):
        if i < y.size-w:
            a = max(0,i-w)
            b = a + 2*w
        else:
            b = y.size-1
            a = b - 2*w
        if method == 'median':
            sy[i] = np.median(y[a:b])
        else:
            sy[i] = np.mean(y[a:b])
    return sy

# Savitzky-Golay filtering
def sgfilter(x,y,N=3,w=21):
    
    if N>=2:
        # Calculate 1st derivative also
        return sx,sy,sdy
    else:
        return sx,sy

# Finite differening derivative
def derivative(x,y,h,N=9,method='lanczos'):
    if method=='lanczos':
        if N==7:
            fn = super_lanc7
        elif N==9:
            fn = super_lanc9
        elif N==11:
            fn = super_lanc11
    elif method=='central':
        if N==3:
            fn = cent_diff3
        elif N==5:
            fn = cent_diff5
        elif N==7:
            fn = cent_diff7
        elif N==9:
            fn = cent_diff9
    else:
        print 'Unknown differentiator'
        return

    # forces n to be odd
    w = N//2
    dy = np.array([fn(y[stencil(i,w,y.size)],h) for i in range(0,y.size)])
    return dy

def stencil(i,w,n):
    """ Return the indices of stencil of w values around i."""
    ind = np.clip(range(i-w,i+w+1),0,n-1)
    return ind

# See Pavel Holoborodko's website
# Central difference
def cent_diff3(s,h): # stencil size 3
    return 0.5*(s[2]-s[0])/h

def cent_diff5(s,h): # stencil size 5
    return (s[0]-8*s[1]+8*s[3]-s[4])/(12*h)

def cent_diff7(s,h): # stencil size 7
    return (-s[0]+9*s[1]-45*s[2]+45*s[4]-9*s[5]+s[6])/(60*h)

def cent_diff9(s,h): # stencil size 9
    return (3*s[0]-32*s[1]+168*s[2]-672*s[3]+672*s[5]-168*s[6]+32*s[7]-3*s[8])/(840*h)

# Super Lanczos differentiators (order 4)
def super_lanc7(s,h): # Stencil size 7
    return (58*(s[4]-s[2])+67*(s[5]-s[1])-22*(s[6]-s[0]))/(252*h)

def super_lanc9(s,h): # Stencil size 9
    return (126*(s[5]-s[3])+193*(s[6]-s[2])+142*(s[7]-s[1])-86*(s[8]-s[0]))/(1188*h)

def super_lanc11(s,h): # Stencil size 11
    return (296*(s[6]-s[4])+503*(s[7]-s[3])+532*(s[8]-s[2])+294*(s[9]-s[1])-300*(s[10]-s[0]))/(5148*h)

