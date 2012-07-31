# Curves force analysis software
# JWA (c)2009

# Statistical functions for fitting

from numpy import *
import scipy.odr as odr
from scipy.interpolate.interpolate import interp1d

sqrt_2_pi = sqrt(2*pi)

# Fit model to x,y data, with choice of model, and returns
# callable function with fitted parameters
def fit(x,y,model='linear',beta0=None):
    # Create model
    #print 'Fitting:',model
    mydata = odr.Data(x,y)
    if model == 'linear':
        mymodel = odr.Model(linear)
        # Set gradient0 to midpoint between min and max
        # Set intercept0 to val of y at x=0
        if beta0==None:
            beta0 = [(max(y)-min(y))/2,y[argmin(abs(x))]]
    elif model == 'gaussian':
        mymodel = odr.Model(gaussian)
        # Set mean0 to x of max y
        # Set stddev0 to half width of data
        # Set offset0 to min y
        if beta0==None:
            beta0 = [x[argmax(y)],(max(x)-min(x))/2]
    elif model == 'lorentzian':
        mymodel = odr.Model(lorentzian)
        # Set location0 to x of max y
        # Set scale0 to half std dev of y
        if beta0==None:
            beta0 = [x[argmax(y)],std(y)/2]
    elif model == 'wormlikechain':
        mymodel = odr.Model(wormlikechain)
                #fjacb=wormlikechain_dB,
                #fjacd=wormlikechain_dx)
        # Set contour length to width of data
        if beta0==None:
            beta0 = [1e-12,(max(x)-min(x))*1.25,0.0001]
        #fixb = [0,1,1]
    elif model == 'offset':
        mymodel = odr.Model(offset)
        if beta0==None:
            beta0 = [median(y)]
    else:
        return None
    myodr = odr.ODR(mydata,mymodel,beta0) 
    myodr.set_job(fit_type=2) # Ordinary Least Squares
    myoutput = myodr.run()
    #myoutput.pprint()
    return FitModel(myoutput.beta,model)
    
# Callable function wrapper for the model, with attached parameters
class FitModel: 
    def __init__(self,beta=None,model='zero'):
        self.beta = beta
        self._model = model
        
    def __call__(self,x):
        return eval(self._model + '(self.beta,x)')
    
# Offset model
# B[0] - constant offset
def offset(B,x):
    return repeat(B[0],len(x))

# Linear model
# B[0] - gradient, B[1] - y intercept
def linear(B,x):
    return B[0]*x + B[1]

# Gaussian model
# B[0] - mean, B[1] - std dev., B[2] - y offset, B[3] - scale parameter
# took out 2,3
def gaussian(B,x):
    # took out B[3] and B[2]
    return exp(-(x-B[0])**2/(2*B[1]**2))/(B[1]*sqrt_2_pi)

# Lorentzian model
# B[0] - location parameter, B[1] - scale parameter, B[2] - y offset
def lorentzian(B,x):
    # took out B[2]
    return (B[1]/((x-B[0])**2 + B[1]**2))/pi

def zero(B,x):
    return 0.0

# Worm-like-chain model
# kT (fixed)
# B[0] - Persistence length
# B[1] - Contour length
# B[2] - offset
# x - extension, nm
kT = 1.38e-23*298*1e9 # Scale to nm?

def wormlikechain(B,x):
    ret = (kT/B[0])*(x/B[1] + .25/(1-x/B[1])**2 - .25) + B[2]
    return ret

def wormlikechain_dx(B,x):
    return (kT/B[0])*(1/B[1] + .5/((1-x/B[1])*B[1]))

def wormlikechain_dB(B,x):
    return vstack([(1/B[1])*(x/B[1] + .25/(1-x/B[1])**2 - .25),
        -(kT/B[0]**2)*(x/B[1] + .25/(1-x/B[1])**2 - .25),
        (kT/B[0])*(.5/((1-x/B[1])*B[1]**2) - x/B[1]**2)])
