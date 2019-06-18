# convert from the Trace3d twiss parameter to the IMPACT input distribution parameters.
# Here, the TRACE3d units are, alpha, beta (m/rad), emittance (mm-mrad, rms, unnormalized)
#                              alpha_z,beta_z (deg/keV), emittance (deg-keV,rms)
#frq(Hz),mass(eV£©,kine(eV)


import math

def Twiss2Sigma(alpha,beta,emittance,freq,mass,kine):
    Egamma = 1 + kine/mass
    Ebeta = math.sqrt(1.0-1.0/Egamma/Egamma)
    xl=2.99792458e8/(2.0*math.pi*freq)     # Length scale
    
    emittance = emittance *1.0e-6 # unnormalized emittance
    
    sig1=math.sqrt(beta*emittance/(1.0+alpha**2))
    sig2=math.sqrt(emittance/beta)
    rho=alpha/math.sqrt(1.0+alpha**2)
    
    sig1 = sig1 / xl
    sig2 = sig2 * Egamma*Ebeta

    #gammabeta is needed because <x'^2>=emittance*gamma instead of <(px/mc)^2>
    
    return sig1, sig2, rho

def Twiss2SigmaZ(alpha,beta,emittance,freq,mass,kine):
    Egamma = 1 + kine/mass
    Ebeta = math.sqrt(1.0-1.0/Egamma/Egamma)
    xl=2.99792458e8/(2.0*math.pi*freq)     # Length scale
    
    beta=1000.0*beta
    eps=1.0e-3*emittance
    sig1=math.sqrt(beta*eps/(1.0+alpha**2))/(180.0/math.pi)
    sig2=math.sqrt(eps/beta)/mass*1.0e6
    rho=-alpha/math.sqrt(1.0+alpha**2)
    
    return sig1, sig2, rho

def Sigma2Twiss(sig1, sig2, rho, freq,mass,kine):
    Egamma = 1 + kine/mass
    Ebeta = math.sqrt(1.0-1.0/Egamma/Egamma)
    xl=2.99792458e8/(2.0*math.pi*freq)     # Length scale
    
    sig1        = sig1 * xl
    #gammabeta is needed because <x'^2>=emittance*gamma instead of <(px/mc)^2>
    sig2        = sig2 / (Egamma*Ebeta)

    alpha       = math.sqrt(rho**2/(1-rho**2)) * (1 if rho>0 else -1) 
    beta        = sig1/sig2 * math.sqrt(1+alpha**2)
    emittance   = sig1*sig2 * math.sqrt(1+alpha**2)

    emittance   = emittance / 1.0e-6 #mm*mrad
    
    return alpha,beta,emittance

def Sigma2TwissZ(sig1, sig2, rho, freq,mass,kine):
    Egamma = 1 + kine/mass
    Ebeta = math.sqrt(1.0-1.0/Egamma/Egamma)
    xl=2.99792458e8/(2.0*math.pi*freq)     # Length scale
    
    sig1        = sig1*(180.0/math.pi)
    #gammabeta is needed because <x'^2>=emittance*gamma instead of <(px/mc)^2>
    sig2        = sig2 *mass/1.0e6
    
    alpha       = math.sqrt(rho**2/(1-rho**2)) * (1 if rho>0 else -1) 
    beta        = sig1/sig2 * math.sqrt(1+alpha**2)  
    emittance   = sig1*sig2 * math.sqrt(1+alpha**2)
    
    beta        = beta      / 1000.0    #deg/keV
    emittance   = emittance / 1.0e-3    #deg*keV
    
    return alpha,beta,emittance