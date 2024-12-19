import numpy as np

def UTot_1D(theta,N,alpha,M):                    #Computes total potential energy of system
    theta_i = np.outer(theta,np.ones(N))
    theta_j = np.outer(np.ones(N),theta)
    arg = theta_i - theta_j
    cos_arg = np.cos(arg)
    U = -(alpha*(M/N)**2)*(np.sum(cos_arg)-N)/2
    return U

def OneDHam(theta,dtheta,N,alpha,M):            #Computes 1D hamiltonian 
    U = UTot_1D(theta,N,alpha,M)                 #Total potential energy
    K = (M*np.dot(dtheta,dtheta))/(2*N)         #Total kinetic energy
    H = K+U 
    return H

def OneD_CM(theta):      #Get theta polar coordinate of center of mass for 1D-HMF
    x = np.cos(theta)    #First covert to cartesian coordinates
    y = np.sin(theta)
    xmean = np.mean(x)   #Get the mean values for x and y
    ymean = np.mean(y)  
    theta_CM = np.arctan2(ymean,xmean)       #Find the polar angle corresponding to the point (xmean,ymean)
    theta_CM = np.mod(theta_CM,2.*np.pi)     #Make sure angle is in interval [0,2pi)
    return theta_CM

def Temp1D(dtheta,M,N):                      #Returns temperature of the system
    vsq = dtheta**2
    T = np.mean(vsq)*M/N                     #T = <dtheta**2>*m
    return T 

def TCrit1D(alpha,M,N):       #Critical temperature
    return alpha*(M**2)/(2*N)

def Rho(theta,alpha,M,T,B,N):       #Spatial density for 1D-HMF
    x = (M*B/N*T)
    r0 = M/(2.*np.pi*i0(x))
    Rho = r0*np.exp(-x*np.cos(theta))
    return Rho

