def UTotal(theta,N,alpha,M):                    #Computes total potential energy of system
    theta_i = np.outer(theta,np.ones(N))
    theta_j = np.outer(np.ones(N),theta)
    arg = theta_i - theta_j
    cos_arg = np.cos(arg)
    U = -(alpha*(M/N)**2)*(np.sum(cos_arg)-N)/2
    return U

def OneDHam(theta,dtheta,N,alpha,M):            #Computes 1D hamiltonian 
    U = UTotal(theta,N,alpha,M)                 #Total potential energy
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
