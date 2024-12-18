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
