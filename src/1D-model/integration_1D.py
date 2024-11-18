def OneDHam(p,dp,N,alpha,M):                #Computes 1D hamiltonian 
    pi = np.outer(p,np.ones(N))
    pj = np.outer(np.ones(N),p)
    arg = pi - pj
    cos_arg = np.cos(arg)
    U = -(alpha*(M/N)**2)*(np.sum(cos_arg)-N)/2     #Total potential energy
    K = (M*np.dot(dp,dp))/(2*N)           #Total kinetic energy
    H = K+U 
    return H

def derivs(y,N,alpha):
    dy = np.zeros_like(y)
    dy[0,:] = y[1,:]           #dtheta
    xi = np.outer(y[0,:],np.ones(N))
    xj = np.outer(np.ones(N),y[0,:])
    arg = xi - xj               
    dy[1,:] = -((alpha*M**2)/N)*np.matmul(np.sin(arg),np.ones(N)) #d^2theta
    return dy

def OneD_Integrate(N,Nt,h,p0,dp0,alpha,M):
    yVals = np.zeros((Nt,2,N))       #To store phase space coordinates for all times
    y = np.zeros((2,N))
    y[0,:] = p0
    y[1,:] = dp0
    yVals[0,:] = y
    time = np.zeros(Nt)
    
    for j in range(Nt-1):
        k1 = h*derivs(y,N,alpha)
        ytmp = y + k1/2.
        k2 = h*derivs(ytmp,N,alpha)
        y = y + k2
        y[0,:] = np.mod(y[0,:],2.*np.pi)
        yVals[j+1,:] = y 
        time[j+1] = time[j]+h
    return yVals[:,0,:],yVals[:,1,:]
