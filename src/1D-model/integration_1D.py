def derivs(y,N,alpha):                     #Returns the time derivative of the vector y = (theta, dtheta) 
    dy = np.zeros_like(y)
    dy[0,:] = y[1,:]                       #First component of dy/dt, dtheta, comes from second component of y vector
    xi = np.outer(y[0,:],np.ones(N))
    xj = np.outer(np.ones(N),y[0,:])
    arg = xi - xj               
    dy[1,:] = -((alpha*M**2)/N)*np.matmul(np.sin(arg),np.ones(N)) #This is the second derivative of theta 
    return dy

def OneD_Integrate(N,Nt,h,theta0,dtheta0,alpha,M):
    yVals = np.zeros((Nt,2,N))       #To store phase space coordinates for all times
    y = np.zeros((2,N))
    y[0,:] = theta_0
    y[1,:] = dtheta_0
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
