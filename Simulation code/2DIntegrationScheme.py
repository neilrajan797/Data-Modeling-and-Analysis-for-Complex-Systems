scnd_ord = False                 #Can choose if the scheme is of second order or not
thresh = np.sin(np.pi/4)         #Threshold value in first cover for which we switch from one cover to the other in the scheme

def integration(N,Nt,h,alpha,M,phi_0,theta_0,dphi_0,dtheta_0): 
    xVals_1 = np.ones((Nt,2,N))                         #An array to store position values at full steps in first cover. The ith array will have values of phi and theta at the time t = i*dt
    uHalfVals_1 = np.ones((Nt,2,N))                     #An array to store velocity values at midpoint steps in first cover. The ith array will have values of dot(phi) and dot(theta) at the time t = i*dt+dt/2
    uVals_1 = np.ones((Nt,2,N))
    xVals_2 = np.ones((Nt,2,N))                         #for storing velocity values at full steps in first cover. The ith array will have values of phi_dot and theta_dot at the time t = i*dt    xVals_2 = np.ones((Nt,2,N))                         #An array to store position values at full steps in second cover. The ith array will have values of phi and theta at the time t = i*dt
    uHalfVals_2 = np.ones((Nt,2,N))                     #An array to store velocity values at midpoint steps in second cover. The ith array will have values of dot(phi) and dot(theta) at the time t = i*dt+dt/2
    uVals_2 = np.ones((Nt,2,N))                         #for storing velocity values at full steps in second cover. The ith array will have values of psi_dot and chi_dot at the time t = i*dt
    H = np.ones(Nt)
    K = np.ones(Nt)
    U = np.ones(Nt)
    time = np.zeros(Nt)
    
    xVals_1[0,0,:],xVals_1[0,1,:] = phi_0,theta_0                #phi and theta at t = 0
    uVals_1[0,0,:],uVals_1[0,1,:] = dphi_0,dtheta_0              #Their derivatives at t = 0
    
    alpha = alpha*M/N
    K[0],U[0] = K_tot(phi_0,theta_0,dphi_0,dtheta_0,M), U_tot(phi_0,theta_0,alpha,M)
    
    #Get cover 2 coordinates at t = 0 
    xVals_2[0,0,:],xVals_2[0,1,:],uVals_2[0,0,:],uVals_2[0,1,:] = Cov1_to_Cov2(xVals_1[0,0,:],xVals_1[0,1,:],uVals_1[0,0,:],uVals_1[0,1,:])
    #Get interaction sums at time zero
    phi_SumVec_0,theta_SumVec_0,psi_SumVec_0,chi_SumVec_0 = d2p_Sum(phi_0,theta_0,N),d2q_Sum(phi_0,theta_0,N),d2p_Sum(xVals_2[0,0,:],xVals_2[0,1,:],N),d2q_Sum(xVals_2[0,0,:],xVals_2[0,1,:],N)  
    #Compute velocities at t = h/2 
    uHalfVals_1[0,0,:],uHalfVals_1[0,1,:] = kick(h/2,xVals_1[0,0,:],xVals_1[0,1,:],uVals_1[0,0,:],uVals_1[0,1,:],phi_SumVec_0,theta_SumVec_0,scnd_ord,alpha)
    uHalfVals_2[0,0,:],uHalfVals_2[0,1,:] = kick(h/2,xVals_2[0,0,:],xVals_2[0,1,:],uVals_2[0,0,:],uVals_2[0,1,:],psi_SumVec_0,chi_SumVec_0,scnd_ord,alpha)

    for t in range(Nt-1):
        phi,theta,dphi_h,dtheta_h,dphi,dtheta = xVals_1[t,0,:],xVals_1[t,1,:],uHalfVals_1[t,0,:],uHalfVals_1[t,1,:],uVals_1[t,0,:],uVals_1[t,1,:]     #Store cover 1 coordinates at initial time step and velocities at half step ahead
        psi,chi,dpsi_h,dchi_h,dpsi,dchi = xVals_2[t,0,:],xVals_2[t,1,:],uHalfVals_2[t,0,:],uHalfVals_2[t,1,:],uVals_2[t,0,:],uVals_2[t,1,:]           #Store cover 2 coordinates at initial time step and velocities at half step ahead
        phi_SumVec,theta_SumVec,psi_SumVec,chi_SumVec = d2p_Sum(phi,theta,N),d2q_Sum(phi,theta,N),d2p_Sum(psi,chi,N),d2q_Sum(psi,chi,N)     #Vectors containing the interaction term sums for all N particles at initial time step
        phi_temp,theta_temp,psi_temp,chi_temp = np.zeros(N),np.ones(N),np.zeros(N),np.ones(N) #Arrays to store computed position coordinates at h(t+1) that have not been corrected
        for n in range(N):  #Get all N particles' non-corrected position coordinates at h(t+1) to be able to compute interaction terms at h(t+1)
            if np.sin(theta[n])<thresh:
                psi_temp[n],chi_temp[n] = drift(h,psi[n],chi[n],dpsi_h[n],dchi_h[n]) #psi,chi at h(t+1)
                phi_temp[n],theta_temp[n] = TwoToOne_x(psi_temp[n],chi_temp[n])      #phi,theta at h(t+1)
            else:
                phi_temp[n],theta_temp[n] = drift(h,phi[n],theta[n],dphi_h[n],dtheta_h[n]) #phi,theta at h(t+1)
                psi_temp[n],chi_temp[n] = OneToTwo_x(phi_temp[n],theta_temp[n])          #psi,chi at h(t+1)
        phi_SumVec_f,theta_SumVec_f,psi_SumVec_f,chi_SumVec_f = d2p_Sum(phi_temp,theta_temp,N),d2q_Sum(phi_temp,theta_temp,N),d2p_Sum(psi_temp,chi_temp,N),d2q_Sum(psi_temp,chi_temp,N)     #Vectors containing the interaction term sums for all N particles at final time step   
        for n in range(N):
            if np.sin(theta_temp[n])<thresh:
                dpsi_f,dchi_f = kick(h/2,psi_temp[n],chi_temp[n],dpsi_h[n],dchi_h[n],psi_SumVec_f[n],chi_SumVec_f[n],scnd_ord,alpha)   #Get psi,chi derivatives at h(t+1) 
                dphi_f,dtheta_f = TwoToOne_u(psi_temp[n],chi_temp[n],dpsi_f,dchi_f)                                              #Convert Cover2 derivatives at time h(t+1) to cover1 derivatives
                phi_f,theta_f,dphi_f,dtheta_f = sphere(phi_temp[n],theta_temp[n],dphi_f,dtheta_f) 
                psi_f,chi_f,dpsi_f,dchi_f = sphere(psi_temp[n],chi_temp[n],dpsi_f,dchi_f)
                xVals_1[t+1,0,n],xVals_1[t+1,1,n],uVals_1[t+1,0,n],uVals_1[t+1,1,n] = phi_f,theta_f,dphi_f,dtheta_f  #Store cover 1 coordinates at h(t+1)
                xVals_2[t+1,0,n],xVals_2[t+1,1,n],uVals_2[t+1,0,n],uVals_2[t+1,1,n] = psi_f,chi_f,dpsi_f,dchi_f  #Store cover 2 coordinates at h(t+1)
            else:
                dphi_f,dtheta_f = kick(h/2,phi_temp[n],theta_temp[n],dphi_h[n],dtheta_h[n],phi_SumVec_f[n],theta_SumVec_f[n],scnd_ord,alpha)   #Get phi,theta derivatives at h(t+1) 
                dpsi_f,dchi_f = OneToTwo_u(phi_temp[n],theta_temp[n],dphi_f,dtheta_f)   #Get psi,chi derivatives at h(t+1) 
                phi_f,theta_f,dphi_f,dtheta_f = sphere(phi_temp[n],theta_temp[n],dphi_f,dtheta_f)                                        #Correct coordinates
                psi_f,chi_f,dpsi_f,dchi_f = sphere(psi_temp[n],chi_temp[n],dpsi_f,dchi_f)
                xVals_1[t+1,0,n],xVals_1[t+1,1,n],uVals_1[t+1,0,n],uVals_1[t+1,1,n] = phi_f,theta_f,dphi_f,dtheta_f  #Store cover 1 coordinates at h(t+1)
                xVals_2[t+1,0,n],xVals_2[t+1,1,n],uVals_2[t+1,0,n],uVals_2[t+1,1,n] = psi_f,chi_f,dpsi_f,dchi_f  #Store cover 2 coordinates at h(t+1)
        if t < Nt-2:
            uHalfVals_2[t+1,0,:],uHalfVals_2[t+1,1,:] = kick(h/2,xVals_2[t+1,0,:],xVals_2[t+1,1,:],uVals_2[t+1,0,:],uVals_2[t+1,1,:],psi_SumVec_f,chi_SumVec_f,scnd_ord,alpha)   #Get psi,chi derivatives at time h(t+1)+h/2
            uHalfVals_1[t+1,0,:],uHalfVals_1[t+1,1,:] = kick(h/2,xVals_1[t+1,0,:],xVals_1[t+1,1,:],uVals_1[t+1,0,:],uVals_1[t+1,1,:],phi_SumVec_f,theta_SumVec_f,scnd_ord,alpha)   #Get phi,theta derivatives at time h(t+1)+h/2 
        K[t+1],U[t+1] = K_tot(xVals_1[t+1,0,:],xVals_1[t+1,1,:],uVals_1[t+1,0,:],uVals_1[t+1,1,:],M), U_tot(xVals_1[t+1,0,:],xVals_1[t+1,1,:],alpha,M)
        time[t+1] = time[t] + h
    
    return xVals_1[:,0,:],xVals_1[:,1,:],uVals_1[:,0,:],uVals_1[:,1,:]
