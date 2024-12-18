def Rotate_1D(thetaVals,RotAngle,N): #Rotate coordinate system for 1D-HMF by angle RotAngle. Note that theta dot is invariant under rotation here
    cA = np.cos(RotAngle)
    sA = np.sin(RotAngle)
    R = np.array([[cA,-sA],[sA,cA]])
    x = np.cos(thetaVals)
    y = np.sin(thetaVals)
    Cart = np.zeros((2,N))
    Cart[0],Cart[1] = x,y
    CartR = np.zeros((2,N))
    for i in range(N):
        CartR[:,i] = np.matmul(R,Cart[:,i]) #ith column of CartR is now the ith rotated vector in cartesian coordinates
    theta_rot = np.arctan2(CartR[1],CartR[0])
    return theta_rot 
