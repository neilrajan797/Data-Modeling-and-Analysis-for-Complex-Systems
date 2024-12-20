"""
This module contains the kick and drift functions used in the
Kick-Drift-Kick integration scheme for the 2D HMF system.

Functions:
- kick: Updates the velocity components based on the interaction forces and optional second-order corrections. 
- drift: Updates the position components based on the velocities and previous positions.
"""

import numpy as np

def kick(h,p,q,dp,dq,d2pSum,d2qSum,scnd_ord,alpha):
  """
    Args:
        h (float): Time step size.
        p (numpy.ndarray): Current longitudes (phi or psi).
        q (numpy.ndarray): Current latitudes (theta or chi).
        dp (numpy.ndarray): Current velocities in longitude direction.
        dq (numpy.ndarray): Current velocities in latitude direction.
        d2pSum (numpy.ndarray): Sum in interaction terms for longitude acceleration.
        d2qSum (numpy.ndarray): Sum in interaction terms for latitude acceleration.
        scnd_ord (bool): Flag to include second-order corrections.
        alpha (float): Coupling constant for the system.

    Returns:
        tuple: Updated dp and dq (velocities after the kick step).
    """
    if alpha == 0:
        if scnd_ord == True:    
            d2phi_a = -2.*h*dp*dq/np.tan(q)
            d2phi_b = h**2*(2.*dp*dq**2/np.tan(q)**2 - dp**3*np.cos(q)**2)
            d2theta_a = (h/2.)*np.sin(2*q)*dp**2
            d2theta_b = -2*h**2*dp*dq**2*np.cos(q)**2
            return dp + d2phi_a + d2phi_b, dq + d2theta_a + d2theta_b
        else: 
            d2phi_a = -2.*h*dp*dq/np.tan(q)
            d2theta_a = (h/2.)*np.sin(2*q)*dp**2
            return dp + d2phi_a, dq + d2theta_a
    else: 
        if scnd_ord == True:    
            d2phi_a = -2.*h*dp*dq/np.tan(q)
            d2phi_b = h**2*(2.*dp*dq**2/np.tan(q)**2 - dp**3*np.cos(q)**2)
            d2theta_a = (h/2.)*np.sin(2*q)*dp**2
            d2theta_b = -2*h**2*dp*dq**2*np.cos(q)**2
            return dp + d2phi_a + d2phi_b, dq + d2theta_a + d2theta_b
        else: 
            d2phi_a = h*(alpha*(d2pSum/np.sin(q))-(2.*dp*dq/np.tan(q)))
            d2theta_a = h*(alpha*d2qSum + (1/2.)*np.sin(2*q)*dp**2)
            return dp + d2phi_a, dq + d2theta_a


def drift(h,p,q,dp,dq):
  """
    Update position components (p, q) based on the current velocities.

    Args:
        h (float): Time step size.
        p (numpy.ndarray): Current longitudes (phi or psi).
        q (numpy.ndarray): Current latitudes (theta or chi).
        dp (numpy.ndarray): Current velocities in longitude direction.
        dq (numpy.ndarray): Current velocities in latitude direction.

    Returns:
        tuple: Updated p and q (positions after the drift step).
    """
    p += h*dp
    q += h*dq
    return p, q
