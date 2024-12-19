"""
This module provides functions for transforming between coordinate systems
used in the 2D HMF system.

Functions:
- Cov1_to_Cov2: Transforms Cover 1 coordinates (phi, theta, phi_dot, theta_dot) to Cover 2 coordinates (psi, chi, psi_dot, chi_dot).
- Cov2_to_Cov1: Transforms Cover 2 coordinates to Cover 1 coordinates.
- SphereToCart: Convert spherical coordinates to Cartesian coordinates.
- TwoToOne_x: Convert Cover 2 positions to Cover 1 positions without using velocities.
- OneToTwo_x: Convert Cover 1 positions to Cover 2 positions without using velocities.
- OneToTwo_u: Convert velocity components from Cover 1 to Cover 2.
- TwoToOne_u: Convert velocity components Cover 2 to Cover 1. 

Note:
    These functions are designed to perform transformations between coordinate systems 
    and do not enforce wrapping of positions to specific intervals 
    (e.g., [0, 2*pi) for longitude or [0, pi] for latitude).
    Wrapping, if required, is handled in the main integration function or other utilities.
"""

import numpy as np

def Cov1_to_Cov2(phi,theta,phi_dot,theta_dot):  
  """
    Args:
        phi (numpy.ndarray): Longitude in Cover 1.
        theta (numpy.ndarray): Latitude in Cover 1.
        phi_dot (numpy.ndarray): Velocity in phi direction (Cover 1).
        theta_dot (numpy.ndarray): Velocity in theta direction (Cover 1).

    Returns:
        tuple: psi, chi, psi_dot, chi_dot (Cover 2 coordinates).
   """
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    x = cphi*stheta
    y = sphi*stheta
    z = ctheta
    xdot = cphi*ctheta*theta_dot - sphi*stheta*phi_dot
    ydot = sphi*ctheta*theta_dot + cphi*stheta*phi_dot
    zdot = -stheta*theta_dot
    psi = np.arctan2(z,y)
    chi = np.arccos(x)
    psi_dot = (y*zdot - z*ydot)/(z**2 + y**2)
    chi_dot = -xdot/(np.sin(chi))
    return psi,chi,psi_dot,chi_dot

def Cov2_to_Cov1(psi,chi,psi_dot,chi_dot):
    """
    Args:
        psi (numpy.ndarray): Longitude in Cover 2.
        chi (numpy.ndarray): Latitude in Cover 2.
        psi_dot (numpy.ndarray): Velocity in psi direction (Cover 2).
        chi_dot (numpy.ndarray): Velocity in chi direction (Cover 2).

    Returns:
        tuple: phi, theta, phi_dot, theta_dot (Cover 1 coordinates).
    """
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    cchi = np.cos(chi)
    schi = np.sin(chi)
    y = cpsi*schi
    z = spsi*schi
    x = cchi
    ydot = cpsi*cchi*chi_dot-spsi*schi*psi_dot
    zdot = spsi*cchi*chi_dot+cpsi*schi*psi_dot
    xdot = -schi*chi_dot
    phi = np.arctan2(y,x)
    theta = np.arccos(z)
    phi_dot = (x*ydot - y*xdot)/(y**2 + x**2)
    theta_dot = -zdot/np.sin(theta)
    return phi,theta,phi_dot,theta_dot

def OneToTwo_x(phi,theta): 
  """
    Args:
        phi (numpy.ndarray): Longitude in Cover 1.
        theta (numpy.ndarray): Latitude in Cover 1.

    Returns:
        tuple: psi, chi (Cover 2 positions).
    """
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    x = cphi*stheta
    y = sphi*stheta
    z = ctheta
    psi = np.arctan2(z,y)
    chi = np.arccos(x)
    return psi,chi

def TwoToOne_x(psi,chi):
    """
    Args:
        psi (numpy.ndarray): Longitude in Cover 2.
        chi (numpy.ndarray): Latitude in Cover 2.

    Returns:
        tuple: phi, theta (Cover 1 positions).
    """
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    cchi = np.cos(chi)
    schi = np.sin(chi)
    y = cpsi*schi
    z = spsi*schi
    x = cchi
    phi = np.arctan2(y,x)
    theta = np.arccos(z)
    return phi,theta

def OneToTwo_u(phi,theta,phi_dot,theta_dot):
    """
    Args:
        phi (numpy.ndarray): Longitude in Cover 1.
        theta (numpy.ndarray): Latitude in Cover 1.
        phi_dot (numpy.ndarray): Velocity in phi direction (Cover 1).
        theta_dot (numpy.ndarray): Velocity in theta direction (Cover 1).

    Returns:
        tuple: psi_dot, chi_dot (Velocity components in Cover 2).
    """
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    x = cphi*stheta
    y = sphi*stheta
    z = ctheta
    xdot = cphi*ctheta*theta_dot - sphi*stheta*phi_dot
    ydot = sphi*ctheta*theta_dot + cphi*stheta*phi_dot
    zdot = -stheta*theta_dot
    chi = np.arccos(x)
    psi_dot = (y*zdot - z*ydot)/(z**2 + y**2)
    chi_dot = -xdot/(np.sin(chi))
    return psi_dot,chi_dot

def TwoToOne_u(psi,chi,psi_dot,chi_dot):
  """
    Args:
        psi (numpy.ndarray): Longitude in Cover 2.
        chi (numpy.ndarray): Latitude in Cover 2.
        psi_dot (numpy.ndarray): Velocity in psi direction (Cover 2).
        chi_dot (numpy.ndarray): Velocity in chi direction (Cover 2).

    Returns:
        tuple: phi_dot, theta_dot (Velocity components in Cover 1).
    """
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    cchi = np.cos(chi)
    schi = np.sin(chi)
    y = cpsi*schi
    z = spsi*schi
    x = cchi
    ydot = cpsi*cchi*chi_dot-spsi*schi*psi_dot
    zdot = spsi*cchi*chi_dot+cpsi*schi*psi_dot
    xdot = -schi*chi_dot
    theta = np.arccos(z)
    phi_dot = (x*ydot - y*xdot)/(y**2 + x**2)
    theta_dot = -zdot/np.sin(theta)
    return phi_dot,theta_dot

def zeroToPi(...): 
    # Code here...
