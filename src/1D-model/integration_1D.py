"""
integration_1D.py

This module implements a 2nd-order Runge-Kutta method to simulate a 1D-HMF system.

Key Notes:
    - The system assumes periodic boundary conditions where positions (theta) 
      are always wrapped to the interval [0, 2*pi).
    - The state vector y consists of:
        - theta (positions)
        - dtheta (velocities)

Functions:
    - derivs: Computes the time derivatives of theta and dtheta.
    - OneD_Integrate: Integrates the equations of motion using a second-order Runge-Kutta scheme.
"""

import numpy as np

def derivs(y, N, alpha):
    """
    Compute the time derivative of the state vector y.

    Args:
        y (numpy.ndarray): State vector with positions (theta) and velocities (dtheta).
        N (int): Number of particles.
        alpha (float): Coupling constant for the system.

    Returns:
        numpy.ndarray: Time derivatives of theta and dtheta.
    """
    dy = np.zeros_like(y)
    dy[0, :] = y[1, :]  # First derivative of theta comes from dtheta
    xi = np.outer(y[0, :], np.ones(N))
    xj = np.outer(np.ones(N), y[0, :])
    arg = xi - xj
    dy[1, :] = -((alpha * M**2) / N) * np.matmul(np.sin(arg), np.ones(N))
    return dy


def OneD_Integrate(N, Nt, h, theta0, dtheta0, alpha, M):
    """
    Perform time integration of the system.

    Args:
        N (int): Number of particles.
        Nt (int): Number of time steps.
        h (float): Time step size.
        theta0 (float): Initial particle positions.
        dtheta0 (float): Initial particle velocities.
        alpha (float): Coupling constant.
        M (float): Total system mass.

    Returns:
        tuple of numpy.ndarray: Positions (theta) and velocities (dtheta) at each time step.
    """
    yVals = np.zeros((Nt, 2, N))
    y = np.zeros((2, N))
    y[0, :] = theta0
    y[1, :] = dtheta0
    yVals[0, :] = y
    time = np.zeros(Nt)

    for j in range(Nt - 1):
        k1 = h * derivs(y, N, alpha)
        ytmp = y + k1 / 2.0
        k2 = h * derivs(ytmp, N, alpha)
        y = y + k2
        y[0, :] = np.mod(y[0, :], 2.0 * np.pi)  # Wrap theta to [0, 2*pi)
        yVals[j + 1, :] = y
        time[j + 1] = time[j] + h

    return yVals[:, 0, :], yVals[:, 1, :]
