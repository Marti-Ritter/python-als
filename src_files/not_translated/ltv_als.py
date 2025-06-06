"""
Open source code of linear time-varying autocovariance least-squares (LTV-ALS)
technique to estimate covariances for nonlinear and time-varying models

This function is used to build up the stacks for the ALS
structure of eq. (17) of Lima & Rawlings (2010). These stacks are used to
calculate (Q,R) using semidefinite programming (SDP)
Please see file case1_als.m for an example of the LTV-ALS implementation for 
a simple reactor system
"""

# Generated with SMOP  0.41-beta
import numpy as np
from scipy.linalg import block_diag
from statsmodels.tsa.tsatools import duplication_matrix

# inst\ltv_als.m


def ltv_als(yinn, N, bgn, Ain, Ak, Ck, Gk, Lk):
    """
    Inputs:
    yinn: matrix with vectors of innovations
    N: number of lags used in the autocovariance matrix
    bgn: ALS starting time (k) 
    Ain = Ak-Ak*Lk*Ck
    Ak, Ck, Gk and Lk (EKF gain): time-varying system matrices

    Outputs:
    Qdet, Rdet: covariance estimates
    LHS: left hand side ALS matrix
    Eyy: data vector
    """

    # Load innovations data from k up to k+N-1
    yinn = yinn[:, bgn - 1:]

    p, n = Ck[0].shape
    g = Gk[0].shape[1]

    Eyyfl = yinn.ravel(order="F") @ yinn[:, 0].T  # should be same as vec in octave if it is done with Fortran order
    Eyy = Eyyfl.ravel(order="F")

    Apr = [np.eye(n)]
    for j in range(1, N):
        Apr.append(Ain[bgn + j - 1] @ Apr[j])

    Gamma = []
    for i in range(bgn, 1, -1):
        tempG = []
        for j in range(N):
            tempG.append(Ck[bgn + j - 1] @ Apr[j])
            Apr[j] = Apr[j] @ Ain[i - 1]
        tempG = np.vstack(tempG)
        Gamma.append(tempG)
    Gamma = np.hstack(Gamma[::-1])

    # Calculate Omega1 and Omega2
    AL = np.zeros(bgn)
    for j in range(bgn):
        AL[j] = - Ak[j] * Lk[j]

    Omega1 = block_diag(*Gk[:bgn])
    Omega2 = block_diag(*AL[:bgn])

    PSI = [np.eye(p)]
    Apr1 = np.eye(n)
    for i in range(N):
        PSI.append(-Ck[bgn + i] @ Apr1 @ Ak[bgn] @ Lk[bgn])
        Apr1 = Ain[bgn + i] @ Apr1
    PSI = np.vstack(PSI)

    # Gamma x Omega1
    Gam1 = Gamma @ Omega1
    Gam11 = Gamma[:p, :] @ Omega1

    Gam2 = Gamma @ Omega2
    Gam22 = Gamma[:p, :] @ Omega2

    LHS_Q = 0
    LHS_R = 0
    for i in range(bgn):
        ee = np.eye(bgn - 1)[:, i]
        LHS_Q += np.kron(Gam11 @ np.kron(ee, np.eye(g)), Gam1 @ np.kron(ee, np.eye(g)))
        LHS_R += np.kron(Gam22 @ np.kron(ee, np.eye(p)), Gam2 @ np.kron(ee, np.eye(p)))

    LHS_R += np.kron(np.eye(p), PSI)
    LHS = np.hstack([LHS_Q * duplication_matrix(g), LHS_R * duplication_matrix(p)])

    X = np.linalg.lstsq(LHS, Eyy, rcond=None)[0]
    Qdet = np.reshape(duplication_matrix(g) @ X[:int(g * (g + 1) / 2)], (g, g), order="F")
    Rdet = np.reshape(duplication_matrix(p) @ X[int(g * (g + 1) / 2):], (p, p), order="F")

    return Qdet, Rdet, LHS, Eyy
