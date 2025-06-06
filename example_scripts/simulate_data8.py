"""
Initially generated with SMOP 0.41-beta from inst\simulate_data8.m

Tests run. Output appears nigh on identical to Octave implementation.
"""

import copy
from time import process_time
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
from control import dlqe
from scipy.linalg import solve_discrete_lyapunov, block_diag

from ..numba_als_funcs import als_sdp_mrQ

# ALS Implementation Example

def simulate_data8():
    np.random.seed(100)
    Aa = np.diag([0.1, 0.2, 0.3])
    Aa[0, 2] = 0.1
    Aa[0, 1] = 0.1
    Ga = np.eye(3)
    Ca = np.eye(3)
    Q_w = np.diag([0.5, 0.2, 0.1])
    R_v = np.diag([0.5, 0.2, 0.8])
    S = solve_discrete_lyapunov(Aa, Ga @ Q_w @ Ga.T)
    vec_Qw, eig_Qw = np.linalg.eig(Q_w)[::-1]  # this actually returns a namedtuple, but their order is reversed
    vec_Rv, eig_Rv = np.linalg.eig(R_v)[::-1]  # also their type hints are wrong and my IDE complains, so no name access
    mult_Qw = vec_Qw * np.sqrt(eig_Qw)
    mult_Rv = vec_Rv * np.sqrt(eig_Rv)


    # initial guesses
    G_hat = np.eye(3)
    Qw_hat = np.diag([1, 2, 3])
    Rv_hat = 0.001 * R_v
    pa, na = Ca.shape
    _na, ga = Ga.shape
    n = na
    p = pa
    g = ga

    datapts = 5000
    L, _P, _E = dlqe(Aa, G_hat, Ca, Qw_hat, Rv_hat)

    P = solve_discrete_lyapunov(Aa - Aa @ L @ Ca,
                                np.hstack([Ga, -Aa @ L]) @ block_diag(Q_w, R_v) @ np.hstack([Ga, -Aa @ L]).T)

    xhat = np.zeros((na, datapts))
    xhat_ = np.zeros((na, datapts + 1))
    x = np.zeros((na, datapts + 1))
    x[:, 0] = 10 * np.ones(na)
    xhat_[:na, 0] = x[:, 0]

    y = np.zeros((pa, datapts))
    for i in range(datapts):
        y[:, [i]] = Ca @ x[:, [i]] + np.sqrt(mult_Rv) @ np.random.randn(pa, 1)
        xhat[:, [i]] = xhat_[:, [i]] + L @ (y[:, [i]] - Ca @ xhat_[:, [i]])
        x[:, [i + 1]] = Aa @ x[:, [i]] + Ga @ (np.sqrt(mult_Qw) @ np.random.randn(ga, 1))
        xhat_[:, [i + 1]] = Aa @ xhat[:, [i]]

    #########################
    # SETUP ALS PROBLEM #####
    #########################

    model = SimpleNamespace()
    model.A = Aa
    model.C = Ca
    model.G = G_hat
    model.xhat0 = xhat_[:, 0]

    data = SimpleNamespace()
    data.datapts = datapts
    data.yk = y
    data.xhatk = xhat_[:, :-1]
    data.start = 100

    N = 15
    estimator = SimpleNamespace()
    estimator.L = L
    # estimator.Q = Qw_hat;
    # estimator.R = Rv_hat;

    # Using updated codes:
    # This is in general how to call the ALS method:
    res = als_sdp_mrQ(data, N, model, estimator)
    Qest_cell, Rest_cell = res[:2]
    Qest1 = Qest_cell[0]
    Rest1 = Rest_cell[0]

    # Without semidefinite constraints
    # you usually should keep the semidefinite constraints, but could remove them to see how good your model is
    res = als_sdp_mrQ(data, N, model, estimator, use_sdp=False, plot=False)
    Qest_celln, Rest_celln = res[:2]
    Qest_indef = Qest_celln[0]
    Rest_indef = Rest_celln[0]

    # With identity weighting
    # Only use identity weighting if you have a good reason;
    # the default of data-based weighting gives lower variance estimates
    res = als_sdp_mrQ(data, N, model, estimator, weight="I", plot=False)
    Qest_celli, Rest_celli = res[:2]
    QestI = Qest_celli[0]
    RestI = Rest_celli[0]

    # Tradeoff curve example
    # This example shows what to do if you have fewer outputs than states
    # here we pretend that we only know the first row of C
    # we generate a tradeoff curve and look for the elbow in the curve
    data2 = copy.copy(data)
    data2.yk = y[[0], :]

    model2 = copy.copy(model)
    model2.C = model.C[[0], :]

    rho_vec = np.logspace(- 6, 6, 25)

    estimator2 = copy.copy(estimator)
    estimator2.L, _P, _E = dlqe(Aa, G_hat, Ca[0, :], Qw_hat, Rv_hat[0, 0])

    t = process_time()
    res = als_sdp_mrQ(data2, N, model2, estimator2, rho_values=rho_vec)
    Qest_cell2, Rest_cell2, trQ, Phi2 = res[:4]

    print('Total cpu time: %f seconds\n' % (process_time() - t))

    fig, ax = plt.subplots()
    ax.plot(Phi2, trQ, '-o', Phi2[9], trQ[9], 'r*')

    Qest2 = Qest_cell2[9]
    Rest2 = Rest_cell2[9]


if __name__ == "__main__":
    simulate_data8()