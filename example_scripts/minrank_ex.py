"""
Initially generated with SMOP 0.41-beta from inst\minrank_ex.m

Example of calculation of minimum rank Q using the TRADE-OFF CURVE method

Tests run. Results look qualitatively identical to Octave output of the corresponding .m script.
"""

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
from control import dlqe

from ..numba_als_funcs import als_sdp_mrQ


def minrank_ex():
    plt.close('all')
    np.random.seed(200)

    #####################################################
    # PLANT SIMULATION ##################################
    #####################################################

    Aa = np.array([[0.53565, 0.23972], [0.75473, 0.40629]])
    Ca = np.array([[1, 1]])

    m, n = Ca.shape
    R_v = 0.1 * np.eye(m)
    Q_w = 0.1 * np.eye(n)
    Ga = np.eye(n)
    pa, na = Ca.shape
    _na, ga = Ga.shape
    n = na
    p = pa
    g = ga
    datapts = 5000
    L, _P, _E = dlqe(Aa, np.eye(n), Ca, Q_w, R_v)

    xhat = np.zeros((na, datapts))
    xhat_ = np.zeros((na, datapts + 1))
    x = np.zeros((na, datapts + 1))
    x[:, 0] = 10 * np.ones(na)
    xhat_[:na, 0] = x[:, 0]

    y = np.zeros((pa, datapts))
    for i in range(datapts):
        y[:, [i]] = Ca @ x[:, [i]] + np.sqrt(R_v) @ np.random.randn(pa, 1)
        xhat[:, [i]] = xhat_[:, [i]] + L @ (y[:, [i]] - Ca @ xhat_[:, [i]])
        x[:, [i + 1]] = Aa @ x[:, [i]] + Ga @ (np.sqrt(Q_w) @ np.random.randn(ga, 1))
        xhat_[:, [i + 1]] = Aa @ xhat[:, [i]]

    model = SimpleNamespace()
    model.A = Aa
    model.C = Ca
    model.G = Ga

    data = SimpleNamespace()
    data.datapts = datapts
    data.yk = y
    data.xhatk = xhat_[:, :-1]
    data.start = 100

    N = 15
    estimator = SimpleNamespace()
    estimator.Q = Q_w
    estimator.R = R_v
    rho = np.logspace(-6, 6, 25).T
    res = als_sdp_mrQ(data, N, model, estimator, rho_values=rho)
    Qest, Rest, trace_Q, phi_Q = res[:4]

    # Build Tradeoff plots using calculated data
    rank_Q = []
    for i in range(len(Qest)):
        rank_Q.append(np.linalg.matrix_rank(Qest[i], 0.0001))

    fig, ax = plt.subplots()
    ax.plot(phi_Q, trace_Q)
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('tr(Q)')

    fig, axs = plt.subplots(3, 1)
    axs[0].set_xscale('log')
    axs[0].plot(rho, phi_Q)
    axs[0].set_ylabel('$\\phi$')

    axs[1].set_xscale('log')
    axs[1].plot(rho, rank_Q, 'ro')
    axs[1].set_ylabel('rank(Q)')

    axs[2].set_xscale('log')
    axs[2].plot(rho, trace_Q, 'k')
    axs[2].set_ylabel('tr(Q)')
    axs[2].set_xlabel('$\\rho$')

    plt.tight_layout()


if __name__ == "__main__":
    minrank_ex()