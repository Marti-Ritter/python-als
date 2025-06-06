"""
Initially generated with SMOP 0.41-beta from inst\simulate_data8_diag.m

Diagonal ALS Implementation Example

Tests run. See als_diag for more information.
"""

from types import SimpleNamespace

# Generated with SMOP  0.41-beta
import numpy as np
from control import dlqe
from scipy.linalg import solve_discrete_lyapunov, block_diag

from ..numba_als_funcs import als_diag

def simulate_data8_diag():
    np.random.seed(100)
    Aa = np.diag([0.1, 0.2, 0.3])
    Aa[0, 2] = 0.1
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
    L, _P, _E = dlqe(Aa, G_hat, Ca, Qw_hat, Rv_hat, method="scipy")

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

    res = als_diag(data, N, model, estimator)
    Qest, Rest, Lest, As, bhat = res[:5]
    return Qest, Rest, Lest, As, bhat


if __name__ == "__main__":
    simulate_data8_diag()