"""
Initially generated with SMOP 0.41-beta from inst\motivating_example.m

Translated by hand from the original code in motivating_example.m

Tested and found to generate a qualitatively identical plot to the original Octave function.
"""

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
from control import dlqe
from scipy.linalg import solve_discrete_lyapunov, block_diag
from tqdm.auto import tqdm

from ..numba_als_funcs import als_sdp_mrQ


def motivating_example():
    np.random.seed(100)

    Qest_vectr = []
    Rest_vectr = []

    for k in tqdm(range(200)):
        Aa = np.diag([0.1, 0.2, 0.3])
        Aa[0, 2] = 0.1

        Ga = np.array([[1, 2, 3]]).T
        Ca = np.array([[0.1, 0.2, 0]])

        Q_w = np.atleast_2d(0.5)
        R_v = np.atleast_2d(0.1)

        S = solve_discrete_lyapunov(Aa, Ga * Q_w * Ga.T)

        vec_Qw, eig_Qw = np.linalg.eig(Q_w)[::-1]
        vec_Rv, eig_Rv = np.linalg.eig(R_v)[::-1]
        mult_Qw = vec_Qw * np.sqrt(eig_Qw)
        mult_Rv = vec_Rv * np.sqrt(eig_Rv)

        # initial guesses
        G_hat = np.array([[1, 2, 3]]).T
        Qw_hat = .4 * Q_w
        Rv_hat = 4 * R_v

        pa, na = Ca.shape
        _na, ga = Ga.shape

        n = na
        p = pa
        g = ga

        datapts = 1000

        L, _P, _E = dlqe(Aa, G_hat, Ca, Qw_hat, Rv_hat)
        # L, _P, _E = dlqe(Aa, Ga, Ca, Q_w, R_v)
        # L = np.zeros((n, p))

        P = solve_discrete_lyapunov(Aa - Aa @ L @ Ca,
                                    np.hstack([Ga, -Aa @ L]) @ block_diag(Q_w, R_v) @ np.hstack([Ga, -Aa @ L]).T)

        xhat = np.zeros((na, datapts))
        xhat_ = np.zeros((na, datapts + 1))
        x = np.zeros((na, datapts + 1))
        x[:, 0] = 10 * np.ones(na)  # x0
        xhat_[:na, 0] = x[:, 0]  # assume initial state perfectly known

        y = np.zeros((pa, datapts))
        for i in range(datapts):
            y[:, [i]] = Ca @ x[:, [i]] + mult_Rv @ np.random.randn(pa, 1)
            xhat[:, [i]] = xhat_[:, [i]] + L @ (y[:, [i]] - Ca @ xhat_[:, [i]])
            x[:, [i + 1]] = Aa @ x[:, [i]] + Ga @ (mult_Qw @ np.random.randn(ga, 1))
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

        res = als_sdp_mrQ(data, N, model, estimator, plot=False)
        Qest, Rest = res[:2]

        Qest_vectr.append(Qest)
        Rest_vectr.append(Rest)

        k += 1

    fig, ax = plt.subplots()
    ax.plot(np.stack(Qest_vectr).squeeze(), np.stack(Rest_vectr).squeeze(), 'ko')
    ax.set_xlabel('Qw')
    ax.set_ylabel('Rv')
    ax.plot([0.3, 0.7], [0.1, 0.1], "r")
    ax.plot([0.5, 0.5], [0.05, 0.15], "r")
    ax.legend(['estimate', 'true value'])
    fig.show()


if __name__ == "__main__":
    motivating_example()