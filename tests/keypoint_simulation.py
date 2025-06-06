"""
Initially generated with SMOP 0.41-beta from custom script keypoint_simulation.m

ALS Speed Test Implementation Example
"""

import os
import pathlib
from time import process_time
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
from control import dlqe

from . import _local_matrices_path
from ..als_funcs import als_sdp_mrQ

def keypoint_simulation():
    plt.close('all')

    n_keypoints = 3
    n_dimensions = 1

    base_path = _local_matrices_path
    file_end = '_k{}_d{}.csv'.format(n_keypoints, n_dimensions)
    my_F = np.genfromtxt(os.path.join(base_path, 'F' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]
    my_H = np.genfromtxt(os.path.join(base_path, 'H' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]
    my_Q = np.genfromtxt(os.path.join(base_path, 'Q' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]
    my_R = np.genfromtxt(os.path.join(base_path, 'R' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]

    # random data
    np.random.seed(100)
    datapts = 50000

    Aa = my_F
    Ca = my_H
    Q_w = my_Q
    R_v = my_R

    vec_Qw, eig_Qw = np.linalg.eig(Q_w)[::-1]
    vec_Rv, eig_Rv = np.linalg.eig(R_v)[::-1]
    mult_Qw = vec_Qw * np.sqrt(eig_Qw)
    mult_Rv = vec_Rv * np.sqrt(eig_Rv)

    # initial guesses
    Qw_hat = Q_w
    Rv_hat = R_v
    pa, na = Ca.shape
    n = na
    p = pa

    # generated data
    L, _P, _E = dlqe(Aa, np.eye(len(Aa)), Ca, Qw_hat, Rv_hat, method="scipy")

    xhat = np.zeros((na, datapts))
    xhat_ = np.zeros((na, datapts + 1))
    x = np.zeros((na, datapts + 1))
    x[:, 0] = 10 * np.ones(na)
    xhat_[:na, 0] = x[:, 0]

    y = np.zeros((pa, datapts))
    for i in range(datapts):
        y[:, [i]] = Ca @ x[:, [i]] + np.sqrt(mult_Rv) @ np.random.randn(pa, 1)
        xhat[:, [i]] = xhat_[:, [i]] + L @ (y[:, [i]] - Ca @ xhat_[:, [i]])
        x[:, [i + 1]] = Aa @ x[:, [i]]
        xhat_[:, [i + 1]] = Aa @ xhat[:, [i]]

    # own data
    data_path = pathlib.Path(os.path.realpath(__file__)).parent / "data" / "test_trace.csv"
    my_data = np.genfromtxt(data_path, delimiter=',', filling_values=np.nan)[3:, 1:].T
    z_dim, datapts = my_data.shape
    xhat_ = np.zeros((na, datapts))
    xhat_[0:z_dim, 0] = my_data[:, 0]

    y = my_data

    #########################
    # SETUP ALS PROBLEM #####
    #########################

    model = SimpleNamespace()
    model.A = my_F
    model.C = my_H
    model.xhat0 = xhat_[:, 0]

    data = SimpleNamespace()
    data.datapts = datapts
    data.yk = y
    data.start = 100

    estimator = SimpleNamespace()
    estimator.Q = my_Q
    estimator.R = my_R
    # estimator.Q = cell2mat(Qest_cell(15));
    # estimator.R = cell2mat(Rest_cell(15));

    N = 15
    rho_vec = np.logspace(- 6, 8, 25)
    t = process_time()
    res = als_sdp_mrQ(data, N, model, estimator, rho_values=rho_vec)
    Qest_cell, Rest_cell, trQ, Phi2 = res[:4]
    print('Total cpu time: {} seconds'.format(process_time() - t))

    # plot results
    fig, ax = plt.subplots()
    ax.plot(Phi2, trQ, marker='o', zorder=0)
    ax.set_xlabel('Phi2')
    ax.set_ylabel('trace(Q)')
    ax.scatter(Phi2[14], trQ[14], color='r', zorder=1)
    fig.show()


if __name__ == "__main__":
    keypoint_simulation()