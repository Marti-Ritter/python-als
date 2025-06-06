import os
from time import process_time
from types import SimpleNamespace

import numpy as np
from control import dlqe
from matplotlib import pyplot as plt

from . import _local_matrices_path
from ..als_funcs import als_sdp_mrQ
from ..numba_als_funcs import als_sdp_mrQ as numba_als_sdp_mrQ


def prepare_als_keypoint_time_measurement(n_keypoints, n_dimensions, datapts):
    """
    Initially generated with SMOP  0.41-beta from custom script prepare_als_keypoint_time_measurement.m
    """

    base_path = _local_matrices_path
    file_end = '_k{}_d{}.csv'.format(n_keypoints, n_dimensions)
    my_F = np.genfromtxt(os.path.join(base_path, 'F' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]
    my_H = np.genfromtxt(os.path.join(base_path, 'H' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]
    my_Q = np.genfromtxt(os.path.join(base_path, 'Q' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]
    my_R = np.genfromtxt(os.path.join(base_path, 'R' + file_end), delimiter=',', filling_values=np.nan)[1:, 1:]

    # random data
    np.random.seed(100)

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

    return data, model, estimator


def measure_als_keypoint_time(data, N, model, estimator, rho_vec):
    """
    Initially generated with SMOP 0.41-beta from custom script autotime_measurement_als.m
    """

    t = process_time()
    _res = als_sdp_mrQ(data, N, model, estimator, rho_values=rho_vec, plot=False)
    seconds_taken = process_time() - t
    plt.close('all')
    return seconds_taken


def measure_numba_als_keypoint_time(data, N, model, estimator, rho_vec):
    """
    Initially generated with SMOP 0.41-beta from custom script autotime_measurement_als.m
    """

    t = process_time()
    _res = numba_als_sdp_mrQ(data, N, model, estimator, rho_values=rho_vec, plot=False)
    seconds_taken = process_time() - t
    plt.close('all')
    return seconds_taken
