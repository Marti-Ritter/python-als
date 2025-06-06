"""
Initially generated with SMOP 0.41-beta from custom script autotime_measurement_als.m
"""

import datetime
import os
import pathlib
from itertools import product

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .autotime_utils import prepare_als_keypoint_time_measurement, measure_als_keypoint_time

def autotime_measurement_als():
    max_n_keypoints = 9
    max_n_dimension = 3

    datapts_steps = [200, 500, 5000]
    repetitions = 5
    N = 15

    rho_vec = np.logspace(-6, 8, 25)
    t_df = pd.DataFrame(columns=['len_x', 'len_z', 'size_Q', 'size_R', 'n_z', 'n_keypoints', 'n_dimensions',
                                'repetition', 'seconds_taken'])

    out_dir = pathlib.Path(os.path.realpath(__file__)).parent / "time_measurement"
    out_dir.mkdir(exist_ok=True, parents=True)

    n_kps = range(1, max_n_keypoints)
    n_dims = range(1, max_n_dimension)
    full_iterator = list(product(n_kps, n_dims, datapts_steps))

    for n_kp, n_dim, n_data in tqdm(full_iterator):
        msg = "Running: {}, {}, {}".format(n_kp, n_dim, n_data)
        print(msg)

        data, model, estimator = prepare_als_keypoint_time_measurement(n_kp, n_dim, n_data)
        len_x = n_kp * n_dim * 2
        len_z = n_kp * n_dim
        size_Q = len_x * len_x
        size_R = len_z * len_z

        for repetition in range(repetitions):
            msg = "Repetition: {}".format(repetition)
            print(msg)
            t = measure_als_keypoint_time(data, N, model, estimator, rho_vec)

            t_df.loc[len(t_df)] = {'len_x': len_x, 'len_z': len_z, 'size_Q': size_Q, 'size_R': size_R, 'n_z': n_data,
                                'n_keypoints': n_kp, 'n_dimensions': n_dim, 'repetition': repetition,
                                'seconds_taken': t}

        out_path = out_dir / ('TimeMeasurement_' + datetime.datetime.now().strftime('%Y%m%d%H%M%S') + '.csv')
        msg = "Writing to {}...".format(out_path)
        print(msg)
        t_df.to_csv(out_path)


if __name__ == "__main__":
    autotime_measurement_als()