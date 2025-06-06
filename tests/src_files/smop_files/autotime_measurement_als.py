# Generated with SMOP  0.41-beta
from libsmop import *
# als/autotime_measurement_als.m

clear('all')
pkg('load', 'control')
pkg('load', 'als')
pkg('load', 'dataframe')
max_n_keypoints = 9
# als/autotime_measurement_als.m:7
max_n_dimension = 3
# als/autotime_measurement_als.m:8
datapts_steps = concat([200, 500, 5000])
# als/autotime_measurement_als.m:9
repetitions = 5
# als/autotime_measurement_als.m:10
N = 15
# als/autotime_measurement_als.m:12
rho_vec = logspace(- 6, 8, 25)
# als/autotime_measurement_als.m:13
t_df_struct = cellarray(
    ['len_x', 'len_z', 'size_Q', 'size_R', 'n_z', 'n_keypoints', 'n_dimensions', 'repetition', 'seconds_taken'])
# als/autotime_measurement_als.m:15
t_df = dataframe([], 'colnames', t_df_struct)
# als/autotime_measurement_als.m:16
out_dir = './time_measurement/'
# als/autotime_measurement_als.m:17
for n_kp in arange(1, max_n_keypoints).reshape(-1):
    for n_dim in arange(1, max_n_dimension).reshape(-1):
        for n_data in datapts_steps.reshape(-1):
            msg = strcat('Running: ', num2str(n_kp), ', ', num2str(n_dim), ', ', num2str(n_data))
            # als/autotime_measurement_als.m:22
            disp(msg)
            data, model, estimator = prepare_als_keypoint_time_measurement(n_kp, n_dim, n_data, nargout=3)
            # als/autotime_measurement_als.m:25
            len_x = dot(dot(n_kp, n_dim), 2)
            # als/autotime_measurement_als.m:27
            len_z = dot(n_kp, n_dim)
            # als/autotime_measurement_als.m:28
            size_Q = dot(len_x, len_x)
            # als/autotime_measurement_als.m:29
            size_R = dot(len_z, len_z)
            # als/autotime_measurement_als.m:30
            for repetition in arange(1, repetitions).reshape(-1):
                msg = strcat('Repetition:', num2str(repetition))
                # als/autotime_measurement_als.m:33
                disp(msg)
                t = measure_als_keypoint_time(data, N, model, estimator, rho_vec)
                # als/autotime_measurement_als.m:36
                row_entries = cellarray([len_x, len_z, size_Q, size_R, n_data, n_kp, n_dim, repetition, t])
                # als/autotime_measurement_als.m:38
                if rows(t_df) == 0:
                    t_df = dataframe(row_entries, 'colnames', t_df_struct)
                # als/autotime_measurement_als.m:40
                else:
                    t_df[end() + 1, arange()] = row_entries
            # als/autotime_measurement_als.m:42
            out_path = strcat(out_dir, 'TimeMeasurement_', datestr(now(), 'yyyymmddHHMMSS'), '.csv')
            # als/autotime_measurement_als.m:46
            msg = strcat('Writing to ', out_path, ' ...')
            # als/autotime_measurement_als.m:47
            disp(msg)
            t_df_mat = cell2mat(struct(t_df).x_data)
            # als/autotime_measurement_als.m:50
            csvwrite(out_path, t_df_mat)
