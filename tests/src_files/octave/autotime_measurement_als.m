clear all

pkg load control
pkg load als
pkg load dataframe

max_n_keypoints = 9;
max_n_dimension = 3;
datapts_steps = [200, 500, 5000];
repetitions = 5;

N = 15;
rho_vec = logspace(-6,8,25);

t_df_struct = {"len_x", "len_z", "size_Q", "size_R", "n_z", "n_keypoints", "n_dimensions", "repetition", "seconds_taken"};
t_df = dataframe([], "colnames", t_df_struct);
out_dir = "./time_measurement/";

for n_kp = 1:max_n_keypoints
  for n_dim = 1:max_n_dimension
    for n_data = datapts_steps
      msg = strcat("Running: ", num2str(n_kp), ", ", num2str(n_dim), ", ", num2str(n_data));
      disp(msg);

      [data,model,estimator] = prepare_als_keypoint_time_measurement(n_kp,n_dim,n_data);

      len_x = n_kp*n_dim*2;
      len_z = n_kp*n_dim;
      size_Q = len_x * len_x;
      size_R = len_z * len_z;

      for repetition = 1:repetitions
        msg = strcat("Repetition:", num2str(repetition));
        disp(msg);

        t = measure_als_keypoint_time(data,N,model,estimator,rho_vec);

        row_entries = {len_x, len_z, size_Q, size_R, n_data, n_kp, n_dim, repetition, t};
        if rows(t_df) == 0
          t_df = dataframe(row_entries, "colnames", t_df_struct);
        else
          t_df(end+1,:) = row_entries;
        endif
      endfor

      out_path = strcat(out_dir, "TimeMeasurement_", datestr(now(), 'yyyymmddHHMMSS'), ".csv");
      msg = strcat("Writing to ", out_path, " ...");
      disp(msg);

      t_df_mat = cell2mat(struct(t_df).x_data);
      csvwrite(out_path, t_df_mat);
    endfor
  endfor
endfor

