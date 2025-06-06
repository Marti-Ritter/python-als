function [data,model,estimator] = prepare_als_keypoint_time_measurement(n_keypoints,n_dimensions,datapts)
  base_path = "./sample_matrices/";
  file_end = strcat("_k", num2str(n_keypoints), "_d", num2str(n_dimensions), ".csv");

  my_F = csvread(strcat(base_path, "F", file_end))(2:end, 2:end);
  my_H = csvread(strcat(base_path, "H", file_end))(2:end, 2:end);
  my_Q = csvread(strcat(base_path, "Q", file_end))(2:end, 2:end);
  my_R = csvread(strcat(base_path, "R", file_end))(2:end, 2:end);

  %%% random data
  randn('seed',100);

  Aa = my_F;
  Ca = my_H;

  Q_w = my_Q;
  R_v = my_R;
  [vec_Qw,eig_Qw] = eig(Q_w);
  [vec_Rv,eig_Rv] = eig(R_v);
  mult_Qw = vec_Qw*sqrt(eig_Qw);
  mult_Rv = vec_Rv*sqrt(eig_Rv);

  %%% initial guesses
  Qw_hat = 1e-2*Q_w;
  Rv_hat= 1e-3*R_v;

  [pa,na]=size(Ca);

  n=na;
  p=pa;

  L=dlqe(Aa, [], Ca, Qw_hat, Rv_hat);
  xhat=zeros(na,datapts);
  xhat_=zeros(na,datapts);

  x(:,1) = 10*ones(na,1);  %% x0
  xhat_(1:na,1) = x(:,1); %% assume initial state perfectly known

  for i = 1:datapts
    y(:,i) = Ca*x(:,i)+mult_Rv*randn(pa,1);
    xhat(:,i) = xhat_(:,i) + L*(y(:,i)-Ca*xhat_(:,i));
    x(:,i+1) = Aa*x(:,i);
    xhat_(:,i+1) = Aa*xhat(:,i);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SETUP ALS PROBLEM %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%

  model.A = my_F;
  model.C = my_H;
  model.xhat0 = xhat_(:,1);

  data.datapts = datapts;
  data.yk = y;
  data.start = 100;

  estimator.Q = my_Q;
  estimator.R = my_R;
endfunction
