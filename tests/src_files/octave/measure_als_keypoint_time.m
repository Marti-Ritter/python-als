function seconds_taken = measure_als_keypoint_time(data,N,model,estimator,rho_vec)
  t = cputime;
  [Qest_cell,Rest_cell,trQ,Phi2] = als_sdp_mrQ(data,N,model,estimator,'rho_values',rho_vec);
  seconds_taken = cputime - t;
  close("all");
endfunction
