# Generated with SMOP  0.41-beta
from libsmop import *
# als\keypoint_simulation.m

## ALS Speed Test Implementation Example

clear('all')
close_('all')
n_keypoints = 2
# als\keypoint_simulation.m:6
n_dimensions = 1
# als\keypoint_simulation.m:7
base_path = './sample_matrices/'
# als\keypoint_simulation.m:9
file_end = strcat('_k', num2str(n_keypoints), '_d', num2str(n_dimensions), '.csv')
# als\keypoint_simulation.m:10
my_F = csvread(strcat(base_path, 'F', file_end))(arange(2, end()), arange(2, end()))
# als\keypoint_simulation.m:12
my_H = csvread(strcat(base_path, 'H', file_end))(arange(2, end()), arange(2, end()))
# als\keypoint_simulation.m:13
my_Q = csvread(strcat(base_path, 'Q', file_end))(arange(2, end()), arange(2, end()))
# als\keypoint_simulation.m:14
my_R = csvread(strcat(base_path, 'R', file_end))(arange(2, end()), arange(2, end()))
# als\keypoint_simulation.m:15
### random data
randn('seed', 100)
datapts = 50000
# als\keypoint_simulation.m:19
Aa = copy(my_F)
# als\keypoint_simulation.m:21
Ca = copy(my_H)
# als\keypoint_simulation.m:22
Q_w = copy(my_Q)
# als\keypoint_simulation.m:24
R_v = copy(my_R)
# als\keypoint_simulation.m:25
vec_Qw, eig_Qw = eig(Q_w, nargout=2)
# als\keypoint_simulation.m:26
vec_Rv, eig_Rv = eig(R_v, nargout=2)
# als\keypoint_simulation.m:27
mult_Qw = dot(vec_Qw, sqrt(eig_Qw))
# als\keypoint_simulation.m:28
mult_Rv = dot(vec_Rv, sqrt(eig_Rv))
# als\keypoint_simulation.m:29
### initial guesses
Qw_hat = copy(Q_w)
# als\keypoint_simulation.m:32
Rv_hat = copy(R_v)
# als\keypoint_simulation.m:33
pa, na = size(Ca, nargout=2)
# als\keypoint_simulation.m:35
n = copy(na)
# als\keypoint_simulation.m:37
p = copy(pa)
# als\keypoint_simulation.m:38
### generatred data
L = dlqe(Aa, [], Ca, Qw_hat, Rv_hat)
# als\keypoint_simulation.m:42
xhat = zeros(na, datapts)
# als\keypoint_simulation.m:43
xhat_ = zeros(na, datapts)
# als\keypoint_simulation.m:44
x[arange(), 1] = dot(10, ones(na, 1))
# als\keypoint_simulation.m:46

xhat_[arange(1, na), 1] = x(arange(), 1)
# als\keypoint_simulation.m:47

for i in arange(1, datapts).reshape(-1):
    y[arange(), i] = dot(Ca, x(arange(), i)) + dot(mult_Rv, randn(pa, 1))
    # als\keypoint_simulation.m:50
    xhat[arange(), i] = xhat_(arange(), i) + dot(L, (y(arange(), i) - dot(Ca, xhat_(arange(), i))))
    # als\keypoint_simulation.m:51
    x[arange(), i + 1] = dot(Aa, x(arange(), i))
    # als\keypoint_simulation.m:52
    xhat_[arange(), i + 1] = dot(Aa, xhat(arange(), i))
# als\keypoint_simulation.m:53

### own data
my_data = csvread('./long_2k_2d_data_y.csv')
# als\keypoint_simulation.m:57
my_data = my_data(arange(4, end()), arange(2, end()))
# als\keypoint_simulation.m:58
my_data = transpose(my_data)
# als\keypoint_simulation.m:59
z_dim, datapts = size(my_data, nargout=2)
# als\keypoint_simulation.m:60
xhat_ = zeros(na, datapts)
# als\keypoint_simulation.m:61
xhat_[arange(1, z_dim), 1] = my_data(arange(), 1)
# als\keypoint_simulation.m:62

y = copy(my_data)
# als\keypoint_simulation.m:63
#########################
### SETUP ALS PROBLEM ###
#########################

model.A = copy(my_F)
# als\keypoint_simulation.m:68
model.C = copy(my_H)
# als\keypoint_simulation.m:69
model.xhat0 = copy(xhat_(arange(), 1))
# als\keypoint_simulation.m:70
data.datapts = copy(datapts)
# als\keypoint_simulation.m:72
data.yk = copy(y)
# als\keypoint_simulation.m:73
data.start = copy(100)
# als\keypoint_simulation.m:74
estimator.Q = copy(my_Q)
# als\keypoint_simulation.m:76
estimator.R = copy(my_R)
# als\keypoint_simulation.m:77
# estimator.Q = cell2mat(Qest_cell(15));
# estimator.R = cell2mat(Rest_cell(15));

N = 15
# als\keypoint_simulation.m:82
rho_vec = logspace(- 6, 8, 25)
# als\keypoint_simulation.m:83
t = copy(cputime)
# als\keypoint_simulation.m:85
Qest_cell, Rest_cell, trQ, Phi2 = als_sdp_mrQ(data, N, model, estimator, 'rho_values', rho_vec, nargout=4)
# als\keypoint_simulation.m:86
printf('Total cpu time: %f seconds\\n', cputime - t)
figure
plot(Phi2, trQ, '-o')
