# Generated with SMOP  0.41-beta
from libsmop import *
# inst\minrank_ex.m

## Example of calculation of minimum rank Q using the TRADE-OFF CURVE method

clear('all')
close_('all')
randn('seed', 200)
#####################################################
## PLANT SIMULATION #################################
#####################################################

Aa = concat([[0.53565, 0.23972], [0.75473, 0.40629]])
# inst\minrank_ex.m:13
Ca = concat([1, 1])
# inst\minrank_ex.m:15
m, n = size(Ca, nargout=2)
# inst\minrank_ex.m:17
R_v = dot(0.1, eye(m))
# inst\minrank_ex.m:18
Q_w = dot(0.1, eye(n))
# inst\minrank_ex.m:19
Ga = eye(n)
# inst\minrank_ex.m:20
pa, na = size(Ca, nargout=2)
# inst\minrank_ex.m:22
na, ga = size(Ga, nargout=2)
# inst\minrank_ex.m:23
n = copy(na)
# inst\minrank_ex.m:25
p = copy(pa)
# inst\minrank_ex.m:26
g = copy(ga)
# inst\minrank_ex.m:27
datapts = 5000
# inst\minrank_ex.m:29
L, P = dlqe(Aa, eye(n), Ca, Q_w, R_v, nargout=2)
# inst\minrank_ex.m:31
xhat = zeros(na, datapts)
# inst\minrank_ex.m:33
xhat_ = zeros(na, datapts)
# inst\minrank_ex.m:34
x[arange(), 1] = dot(10, ones(na, 1))
# inst\minrank_ex.m:36

xhat_[arange(1, na), 1] = x(arange(), 1)
# inst\minrank_ex.m:38

for i in arange(1, datapts).reshape(-1):
    y[arange(), i] = dot(Ca, x(arange(), i)) + dot(sqrt(R_v), randn(pa, 1))
    # inst\minrank_ex.m:42
    xhat[arange(), i] = xhat_(arange(), i) + dot(L, (y(arange(), i) - dot(Ca, xhat_(arange(), i))))
    # inst\minrank_ex.m:43
    x[arange(), i + 1] = dot(Aa, x(arange(), i)) + dot(Ga, (dot(sqrt(Q_w), randn(ga, 1))))
    # inst\minrank_ex.m:44
    xhat_[arange(), i + 1] = dot(Aa, xhat(arange(), i))
# inst\minrank_ex.m:45

model.A = copy(Aa)
# inst\minrank_ex.m:49
model.C = copy(Ca)
# inst\minrank_ex.m:50
model.G = copy(Ga)
# inst\minrank_ex.m:51
data.datapts = copy(datapts)
# inst\minrank_ex.m:52
data.yk = copy(y)
# inst\minrank_ex.m:53
data.xhatk = copy(xhat_(arange(), arange(1, end() - 1)))
# inst\minrank_ex.m:54
data.start = copy(100)
# inst\minrank_ex.m:55
N = 15
# inst\minrank_ex.m:57
estimator.Q = copy(Q_w)
# inst\minrank_ex.m:59
estimator.R = copy(R_v)
# inst\minrank_ex.m:60
rho = logspace(- 6, 6, 25).T
# inst\minrank_ex.m:61
Qest, Rest, trace_Q, phi_Q = als_sdp_mrQ(data, N, model, estimator, 'rho_values', rho, nargout=4)
# inst\minrank_ex.m:62
# Build Tradeoff plots using calculated data
for i in arange(1, length(Qest)).reshape(-1):
    rank_Q[i] = rank(Qest[i], 0.0001)
# inst\minrank_ex.m:66

figure(1)
plot(phi_Q, trace_Q)
xlabel('\\phi')
ylabel('tr(Q)')
figure(2)
subplot(3, 1, 1)
semilogx(rho, phi_Q)
ylabel('\\phi')
subplot(3, 1, 2)
semilogx(rho, rank_Q, 'ro')
ylabel('rank(Q)')
subplot(3, 1, 3)
semilogx(rho, trace_Q, 'k')
ylabel('tr(Q)')
