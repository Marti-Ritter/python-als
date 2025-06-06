# Generated with SMOP  0.41-beta
from libsmop import *
# inst\autocov_calc_thry.m

    
@function
def autocov_calc_thry(A=None,C=None,L=None,Q=None,R=None,N=None,*args,**kwargs):
    varargin = autocov_calc_thry.varargin
    nargin = autocov_calc_thry.nargin

    ## function [b_mtr b] = autocov_calc_thry(A,C,L,Q,R,N)
## Calculates theoretical autocovariances given system matrices
## and noise covariances
##
## Function inputs :
## System matrices A and C
## Estimator gain L
## Process noise covariance Q (note G is assumed identity)
## Measurement noise covariance R
## Number of lags N
##
## Function outputs :
## b_mtr - N x p matrix of autocovariances
## b_mtr = 
## [ E(y_1(k)y_1(k))   E(y_1(k)y_2(k))   ...   E(y_p-1(k)y_p(k))     E(y_p(k)y_p(k))
##	...	              ...                     ...                 ...  
## E(y_1(k)y_1(k-N+1)) E(y_1(k)y_2(k-N+1)) ... E(y_p-1(k)y_p(k-N+1)) E(y_p(k)y_p(k-N+1))];
## cov_bound - estimated standard deviations of autocovariances
## b - vector of autocovariances as used in ALS
    
    Abar=A - dot(dot(A,L),C)
# inst\autocov_calc_thry.m:22
    P=dlyap(Abar,Q + dot(dot(dot(dot(A,L),R),L.T),A.T))
# inst\autocov_calc_thry.m:23
    Eyy=dot(dot(C,P),C.T) + R
# inst\autocov_calc_thry.m:24
    for i in arange(1,N - 1,1).reshape(-1):
        Eyy=concat([[Eyy],[dot(dot(dot(C,Abar ** i),P),C.T) - dot(dot(dot(dot(C,Abar ** (i - 1)),A),L),R)]])
# inst\autocov_calc_thry.m:26
    
    b=ravel(Eyy)
# inst\autocov_calc_thry.m:28
    p=size(C,1)
# inst\autocov_calc_thry.m:30
    for i in arange(1,p,1).reshape(-1):
        for j in arange(1,p,1).reshape(-1):
            for k in arange(1,N,1).reshape(-1):
                b_mtr[k,dot(p,(i - 1)) + j]=Eyy(dot(p,(k - 1)) + i,j)
# inst\autocov_calc_thry.m:34
    