# Generated with SMOP  0.41-beta
from libsmop import *
# inst\autocov_calc.m

    
@function
def autocov_calc(Yi=None,N=None,*args,**kwargs):
    varargin = autocov_calc.varargin
    nargin = autocov_calc.nargin

    ## function [bhat_mtr cov_bound bhat] = autocov_calc(Yi,N)
## Given data Yi, estimates autocovariances up to lag N-1
##
## Function inputs :
## Yi - Nd x p matrix of data, where Nd is number of samples
## and p is number of measured variables
## N - number of lags
##
## Function outputs :
## bhat_mtr - N x p matrix of autocovariances
## bhat_mtr = 
## [ <y_1(k)y_1(k)>   <y_1(k)y_2(k)>   ...   <y_p-1(k)y_p(k)>     <y_p(k)y_p(k)>
##	...	              ...                     ...                 ...  
## <y_1(k)y_1(k-N+1)> <y_1(k)y_2(k-N+1)> ... <y_p-1(k)y_p(k-N+1)> <y_p(k)y_p(k-N+1)>];
## cov_bound - estimated standard deviations of autocovariances
## cov_bound = [sigma_11 sigma_12 ... sigma_1p sigpa_pp]'
## bhat - vector of autocovariances as used in ALS
    Eyy=[]
# inst\autocov_calc.m:19
    p=size(Yi,1)
# inst\autocov_calc.m:20
    for i in arange(0,N - 1).reshape(-1):
        temp=dot(Yi(arange(),arange(i + 1,end())),Yi(arange(),arange(1,end() - i)).T)
# inst\autocov_calc.m:22
        temp=temp / (length(Yi) - i)
# inst\autocov_calc.m:23
        Eyy=concat([[Eyy],[temp]])
# inst\autocov_calc.m:24
    
    bhat=ravel(Eyy)
# inst\autocov_calc.m:26
    for i in arange(1,p,1).reshape(-1):
        for j in arange(1,p,1).reshape(-1):
            for k in arange(1,N,1).reshape(-1):
                bhat_mtr[k,dot(p,(i - 1)) + j]=Eyy(dot(p,(k - 1)) + i,j)
# inst\autocov_calc.m:31
    
    npts=size(Yi,2) - N
# inst\autocov_calc.m:36
    for i in arange(1,p,1).reshape(-1):
        for j in arange(1,p,1).reshape(-1):
            cov_bound[i,j]=sqrt(dot(bhat_mtr(1,dot(p,(i - 1)) + i),bhat_mtr(1,dot(p,(j - 1)) + j)) / (npts))
# inst\autocov_calc.m:38
    
    cov_bound=ravel(cov_bound)
# inst\autocov_calc.m:40