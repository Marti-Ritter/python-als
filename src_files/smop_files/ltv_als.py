# Generated with SMOP  0.41-beta
from libsmop import *
# inst\ltv_als.m

    # Open source code of linear time-varying autocovariance least-squares (LTV-ALS)
# technique to estimate covariances for nonlinear and time-varying models
    
    # This function is used to build up the stacks for the ALS
# structure of eq. (17) of Lima & Rawlings (2010). These stacks are used to
# calculate (Q,R) using semidefinite programming (SDP)
# Please see file case1_als.m for an example of the LTV-ALS implementation for 
# a simple reactor system
    
    
@function
def ltv_als(yinn=None,N=None,bgn=None,Ain=None,Ak=None,Ck=None,Gk=None,Lk=None,*args,**kwargs):
    varargin = ltv_als.varargin
    nargin = ltv_als.nargin

    # Inputs:
    # yinn: matrix with vectors of innovations
    # N: number of lags used in the autocovariance matrix
    # bgn: ALS starting time (k) 
    # Ain = Ak-Ak*Lk*Ck
    # Ak, Ck, Gk and Lk (EKF gain): time-varying system matrices
    
    # Outputs:
    # Qdet, Rdet: covariance estimates
    # LHS: left hand side ALS matrix
    # Eyy: data vector
    
    # Load innovations data from k up to k+N-1
    yinn=yinn(arange(),arange(bgn,end()))
# inst\ltv_als.m:25
    
    p,n=size(Ck[1],nargout=2)
# inst\ltv_als.m:28
    g=columns(Gk[1])
# inst\ltv_als.m:29
    
    Eyyfl=dot(vec(yinn),yinn(arange(),1).T)
# inst\ltv_als.m:32
    Eyy=vec(Eyyfl)
# inst\ltv_als.m:33
    
    Apr[1]=eye(n)
# inst\ltv_als.m:36
    for j in arange(2,N).reshape(-1):
        Apr[j]=dot(Ain[bgn + j - 2],Apr[j - 1])
# inst\ltv_als.m:37
    
    Gamma=[]
# inst\ltv_als.m:39
    for i in arange(bgn,2,- 1).reshape(-1):
        tempG=[]
# inst\ltv_als.m:41
        for j in arange(1,N).reshape(-1):
            tempG=concat([[tempG],[dot(Ck[bgn + j - 1],Apr[j])]])
# inst\ltv_als.m:43
            Apr[j]=dot(Apr[j],Ain[i - 1])
# inst\ltv_als.m:44
        Gamma=concat([tempG,Gamma])
# inst\ltv_als.m:46
    
    # Calculate Omega1 and Omega2
    for j in arange(1,bgn - 1).reshape(-1):
        AL[j]=dot(- Ak[j],Lk[j])
# inst\ltv_als.m:50
    
    
    Omega1=blkdiag(Gk[arange(1,bgn - 1)])
# inst\ltv_als.m:52
    Omega2=blkdiag(AL[arange(1,bgn - 1)])
# inst\ltv_als.m:53
    
    PSI=eye(p)
# inst\ltv_als.m:56
    Apr1=eye(n)
# inst\ltv_als.m:57
    for i in arange(1,N - 1).reshape(-1):
        PSI=concat([[PSI],[dot(dot(dot(- Ck[bgn + i],Apr1),Ak[bgn]),Lk[bgn])]])
# inst\ltv_als.m:60
        Apr1=dot(Ain[bgn + i],Apr1)
# inst\ltv_als.m:61
    
    # Gamma x Omega1
    Gam1=dot(Gamma,Omega1)
# inst\ltv_als.m:65
    
    Gam11=dot(Gamma(arange(1,p),arange()),Omega1)
# inst\ltv_als.m:67
    
    Gam2=dot(Gamma,Omega2)
# inst\ltv_als.m:69
    
    Gam22=dot(Gamma(arange(1,p),arange()),Omega2)
# inst\ltv_als.m:71
    LHS_Q=0
# inst\ltv_als.m:73
    LHS_R=0
# inst\ltv_als.m:74
    for i in arange(1,bgn - 1).reshape(-1):
        ee=eye(bgn - 1)(arange(),i)
# inst\ltv_als.m:78
        LHS_Q=LHS_Q + kron(dot(Gam11,kron(ee,eye(g))),dot(Gam1,kron(ee,eye(g))))
# inst\ltv_als.m:80
        LHS_R=LHS_R + kron(dot(Gam22,kron(ee,eye(p))),dot(Gam2,kron(ee,eye(p))))
# inst\ltv_als.m:81
    
    LHS_R=LHS_R + kron(eye(p),PSI)
# inst\ltv_als.m:85
    LHS=concat([dot(LHS_Q,duplication_matrix(g)),dot(LHS_R,duplication_matrix(p))])
# inst\ltv_als.m:87
    X=ols(Eyy,LHS)
# inst\ltv_als.m:89
    Qdet=reshape(dot(duplication_matrix(g),X(arange(1,dot(g,(g + 1)) / 2),1)),g,g)
# inst\ltv_als.m:91
    Rdet=reshape(dot(duplication_matrix(p),X(arange(dot(g,(g + 1)) / 2 + 1,end()),1)),p,p)
# inst\ltv_als.m:93
    return Qdet,Rdet,LHS,Eyy
    
if __name__ == '__main__':
    pass
    