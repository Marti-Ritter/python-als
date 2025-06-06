# Generated with SMOP  0.41-beta
from libsmop import *
# inst\ECM_iter.m

    
    ## Function to calculate the state error covariance matrix for the DLQE
## using iterations
    
    
@function
def ECM_iter(A=None,G=None,C=None,Q=None,R=None,*args,**kwargs):
    varargin = ECM_iter.varargin
    nargin = ECM_iter.nargin

    p,n=size(C,nargout=2)
# inst\ECM_iter.m:7
    Pold=zeros(n)
# inst\ECM_iter.m:8
    P=eye(n)
# inst\ECM_iter.m:9
    tol=1e-12
# inst\ECM_iter.m:10
    iter=0
# inst\ECM_iter.m:11
    while (norm(Pold - P,'fro') > tol):

        iter=iter + 1
# inst\ECM_iter.m:14
        Pold=copy(P)
# inst\ECM_iter.m:15
        P=dot(dot(A,Pold),A.T) + dot(dot(G,Q),G.T) - dot(dot(dot(dot(dot(dot(A,Pold),C.T),inv(dot(dot(C,Pold),C.T) + R)),C),Pold),A.T)
# inst\ECM_iter.m:16
        #    norm(Pold-P,"fro")

    
    L=dot(dot(P,C.T),inv(dot(dot(C,P),C.T) + R))
# inst\ECM_iter.m:20