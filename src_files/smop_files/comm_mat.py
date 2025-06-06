# Generated with SMOP  0.41-beta
from libsmop import *
# inst\comm_mat.m

    # Creates mn x mn commutation matrix Kmn so that for any m x n matrix A, 
# Kmn*vec(A) = vec(A')
    
@function
def comm_mat(m=None,n=None,*args,**kwargs):
    varargin = comm_mat.varargin
    nargin = comm_mat.nargin

    Kmn=zeros(dot(m,n),dot(m,n))
# inst\comm_mat.m:4
    for i in arange(1,m,1).reshape(-1):
        for j in arange(1,n,1).reshape(-1):
            Kmn[dot((i - 1),n) + j,dot((j - 1),m) + i]=1
# inst\comm_mat.m:7
    