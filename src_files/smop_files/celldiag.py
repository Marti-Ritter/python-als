# Generated with SMOP  0.41-beta
from libsmop import *
# inst\celldiag.m

    
@function
def celldiag(F=None,n=None,*args,**kwargs):
    varargin = celldiag.varargin
    nargin = celldiag.nargin

    ## Function for diagonalizing matrices
## function A = celldiag(F,n)
## Function inputs:
## F - cell containing submatrices  
## n - number of rows of zeros to add
## if n < 0, zeros are added to upper right
## if n > 0, zeros are added to lower left
## Function outputs:
## A - block diagonal matrix containing elements of F
## (augmented with zeros if n specified)
## A = [0 0; blockdiag(F) 0] if n < 0
## A = [0 blockdiag(F); 0 0] if n > 0
    if nargin == 1:
        n=0
# inst\celldiag.m:14
    
    len_=prod(size(F))
# inst\celldiag.m:15
    A=F[1]
# inst\celldiag.m:16
    x,y=size(A,nargout=2)
# inst\celldiag.m:17
    for i in arange(2,len_).reshape(-1):
        a,b=size(F[i],nargout=2)
# inst\celldiag.m:19
        A=concat([[A,zeros(x,b)],[zeros(a,y),F[i]]])
# inst\celldiag.m:20
        x=x + a
# inst\celldiag.m:21
        y=y + b
# inst\celldiag.m:22
    
    up=(n < 0)
# inst\celldiag.m:24
    dn=(n > 0)
# inst\celldiag.m:25
    n=abs(n)
# inst\celldiag.m:26
    A=concat([[zeros(dot(up,n),y + n)],[zeros(x,dot(dn,n)),A,zeros(x,dot(up,n))],[zeros(dot(dn,n),y + n)]])
# inst\celldiag.m:27