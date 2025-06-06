# Generated with SMOP  0.41-beta
from libsmop import *
# inst\obj_ls.m

    ## Evaluates only least-squares portion of objective 
## (for use in choosing Q from tradeoff curve)
    
@function
def obj_ls(Q=None,A=None,b=None,W=None,*args,**kwargs):
    varargin = obj_ls.varargin
    nargin = obj_ls.nargin

    y=dot(dot((dot(A,ravel(Q)) - b).T,W),(dot(A,ravel(Q)) - b))
# inst\obj_ls.m:4
    if imag(y) != 0:
        y=10000000000.0
# inst\obj_ls.m:5
    