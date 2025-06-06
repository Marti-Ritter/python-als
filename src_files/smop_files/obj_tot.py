# Generated with SMOP  0.41-beta
from libsmop import *
# inst\obj_tot.m

    ## Evaluate objection function value, including penalty on tr(Q) and log-barrier term
    
@function
def obj_tot(Q=None,A=None,b=None,W=None,nu=None,vlam_mtr=None,*args,**kwargs):
    varargin = obj_tot.varargin
    nargin = obj_tot.nargin

    r,p=chol(Q,nargout=2)
# inst\obj_tot.m:3
    # p = 0 if Q>0
    if p == 0:
        y=dot(dot((dot(A,ravel(Q)) - b).T,W),(dot(A,ravel(Q)) - b)) + trace(dot(Q,vlam_mtr)) - dot(dot(nu,2),sum(log(diag(r))))
# inst\obj_tot.m:6
    else:
        # if Q is infeasible, objective function is infinity
        y=1e+100
# inst\obj_tot.m:10
    