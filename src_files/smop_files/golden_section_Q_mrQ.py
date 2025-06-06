# Generated with SMOP  0.41-beta
from libsmop import *
# inst\golden_section_Q_mrQ.m

    
@function
def golden_section_Q_mrQ(x00=None,x33=None,Q=None,delQ=None,A=None,b=None,W=None,nu=None,vlam_mtr=None,*args,**kwargs):
    varargin = golden_section_Q_mrQ.varargin
    nargin = golden_section_Q_mrQ.nargin

    ## Uses golden section method to calculate optimal step size for given search direction
## Objective function: y =  norm(A*Q(:) - b)^2_w + trace(Q*vlam_mtr) - nu*log(det(Q));
    
    tol=0.001
# inst\golden_section_Q_mrQ.m:5
    x0=copy(x00)
# inst\golden_section_Q_mrQ.m:6
    x3=copy(x33)
# inst\golden_section_Q_mrQ.m:7
    alpha=(3 - sqrt(5)) / 2
# inst\golden_section_Q_mrQ.m:8
    x1=x0 + dot(alpha,(x3 - x0))
# inst\golden_section_Q_mrQ.m:9
    x2=x3 - dot(alpha,(x3 - x0))
# inst\golden_section_Q_mrQ.m:10
    Q0=Q + dot(x0,delQ)
# inst\golden_section_Q_mrQ.m:12
    J0=obj_tot(Q0,A,b,W,nu,vlam_mtr)
# inst\golden_section_Q_mrQ.m:13
    Q1=Q + dot(x1,delQ)
# inst\golden_section_Q_mrQ.m:14
    J1=obj_tot(Q1,A,b,W,nu,vlam_mtr)
# inst\golden_section_Q_mrQ.m:15
    Q2=Q + dot(x2,delQ)
# inst\golden_section_Q_mrQ.m:16
    J2=obj_tot(Q2,A,b,W,nu,vlam_mtr)
# inst\golden_section_Q_mrQ.m:17
    Q3=Q + dot(x3,delQ)
# inst\golden_section_Q_mrQ.m:18
    J3=obj_tot(Q3,A,b,W,nu,vlam_mtr)
# inst\golden_section_Q_mrQ.m:19
    while (norm(x3 - x0)) > tol:

        if (J2 < J1):
            x0=copy(x1)
# inst\golden_section_Q_mrQ.m:22
            J0=copy(J1)
# inst\golden_section_Q_mrQ.m:23
            x1=copy(x2)
# inst\golden_section_Q_mrQ.m:24
            J1=copy(J2)
# inst\golden_section_Q_mrQ.m:25
            x2=x3 - dot(alpha,(x3 - x0))
# inst\golden_section_Q_mrQ.m:26
            Q2=Q + dot(x2,delQ)
# inst\golden_section_Q_mrQ.m:27
            J2=obj_tot(Q2,A,b,W,nu,vlam_mtr)
# inst\golden_section_Q_mrQ.m:28
        else:
            x3=copy(x2)
# inst\golden_section_Q_mrQ.m:30
            J3=copy(J2)
# inst\golden_section_Q_mrQ.m:31
            x2=copy(x1)
# inst\golden_section_Q_mrQ.m:32
            J2=copy(J1)
# inst\golden_section_Q_mrQ.m:33
            x1=x0 + dot(alpha,(x3 - x0))
# inst\golden_section_Q_mrQ.m:34
            Q1=Q + dot(x1,delQ)
# inst\golden_section_Q_mrQ.m:35
            J1=obj_tot(Q1,A,b,W,nu,vlam_mtr)
# inst\golden_section_Q_mrQ.m:36

    
    if (J1 < J2):
        h=copy(x1)
# inst\golden_section_Q_mrQ.m:41
    else:
        h=copy(x2)
# inst\golden_section_Q_mrQ.m:43
    