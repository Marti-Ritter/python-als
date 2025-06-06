# Generated with SMOP  0.41-beta
from libsmop import *
# inst\cvode_sens24.m

    
@function
def cvode_sens24(x=None,u=None,t0=None,t1=None,*args,**kwargs):
    varargin = cvode_sens24.varargin
    nargin = cvode_sens24.nargin

    ## Call the cvode integrator and calculate states, x, and 
## sensitivities, sx, at the requested times.
    
    global model
    diff=copy(cputime)
# inst\cvode_sens24.m:8
    
    nx=size(x,1)
# inst\cvode_sens24.m:10
    nu=size(u,1)
# inst\cvode_sens24.m:11
    nt=2
# inst\cvode_sens24.m:12
    
    np=nx + nu
# inst\cvode_sens24.m:14
    sx0=eye(nx,np)
# inst\cvode_sens24.m:15
    param=concat([[x],[u]])
# inst\cvode_sens24.m:17
    dataCV.theta = copy(param + 1e-10)
# inst\cvode_sens24.m:18
    estflag=arange(1,np)
# inst\cvode_sens24.m:20
    
    # CVODES initialization
  # ---------------------
    
    ## Set options for integrator.
    
    options=CVodeSetOptions('UserData',dataCV,'RelTol',model.rtol,'AbsTol',model.atol,'MaxNumSteps',model.odesteps)
# inst\cvode_sens24.m:28
    
    #  if (~isfield(model,'cvodesfun'))
#    disp('a')
    CVodeInit(oderhs1,'BDF','Newton',t0,x,options)
    
    #  else
#    CVodeMalloc (model.cvodesfun, t0, x, options, dataCV);
#  end
    
    ## Set options for forward sensitivity problem.
    
    fsa_options=CVodeSensSetOptions('method','Simultaneous','ErrControl',true,'ParamField','theta','ParamList',estflag)
# inst\cvode_sens24.m:46
    if (logical_not(isfield(model,'dodedu')) or logical_not(isfield(model,'dodedx'))):
        CVodeSensInit(np,[],sx0,fsa_options)
    else:
        CVodeSensInit(np,sensrhs,sx0,fsa_options)
    
    ## Allocate storage for forward sensitivity problem.
#  CVodeSensMalloc (np, 'Simultaneous', sx0, fsa_options);
    
    status,t,x_step,sx_step=CVode(t1,'Normal',nargout=4)
# inst\cvode_sens24.m:60
    if (status == 0):
        x[arange(),2]=x_step
# inst\cvode_sens24.m:63
    else:
        if (status < 0):
            warning('CVode failed with status = %d',status)
            break
    
    nlin=x(arange(),2)
# inst\cvode_sens24.m:69
    Sk=[]
# inst\cvode_sens24.m:70
    Sk=copy(sx_step)
# inst\cvode_sens24.m:71
    A=Sk(arange(),arange(1,nx))
# inst\cvode_sens24.m:72
    B=Sk(arange(),arange(nx + 1,nx + nu))
# inst\cvode_sens24.m:73
    CVodeFree()
    return nlin,A,B
    
if __name__ == '__main__':
    pass
    
    #########################################################
    
@function
def sensrhs(t=None,x=None,xdot=None,xs=None,dataCV=None,*args,**kwargs):
    varargin = sensrhs.varargin
    nargin = sensrhs.nargin

    ## Right hand side of forward sensitivity equation
## sx = dx/dp
## dsx / dt = df/dx * sx + df/dp
    global model
    nx=length(x)
# inst\cvode_sens24.m:84
    p=dataCV.theta
# inst\cvode_sens24.m:85
    u=p(arange(nx + 1,end()))
# inst\cvode_sens24.m:86
    dfdp=concat([model.dodedx(x,t,u),model.dodedu(x,t,u)])
# inst\cvode_sens24.m:88
    
    xsd=dot(model.dodedx(x,t,u),xs) + dfdp
# inst\cvode_sens24.m:89
    flag=0
# inst\cvode_sens24.m:90
    new_data=[]
# inst\cvode_sens24.m:91
    return xsd,flag,new_data
    
if __name__ == '__main__':
    pass
    
    #########################################################
    
    
@function
def oderhs1(t=None,x=None,dataCV=None,*args,**kwargs):
    varargin = oderhs1.varargin
    nargin = oderhs1.nargin

    global model
    u=dataCV.theta(arange(length(x) + 1,end()))
# inst\cvode_sens24.m:97
    xdot=model.odefun(x,t,u)
# inst\cvode_sens24.m:98
    flag=0
# inst\cvode_sens24.m:99
    new_data=[]
# inst\cvode_sens24.m:100
    return xdot,flag,new_data
    
if __name__ == '__main__':
    pass
    