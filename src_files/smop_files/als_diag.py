# Generated with SMOP  0.41-beta
from libsmop import *
# inst\als_diag.m

    
    ## The Diagonal ALS function.  Uses Diagonal ALS form to estimate only
## the diagonal elements of Qw and Rv.
## Function Inputs : data.yk (measurements), data.uk (inputs),
## data.xhatk (state estimates- optional), data.datapts (number of data
## points considered), data.start (data to be ignored in the beginning 
## till initial condition is negligible),
## model.A, model.B (optional), model.C, model.G , N (window size),
## estimator.Q (Qw initial guess), estimator.R (Rv initial guess),
## estimator.L (initial estimator gain - optional).
## Function outputs : estimated Qw, estimated Rv, estimated filter gain
## L, and ALS LHS matrix As and RHS vector bhat.
## Please see file simulate_data8_diag.m for an example of the 
## Diagonal ALS implementation for a simple system.
## For Matlab implementation, the call of qp must be replaced by
## quadprog - see line 133
    
    
@function
def als_diag(data=None,N=None,model=None,estimator=None,*args,**kwargs):
    varargin = als_diag.varargin
    nargin = als_diag.nargin

    if nargin != 4:
        error('als: invalid number of arguments')
    
    
    datapts=data.datapts
# inst\als_diag.m:24
    Aa=model.A
# inst\als_diag.m:25
    Ca=model.C
# inst\als_diag.m:26
    Ga=model.G
# inst\als_diag.m:27
    n,g=size(Ga,nargout=2)
# inst\als_diag.m:28
    p=size(Ca,1)
# inst\als_diag.m:29
    if isfield(model,'B'):
        Ba=model.B
# inst\als_diag.m:30
        m=columns(Ba)
# inst\als_diag.m:30
    else:
        Ba=zeros(n)
# inst\als_diag.m:31
        m=copy(n)
# inst\als_diag.m:31
        data.uk = copy(zeros(m,datapts))
# inst\als_diag.m:31
    
    start=data.start
# inst\als_diag.m:33
    na=copy(n)
# inst\als_diag.m:35
    ga=copy(g)
# inst\als_diag.m:35
    pa=copy(p)
# inst\als_diag.m:35
    ## Estimator Simulation
    y=data.yk
# inst\als_diag.m:37
    u=data.uk
# inst\als_diag.m:38
    if isfield(estimator,'L'):
        L=estimator.L
# inst\als_diag.m:40
    else:
        L=dlqe(Aa,Ga,Ca,estimator.Q,estimator.R)
# inst\als_diag.m:41
    
    if (logical_not(isfield(data,'xhatk'))):
        xhat=zeros(na,datapts)
# inst\als_diag.m:45
        xhat_=zeros(n,datapts)
# inst\als_diag.m:46
        xhat_[arange(1,n),1]=model.xhat0
# inst\als_diag.m:47
        for i in arange(1,datapts).reshape(-1):
            xhat[arange(),i]=xhat_(arange(),i) + dot(L,(y(arange(),i) - dot(Ca,xhat_(arange(),i))))
# inst\als_diag.m:49
            xhat_[arange(),i + 1]=dot(Aa,xhat(arange(),i)) + dot(Ba,u(arange(),i))
# inst\als_diag.m:50
        xhat_=xhat_(arange(),arange(1,end() - 1))
# inst\als_diag.m:52
    else:
        xhat_=data.xhatk
# inst\als_diag.m:54
    
    inntrun=y(arange(),arange(start + 1,end())) - dot(Ca,xhat_(arange(),arange(start + 1,end())))
# inst\als_diag.m:57
    datatrun=datapts - start
# inst\als_diag.m:58
    ## Calculation of Autocorrelations for one column ALS
    Eyy=[]
# inst\als_diag.m:61
    for i in arange(0,N - 1).reshape(-1):
        temp=dot(inntrun(arange(),arange(i + 1,end())),inntrun(arange(),arange(1,end() - i)).T)
# inst\als_diag.m:63
        temp=temp / (datatrun - i)
# inst\als_diag.m:64
        Eyy=concat([[Eyy],[temp]])
# inst\als_diag.m:65
    
    Eyy=ravel(Eyy)
# inst\als_diag.m:67
    ##################################################
## Building the constant matrix for the LS problem
##################################################
    
    Ain=Aa - dot(dot(Aa,L),Ca)
# inst\als_diag.m:73
    OO=[]
# inst\als_diag.m:75
    temp=eye(n)
# inst\als_diag.m:76
    for i in arange(1,N).reshape(-1):
        OO=concat([[OO],[dot(Ca,temp)]])
# inst\als_diag.m:78
        temp=dot(temp,Ain)
# inst\als_diag.m:79
    
    ## temporary variables
    i=1
# inst\als_diag.m:83
    for j in arange(1,ga).reshape(-1):
        for k in arange(1,ga).reshape(-1):
            II=zeros(ga)
# inst\als_diag.m:86
            II[k,j]=1
# inst\als_diag.m:87
            t1=dlyap(Ain,dot(dot(Ga,II),Ga.T))
# inst\als_diag.m:88
            M1[arange(),i]=ravel(t1)
# inst\als_diag.m:89
            i=i + 1
# inst\als_diag.m:90
    
    i=1
# inst\als_diag.m:93
    for j in arange(1,pa).reshape(-1):
        for k in arange(1,pa).reshape(-1):
            II=zeros(pa)
# inst\als_diag.m:96
            II[k,j]=1
# inst\als_diag.m:97
            t2=dlyap(Ain,dot(dot(dot(dot(Aa,L),II),L.T),Aa.T))
# inst\als_diag.m:98
            M2[arange(),i]=ravel(t2)
# inst\als_diag.m:99
            i=i + 1
# inst\als_diag.m:100
    
    #############################
## Diagonal ALS method
#############################
    
    PSI=eye(pa)
# inst\als_diag.m:108
    for i in arange(1,N - 1).reshape(-1):
        PSI=concat([[PSI],[dot(dot(dot(- Ca,Ain ** (i - 1)),Aa),L)]])
# inst\als_diag.m:110
    
    
    OOtemp=kron(Ca,OO)
# inst\als_diag.m:113
    PSItemp=kron(eye(pa),PSI)
# inst\als_diag.m:114
    LH1=dot(OOtemp,M1)
# inst\als_diag.m:116
    LH2=dot(OOtemp,M2) + PSItemp
# inst\als_diag.m:117
    As_diag=concat([LH1(arange(),arange(1,ga ** 2,ga + 1)),LH2(arange(),arange(1,pa ** 2,pa + 1))])
# inst\als_diag.m:119
    # Testing the uniqueness of covariance estimates
    Arank=rank(As_diag,0.0001)
# inst\als_diag.m:122
    nr,nc=size(As_diag,nargout=2)
# inst\als_diag.m:124
    if nc > Arank:
        printf('Warning: Covariance estimates are not unique!\\n')
        pause
    
    Xest_diag=qp(ones(ga + pa,1),dot(As_diag.T,As_diag),dot(- As_diag.T,Eyy),[],[],zeros(ga + pa,1),[],[],eye(pa + ga),[])
# inst\als_diag.m:131
    # For Matlab, call line 134 instead of 131:
#Xest_diag = quadprog(As_diag'*As_diag, -As_diag'*Eyy, -eye(pa+ga), zeros(ga+pa,1));
    
    if prod(Xest_diag) == 0:
        printf('Warning: Covariance estimate(s) is (are) at constraints! You may have bad data! \\n')
    
    Qest=diag(Xest_diag(arange(1,ga)))
# inst\als_diag.m:140
    Rest=diag(Xest_diag(arange(ga + 1,end())))
# inst\als_diag.m:141
    Lest=ECM_iter(Aa,Ga,Ca,Qest,Rest)
# inst\als_diag.m:143
    As=copy(As_diag)
# inst\als_diag.m:145
    bhat=copy(Eyy)
# inst\als_diag.m:146