# Generated with SMOP  0.41-beta
from libsmop import *
# inst\als_condno.m

    
@function
def als_condno(N=None,model=None,estimator=None,Rform=None,*args,**kwargs):
    varargin = als_condno.varargin
    nargin = als_condno.nargin

    ## Returns condition number and largest and smaller singular values for
## the ALS matrix (script A) formed from N, model, estimator.  Does not
## solve ALS problem
##
## [cond_no svals Mcond] = als_condno(N,model,estimator,Rform)
##
## Function Inputs :
## N (window size),
## model.A,
## model.C, 
## model.G (optional, default is identity matrix)
## estimator.L (initial estimator gain - optional)
## or estimator.Q and estimator.R (process and measurement noise
## covariances, used to calculate estimator gain),
## Rform (optional input to specify structure of R_v.  If
## Rform = "sym", R_v is symmetric, if Rform = "diag", R_v is
## diagonal.  Default is diagonal)
## 
## Function outputs : 
## cond_no (condition number of script A)
## svals (singular values of script A)
## Mcond (condition number of M = kron(C,I)*(I-kron(A,A))^-1*kron(G,G)*D_g
    
    Aa=model.A
# inst\als_condno.m:27
    Ca=model.C
# inst\als_condno.m:28
    if isfield(model,'G'):
        Ga=model.G
# inst\als_condno.m:29
    else:
        Ga=eye(size(Aa))
# inst\als_condno.m:30
    
    n,g=size(Ga,nargout=2)
# inst\als_condno.m:32
    p=size(Ca,1)
# inst\als_condno.m:33
    if isfield(estimator,'L'):
        L=estimator.L
# inst\als_condno.m:35
    else:
        L=dlqe(Aa,Ga,Ca,estimator.Q,estimator.R)
# inst\als_condno.m:36
    
    na=copy(n)
# inst\als_condno.m:39
    ga=copy(g)
# inst\als_condno.m:39
    pa=copy(p)
# inst\als_condno.m:39
    if nargin == 3:
        Rsym=0
# inst\als_condno.m:42
    else:
        if strmatch(lower(Rform),'sym'):
            Rsym=1
# inst\als_condno.m:44
        else:
            if strmatch(lower(Rform),'diag'):
                Rsym=0
# inst\als_condno.m:46
            else:
                warning('als_sdp_mrQ: Unknown structure type for R; defaulting to diagonal')
                Rsym=0
# inst\als_condno.m:49
    
    
    ##################################################
## Building the constant matrix for the LS problem
##################################################
    
    Ain=Aa - dot(dot(Aa,L),Ca)
# inst\als_condno.m:57
    OO=[]
# inst\als_condno.m:59
    temp=eye(n)
# inst\als_condno.m:60
    for i in arange(1,N).reshape(-1):
        OO=concat([[OO],[dot(Ca,temp)]])
# inst\als_condno.m:62
        temp=dot(temp,Ain)
# inst\als_condno.m:63
    
    M1=zeros(na ** 2,ga ** 2)
# inst\als_condno.m:65
    i=1
# inst\als_condno.m:66
    for j in arange(1,ga).reshape(-1):
        for k in arange(1,ga).reshape(-1):
            II=zeros(ga)
# inst\als_condno.m:69
            II[k,j]=1
# inst\als_condno.m:70
            t1=dlyap(Ain,dot(dot(Ga,II),Ga.T))
# inst\als_condno.m:71
            M1[arange(),i]=ravel(t1)
# inst\als_condno.m:72
            i=i + 1
# inst\als_condno.m:73
    
    M2=zeros(na ** 2,pa ** 2)
# inst\als_condno.m:76
    i=1
# inst\als_condno.m:77
    for j in arange(1,pa).reshape(-1):
        for k in arange(1,pa).reshape(-1):
            II=zeros(pa)
# inst\als_condno.m:80
            II[k,j]=1
# inst\als_condno.m:81
            t2=dlyap(Ain,dot(dot(dot(dot(Aa,L),II),L.T),Aa.T))
# inst\als_condno.m:82
            M2[arange(),i]=ravel(t2)
# inst\als_condno.m:83
            i=i + 1
# inst\als_condno.m:84
    
    
    ##########################################
## Single column ALS method
##########################################
    
    PSI=eye(pa)
# inst\als_condno.m:92
    for i in arange(1,N - 1).reshape(-1):
        PSI=concat([[PSI],[dot(dot(dot(- Ca,Ain ** (i - 1)),Aa),L)]])
# inst\als_condno.m:94
    
    
    OOtemp=kron(Ca,OO)
# inst\als_condno.m:97
    PSItemp=kron(eye(pa),PSI)
# inst\als_condno.m:98
    LH1=dot(OOtemp,M1)
# inst\als_condno.m:99
    LH2=dot(OOtemp,M2) + PSItemp
# inst\als_condno.m:100
    if Rsym == 1:
        LHSsingc=concat([dot(LH1,symtran(ga)),dot(LH2,symtran(pa))])
# inst\als_condno.m:103
    else:
        # Adds symmetric constraint to Q_w and diagonal constraint to R_v
        LHSsingc=concat([dot(LH1,symtran(ga)),LH2(arange(),arange(1,pa ** 2,pa + 1))])
# inst\als_condno.m:106
    
    M=dot(dot(kron(Ca,eye(n)),M1),symtran(ga))
# inst\als_condno.m:108
    Mcond=cond(M)
# inst\als_condno.m:109
    cond_no=cond(LHSsingc)
# inst\als_condno.m:111
    svals=svd(LHSsingc)
# inst\als_condno.m:112