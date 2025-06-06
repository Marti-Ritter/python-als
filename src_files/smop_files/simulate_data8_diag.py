# Generated with SMOP  0.41-beta
from libsmop import *
# inst\simulate_data8_diag.m

    ## Diagonal ALS Implementation Example
    
    clear('all')
    randn('seed',100)
    Aa=diag(concat([0.1,0.2,0.3]))
# inst\simulate_data8_diag.m:6
    Aa[1,3]=0.1
# inst\simulate_data8_diag.m:6
    #Ga = [1;0.2;0.3];
    Ga=eye(3)
# inst\simulate_data8_diag.m:9
    #Ca = [0.1 0.2 0];
    Ca=eye(3)
# inst\simulate_data8_diag.m:12
    Q_w=diag(concat([0.5,0.2,0.1]))
# inst\simulate_data8_diag.m:14
    #Q_w = 0.5;
#R_v = 0.1;
    R_v=diag(concat([0.5,0.2,0.8]))
# inst\simulate_data8_diag.m:17
    S=dlyap(Aa,dot(dot(Ga,Q_w),Ga.T))
# inst\simulate_data8_diag.m:19
    vec_Qw,eig_Qw=eig(Q_w,nargout=2)
# inst\simulate_data8_diag.m:21
    vec_Rv,eig_Rv=eig(R_v,nargout=2)
# inst\simulate_data8_diag.m:22
    mult_Qw=dot(vec_Qw,sqrt(eig_Qw))
# inst\simulate_data8_diag.m:23
    mult_Rv=dot(vec_Rv,sqrt(eig_Rv))
# inst\simulate_data8_diag.m:24
    ### initial guesses
    
    G_hat=eye(3)
# inst\simulate_data8_diag.m:28
    Qw_hat=diag(concat([1,2,3]))
# inst\simulate_data8_diag.m:29
    #Qw_hat= 5*Q_w;
    Rv_hat=dot(0.001,R_v)
# inst\simulate_data8_diag.m:31
    pa,na=size(Ca,nargout=2)
# inst\simulate_data8_diag.m:33
    na,ga=size(Ga,nargout=2)
# inst\simulate_data8_diag.m:34
    n=copy(na)
# inst\simulate_data8_diag.m:36
    p=copy(pa)
# inst\simulate_data8_diag.m:37
    g=copy(ga)
# inst\simulate_data8_diag.m:38
    datapts=5000
# inst\simulate_data8_diag.m:40
    L=dlqe(Aa,G_hat,Ca,Qw_hat,Rv_hat)
# inst\simulate_data8_diag.m:42
    #[L,P]=dlqe(Aa,Ga,Ca,Q_w,R_v);
#L = zeros(n,p);
    
    P=dlyap((Aa - dot(dot(Aa,L),Ca)),dot(dot(concat([Ga,dot(- Aa,L)]),concat([[Q_w,zeros(g,p)],[zeros(p,g),R_v]])),concat([Ga,dot(- Aa,L)]).T))
# inst\simulate_data8_diag.m:46
    xhat=zeros(na,datapts)
# inst\simulate_data8_diag.m:48
    xhat_=zeros(na,datapts)
# inst\simulate_data8_diag.m:49
    x[arange(),1]=dot(10,ones(na,1))
# inst\simulate_data8_diag.m:51
    
    xhat_[arange(1,na),1]=x(arange(),1)
# inst\simulate_data8_diag.m:53
    
    for i in arange(1,datapts).reshape(-1):
        y[arange(),i]=dot(Ca,x(arange(),i)) + dot(mult_Rv,randn(pa,1))
# inst\simulate_data8_diag.m:57
        xhat[arange(),i]=xhat_(arange(),i) + dot(L,(y(arange(),i) - dot(Ca,xhat_(arange(),i))))
# inst\simulate_data8_diag.m:58
        x[arange(),i + 1]=dot(Aa,x(arange(),i)) + dot(Ga,(dot(mult_Qw,randn(ga,1))))
# inst\simulate_data8_diag.m:59
        xhat_[arange(),i + 1]=dot(Aa,xhat(arange(),i))
# inst\simulate_data8_diag.m:60
    
    #########################
### SETUP ALS PROBLEM ###
#########################
    
    model.A = copy(Aa)
# inst\simulate_data8_diag.m:68
    model.C = copy(Ca)
# inst\simulate_data8_diag.m:69
    model.G = copy(G_hat)
# inst\simulate_data8_diag.m:70
    model.xhat0 = copy(xhat_(arange(),1))
# inst\simulate_data8_diag.m:71
    data.datapts = copy(datapts)
# inst\simulate_data8_diag.m:73
    data.yk = copy(y)
# inst\simulate_data8_diag.m:74
    data.xhatk = copy(xhat_(arange(),arange(1,end() - 1)))
# inst\simulate_data8_diag.m:75
    data.start = copy(100)
# inst\simulate_data8_diag.m:76
    N=15
# inst\simulate_data8_diag.m:78
    estimator.L = copy(L)
# inst\simulate_data8_diag.m:80
    #estimator.Q = Qw_hat;
#estimator.R = Rv_hat;
    
    Qest,Rest,Lest,As,bhat=als_diag(data,N,model,estimator,nargout=5)
# inst\simulate_data8_diag.m:84