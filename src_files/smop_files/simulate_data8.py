# Generated with SMOP  0.41-beta
from libsmop import *
# inst\simulate_data8.m

    ## ALS Implementation Example
    
    clear('all')
    randn('seed',100)
    Aa=diag(concat([0.1,0.2,0.3]))
# inst\simulate_data8.m:6
    Aa[1,3]=0.1
# inst\simulate_data8.m:6
    Aa[1,2]=0.1
# inst\simulate_data8.m:6
    #Ga = [1;0.2;0.3];
    Ga=eye(3)
# inst\simulate_data8.m:10
    #Ca = [0.1 0.2 0];
    Ca=eye(3)
# inst\simulate_data8.m:13
    Q_w=diag(concat([0.5,0.2,0.1]))
# inst\simulate_data8.m:15
    #Q_w = 0.5;
#R_v = 0.1;
    R_v=diag(concat([0.5,0.2,0.8]))
# inst\simulate_data8.m:18
    S=dlyap(Aa,dot(dot(Ga,Q_w),Ga.T))
# inst\simulate_data8.m:20
    vec_Qw,eig_Qw=eig(Q_w,nargout=2)
# inst\simulate_data8.m:22
    vec_Rv,eig_Rv=eig(R_v,nargout=2)
# inst\simulate_data8.m:23
    mult_Qw=dot(vec_Qw,sqrt(eig_Qw))
# inst\simulate_data8.m:24
    mult_Rv=dot(vec_Rv,sqrt(eig_Rv))
# inst\simulate_data8.m:25
    ### initial guesses
    
    G_hat=eye(3)
# inst\simulate_data8.m:29
    Qw_hat=diag(concat([1,2,3]))
# inst\simulate_data8.m:30
    #Qw_hat= 5*Q_w;
    Rv_hat=dot(0.001,R_v)
# inst\simulate_data8.m:32
    pa,na=size(Ca,nargout=2)
# inst\simulate_data8.m:34
    na,ga=size(Ga,nargout=2)
# inst\simulate_data8.m:35
    n=copy(na)
# inst\simulate_data8.m:37
    p=copy(pa)
# inst\simulate_data8.m:38
    g=copy(ga)
# inst\simulate_data8.m:39
    datapts=5000
# inst\simulate_data8.m:41
    L=dlqe(Aa,G_hat,Ca,Qw_hat,Rv_hat)
# inst\simulate_data8.m:43
    #[L,P]=dlqe(Aa,Ga,Ca,Q_w,R_v);
#L = zeros(n,p);
    
    P=dlyap((Aa - dot(dot(Aa,L),Ca)),dot(dot(concat([Ga,dot(- Aa,L)]),concat([[Q_w,zeros(g,p)],[zeros(p,g),R_v]])),concat([Ga,dot(- Aa,L)]).T))
# inst\simulate_data8.m:47
    xhat=zeros(na,datapts)
# inst\simulate_data8.m:49
    xhat_=zeros(na,datapts)
# inst\simulate_data8.m:50
    x[arange(),1]=dot(10,ones(na,1))
# inst\simulate_data8.m:52
    
    xhat_[arange(1,na),1]=x(arange(),1)
# inst\simulate_data8.m:54
    
    for i in arange(1,datapts).reshape(-1):
        y[arange(),i]=dot(Ca,x(arange(),i)) + dot(mult_Rv,randn(pa,1))
# inst\simulate_data8.m:58
        xhat[arange(),i]=xhat_(arange(),i) + dot(L,(y(arange(),i) - dot(Ca,xhat_(arange(),i))))
# inst\simulate_data8.m:59
        x[arange(),i + 1]=dot(Aa,x(arange(),i)) + dot(Ga,(dot(mult_Qw,randn(ga,1))))
# inst\simulate_data8.m:60
        xhat_[arange(),i + 1]=dot(Aa,xhat(arange(),i))
# inst\simulate_data8.m:61
    
    #########################
### SETUP ALS PROBLEM ###
#########################
    
    model.A = copy(Aa)
# inst\simulate_data8.m:69
    model.C = copy(Ca)
# inst\simulate_data8.m:70
    model.G = copy(G_hat)
# inst\simulate_data8.m:71
    model.xhat0 = copy(xhat_(arange(),1))
# inst\simulate_data8.m:72
    data.datapts = copy(datapts)
# inst\simulate_data8.m:74
    data.yk = copy(y)
# inst\simulate_data8.m:75
    # data.xhatk = xhat_(:,1:end-1);
    data.start = copy(100)
# inst\simulate_data8.m:77
    N=15
# inst\simulate_data8.m:79
    estimator.L = copy(L)
# inst\simulate_data8.m:81
    #estimator.Q = Qw_hat;
#estimator.R = Rv_hat;
    
    # Using updated codes:
# This is in general how to call the ALS method:
    Qest_cell,Rest_cell=als_sdp_mrQ(data,N,model,estimator,nargout=2)
# inst\simulate_data8.m:87
    Qest1=Qest_cell[1]
# inst\simulate_data8.m:88
    Rest1=Rest_cell[1]
# inst\simulate_data8.m:89
    # Without semidefinite constraints
# you usually should keep the semidefinite constraints, but could remove them to see how good your model is
    Qest_celln,Rest_celln=als_sdp_mrQ(data,N,model,estimator,'sdp',0,'plot',0,nargout=2)
# inst\simulate_data8.m:92
    Qest_indef=Qest_celln[1]
# inst\simulate_data8.m:93
    Rest_indef=Rest_celln[1]
# inst\simulate_data8.m:94
    # With identity weighting
# Only use identity weighting if you have a good reason; the default of data-based weighting gives lower variance estimates
    Qest_celli,Rest_celli=als_sdp_mrQ(data,N,model,estimator,'weight','I','plot',0,nargout=2)
# inst\simulate_data8.m:97
    QestI=Qest_celli[1]
# inst\simulate_data8.m:98
    RestI=Rest_celli[1]
# inst\simulate_data8.m:99
    # Tradeoff curve example
# This example shows what to do if you have fewer outputs than states
# here we pretend that we only know the first row of C
# we generate a tradeoff curve and look for the elbow in the curve
    data2=copy(data)
# inst\simulate_data8.m:104
    data2.yk = copy(y(1,arange()))
# inst\simulate_data8.m:105
    model2=copy(model)
# inst\simulate_data8.m:106
    model2.C = copy(model.C(1,arange()))
# inst\simulate_data8.m:107
    rho_vec=logspace(- 6,6,25)
# inst\simulate_data8.m:108
    estimator2.L = copy(dlqe(Aa,G_hat,Ca(1,arange()),Qw_hat,Rv_hat(1,1)))
# inst\simulate_data8.m:109
    t=copy(cputime)
# inst\simulate_data8.m:111
    Qest_cell2,Rest_cell2,trQ,Phi2=als_sdp_mrQ(data2,N,model2,estimator2,'rho_values',rho_vec,nargout=4)
# inst\simulate_data8.m:112
    printf('Total cpu time: %f seconds\\n',cputime - t)
    figure
    plot(Phi2,trQ,'-o',Phi2(9),trQ(9),'r*')
    Qest2=Qest_cell2[9]
# inst\simulate_data8.m:116
    Rest2=Rest_cell2[9]
# inst\simulate_data8.m:117