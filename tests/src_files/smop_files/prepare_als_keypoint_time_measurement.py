# Generated with SMOP  0.41-beta
from libsmop import *
# als/prepare_als_keypoint_time_measurement.m

    
@function
def prepare_als_keypoint_time_measurement(n_keypoints=None,n_dimensions=None,datapts=None,*args,**kwargs):
    varargin = prepare_als_keypoint_time_measurement.varargin
    nargin = prepare_als_keypoint_time_measurement.nargin

    base_path='./sample_matrices/'
# als/prepare_als_keypoint_time_measurement.m:2
    file_end=strcat('_k',num2str(n_keypoints),'_d',num2str(n_dimensions),'.csv')
# als/prepare_als_keypoint_time_measurement.m:3
    my_F=csvread(strcat(base_path,'F',file_end))(arange(2,end()),arange(2,end()))
# als/prepare_als_keypoint_time_measurement.m:5
    my_H=csvread(strcat(base_path,'H',file_end))(arange(2,end()),arange(2,end()))
# als/prepare_als_keypoint_time_measurement.m:6
    my_Q=csvread(strcat(base_path,'Q',file_end))(arange(2,end()),arange(2,end()))
# als/prepare_als_keypoint_time_measurement.m:7
    my_R=csvread(strcat(base_path,'R',file_end))(arange(2,end()),arange(2,end()))
# als/prepare_als_keypoint_time_measurement.m:8
    
    randn('seed',100)
    Aa=copy(my_F)
# als/prepare_als_keypoint_time_measurement.m:13
    Ca=copy(my_H)
# als/prepare_als_keypoint_time_measurement.m:14
    Q_w=copy(my_Q)
# als/prepare_als_keypoint_time_measurement.m:16
    R_v=copy(my_R)
# als/prepare_als_keypoint_time_measurement.m:17
    vec_Qw,eig_Qw=eig(Q_w,nargout=2)
# als/prepare_als_keypoint_time_measurement.m:18
    vec_Rv,eig_Rv=eig(R_v,nargout=2)
# als/prepare_als_keypoint_time_measurement.m:19
    mult_Qw=dot(vec_Qw,sqrt(eig_Qw))
# als/prepare_als_keypoint_time_measurement.m:20
    mult_Rv=dot(vec_Rv,sqrt(eig_Rv))
# als/prepare_als_keypoint_time_measurement.m:21
    
    Qw_hat=dot(0.01,Q_w)
# als/prepare_als_keypoint_time_measurement.m:24
    Rv_hat=dot(0.001,R_v)
# als/prepare_als_keypoint_time_measurement.m:25
    pa,na=size(Ca,nargout=2)
# als/prepare_als_keypoint_time_measurement.m:27
    n=copy(na)
# als/prepare_als_keypoint_time_measurement.m:29
    p=copy(pa)
# als/prepare_als_keypoint_time_measurement.m:30
    L=dlqe(Aa,[],Ca,Qw_hat,Rv_hat)
# als/prepare_als_keypoint_time_measurement.m:32
    xhat=zeros(na,datapts)
# als/prepare_als_keypoint_time_measurement.m:33
    xhat_=zeros(na,datapts)
# als/prepare_als_keypoint_time_measurement.m:34
    x[arange(),1]=dot(10,ones(na,1))
# als/prepare_als_keypoint_time_measurement.m:36
    
    xhat_[arange(1,na),1]=x(arange(),1)
# als/prepare_als_keypoint_time_measurement.m:37
    
    for i in arange(1,datapts).reshape(-1):
        y[arange(),i]=dot(Ca,x(arange(),i)) + dot(mult_Rv,randn(pa,1))
# als/prepare_als_keypoint_time_measurement.m:40
        xhat[arange(),i]=xhat_(arange(),i) + dot(L,(y(arange(),i) - dot(Ca,xhat_(arange(),i))))
# als/prepare_als_keypoint_time_measurement.m:41
        x[arange(),i + 1]=dot(Aa,x(arange(),i))
# als/prepare_als_keypoint_time_measurement.m:42
        xhat_[arange(),i + 1]=dot(Aa,xhat(arange(),i))
# als/prepare_als_keypoint_time_measurement.m:43
    
    #########################
  ### SETUP ALS PROBLEM ###
  #########################
    
    model.A = copy(my_F)
# als/prepare_als_keypoint_time_measurement.m:50
    model.C = copy(my_H)
# als/prepare_als_keypoint_time_measurement.m:51
    model.xhat0 = copy(xhat_(arange(),1))
# als/prepare_als_keypoint_time_measurement.m:52
    data.datapts = copy(datapts)
# als/prepare_als_keypoint_time_measurement.m:54
    data.yk = copy(y)
# als/prepare_als_keypoint_time_measurement.m:55
    data.start = copy(100)
# als/prepare_als_keypoint_time_measurement.m:56
    estimator.Q = copy(my_Q)
# als/prepare_als_keypoint_time_measurement.m:58
    estimator.R = copy(my_R)
# als/prepare_als_keypoint_time_measurement.m:59
    return data,model,estimator
    
if __name__ == '__main__':
    pass
    