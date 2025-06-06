# Generated with SMOP  0.41-beta
from libsmop import *
# als/measure_als_keypoint_time.m

    
@function
def measure_als_keypoint_time(data=None,N=None,model=None,estimator=None,rho_vec=None,*args,**kwargs):
    varargin = measure_als_keypoint_time.varargin
    nargin = measure_als_keypoint_time.nargin

    t=copy(cputime)
# als/measure_als_keypoint_time.m:2
    Qest_cell,Rest_cell,trQ,Phi2=als_sdp_mrQ(data,N,model,estimator,'rho_values',rho_vec,nargout=4)
# als/measure_als_keypoint_time.m:3
    seconds_taken=cputime - t
# als/measure_als_keypoint_time.m:4
    close_('all')
    return seconds_taken
    
if __name__ == '__main__':
    pass
    