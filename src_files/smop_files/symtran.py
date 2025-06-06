# Generated with SMOP  0.41-beta
from libsmop import *
# als\inst\symtran.m

    ## Function to give the duplication matrix that will impose symmetry constraints on
## the covariance matrices, (Q)_s = tran*(Q)_ss
    
@function
def symtran(n=None,*args,**kwargs):
    varargin = symtran.varargin
    nargin = symtran.nargin

    r=dot(n,(n + 1)) / 2
# als\inst\symtran.m:5
    tran=zeros(n ** 2,r)
# als\inst\symtran.m:7
    k=1
# als\inst\symtran.m:8
    for i in arange(1,n).reshape(-1):
        for j in arange(i,n).reshape(-1):
            temp=zeros(n)
# als\inst\symtran.m:11
            temp[i,j]=1
# als\inst\symtran.m:12
            if i == j:
                div=2
# als\inst\symtran.m:13
            else:
                div=1
# als\inst\symtran.m:13
            t2=(temp + temp.T) / div
# als\inst\symtran.m:14
            tran[arange(),k]=ravel(t2)
# als\inst\symtran.m:15
            k=k + 1
# als\inst\symtran.m:16
    