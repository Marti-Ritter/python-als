# Generated with SMOP  0.41-beta
from libsmop import *
# als\inst\symtranT.m

    ## Function to "undo" the symmetry imposing of the duplication matrix
## tran*(x)s = (x)ss
## Equivalent to pinv(symtran(n)), but much faster
    
@function
def symtranT(n=None,*args,**kwargs):
    varargin = symtranT.varargin
    nargin = symtranT.nargin

    tran=symtran(n).T
# als\inst\symtranT.m:5
    for i in arange(1,size(tran,1),1).reshape(-1):
        tran[i,arange()]=tran(i,arange()) / sum(tran(i,arange()))
# als\inst\symtranT.m:7
    
    return tran
    
if __name__ == '__main__':
    pass
    
    
    