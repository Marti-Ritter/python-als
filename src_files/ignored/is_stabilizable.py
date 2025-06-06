# Generated with SMOP  0.41-beta
from libsmop import *
# inst\is_stabilizable.m

    # Copyright (C) 1993, 1994, 1995 Auburn University.  All rights reserved.
    
    # This file is part of Octave.
    
    # Octave is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2, or (at your option) any
# later version.
    
    # Octave is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
    
    # You should have received a copy of the GNU General Public License
# along with Octave; see the file COPYING.  If not, write to the Free
# Software Foundation, 59 Temple Place, Suite 330, Boston, MA 02111 USA.
    
    # -*- texinfo -*-
# @deftypefn {Function File} {[@var{retval}, @var{u}] =} is_stabilizable (@var{sys}, @var{tol})
# @deftypefnx {Function File} {[@var{retval}, @var{u}] =} is_stabilizable (@var{a}, @var{b}, @var{tol})
# Logical check for system stabilizability (i.e., all unstable modes are controllable).
    
    # Test for stabilizability is performed via an ordered Schur decomposition
# that reveals the unstable subspace of the system @var{a} matrix.
    
    # Returns @code{retval} = 1 if the system, @var{a}, is stabilizable,
# if the pair  (@var{a}, @var{b}) is stabilizable, or 0 if not.
# @var{u} = orthogonal basis of controllable subspace.
    
    # Controllable subspace is determined by applying Arnoldi iteration with
# complete re-orthogonalization to obtain an orthogonal basis of the
# Krylov subspace.
# @example
#   span ([b,a*b,...,a^   b]).
# @end example
# tol is a roundoff paramter, set to 200*eps if omitted.
# @end deftypefn
    
    # See also: size, rows, columns, length, ismatrix, isscalar, isvector
#     is_observable, is_stabilizable, is_detectable
    
    # Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
# Created: August 1993
# Updated by A. S. Hodel (scotte@eng.auburn.edu) Aubust, 1995 to use krylovb
# Updated by John Ingram (ingraje@eng.auburn.edu) July, 1996 to accept systems
    
    
@function
def is_stabilizable(a=None,b=None,tol=None,*args,**kwargs):
    varargin = is_stabilizable.varargin
    nargin = is_stabilizable.nargin

    if (nargin < 1):
        usage('[retval,U] = is_stabilizable(a {, b ,tol})')
    else:
        if (isstruct(a)):
            # sustem passed.
            if (nargin == 2):
                tol=copy(b)
# inst\is_stabilizable.m:55
            else:
                if (nargin > 2):
                    usage('[retval,U] = is_stabilizable(sys{,tol})')
            a,b=sys2ss(a,nargout=2)
# inst\is_stabilizable.m:59
        else:
            # a,b arguments sent directly.
            if (nargin > 3):
                usage('[retval,U] = is_stabilizable(a {, b ,tol})')
    
    if exist('tol'):
        retval,U=is_controllable(a,b,tol,nargout=2)
# inst\is_stabilizable.m:68
    else:
        retval,U=is_controllable(a,b,nargout=2)
# inst\is_stabilizable.m:70
        tol=dot(dot(100.0,rows(b)),eps)
# inst\is_stabilizable.m:71
    
    if logical_not(retval) and size(U,2) > 0:
        # now use an ordered Schur decomposition to get an orthogonal
    # basis of the unstable subspace...
        n=rows(a)
# inst\is_stabilizable.m:77
        ua,__=schur(- (a + dot(eye(n),tol)),'A',nargout=2)
# inst\is_stabilizable.m:78
        k=sum(real(eig(a)) >= 0)
# inst\is_stabilizable.m:79
        if (k > 0):
            ua=ua(arange(),arange(1,k))
# inst\is_stabilizable.m:82
            retval=(norm(ua - dot(dot(U,U.T),ua)) < tol)
# inst\is_stabilizable.m:84
        else:
            retval=1
# inst\is_stabilizable.m:86
    
    return retval,U
    
if __name__ == '__main__':
    pass
    