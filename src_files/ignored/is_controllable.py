# Generated with SMOP  0.41-beta
from libsmop import *
# inst\is_controllable.m

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
# @deftypefn {Function File} {[@var{retval}, @var{u}] =} is_controllable (@var{sys}, @var{tol})
# @deftypefnx {Function File} {[@var{retval}, @var{u}] =} is_controllable (@var{a}, @var{b}, @var{tol})
# Logical check for system controllability.
    
    # @strong{Inputs}
# @table @var
# @item sys
# system data structure
# @item a
# @itemx b
# @var{n} by @var{n}, @var{n} by @var{m} matrices, respectively
# @item tol
# optional roundoff paramter.  default value: @code{10*eps}
# @end table
    
    # @strong{Outputs}
# @table @var
# @item retval
# Logical flag; returns true (1) if the system @var{sys} or the
# pair (@var{a},@var{b}) is controllable, whichever was passed as input
# arguments.
# @item U
#  U is an orthogonal basis of the controllable subspace.
# @end table
    
    # @strong{Method}
# Controllability is determined by applying Arnoldi iteration with
# complete re-orthogonalization to obtain an orthogonal basis of the
# Krylov subspace
# @example
# span ([b,a*b,...,a^@{n-1@}*b]).
# @end example
# The Arnoldi iteration is executed with @code{krylov} if the system
# has a single input; otherwise a block Arnoldi iteration is performed
# with @code{krylovb}.
# @end deftypefn
# @seealso{size, rows, columns, length, ismatrix, isscalar, isvector
# is_observable, is_stabilizable, is_detectable, krylov, and krylovb}
    
    # Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
# Created: August 1993
# Updated by A. S. Hodel (scotte@eng.auburn.edu) Aubust, 1995 to use krylovb
# Updated by John Ingram (ingraje@eng.auburn.edu) July, 1996 for packed systems
    
    
@function
def is_controllable(a=None,b=None,tol=None,*args,**kwargs):
    varargin = is_controllable.varargin
    nargin = is_controllable.nargin

    deftol=1
# inst\is_controllable.m:66
    
    if (nargin < logical_or(1,nargin) > 3):
        usage('[retval,U] = %s\n\t%s','is_controllable(a {, b, tol})','is_controllable(sys{,tol})')
    else:
        if (isstruct(a)):
            # system structure passed.
            sys=sysupdate(a,'ss')
# inst\is_controllable.m:72
            a,bs=sys2ss(sys,nargout=2)
# inst\is_controllable.m:73
            if (nargin > 2):
                usage('[retval,U] = is_controllable(sys{,tol})')
            else:
                if (nargin == 2):
                    tol=copy(b)
# inst\is_controllable.m:77
                    deftol=0
# inst\is_controllable.m:78
            b=copy(bs)
# inst\is_controllable.m:80
        else:
            # a,b arguments sent directly.
            if (nargin < 2):
                usage('[retval,U] = is_controllable(a {, b ,tol})')
            else:
                deftol=1
# inst\is_controllable.m:86
    
    # check for default tolerance
    if (deftol):
        tol=dot(1000,eps)
# inst\is_controllable.m:92
    
    # check tol dimensions
    if length(tol) != 1:
        error('is_controllable: tol(%dx%d) must be a scalar',rows(tol),columns(tol))
    else:
        if logical_not(is_sample(tol)):
            error('is_controllable: tol=%e must be positive',tol)
    
    # check dimensions compatibility
    n=issquare(a)
# inst\is_controllable.m:104
    nr,nc=size(b,nargout=2)
# inst\is_controllable.m:105
    if (n == 0 or n != nr or nc == 0):
        warning('is_controllable: a=(%dx%d), b(%dx%d)',rows(a),columns(a),nr,nc)
        retval=0
# inst\is_controllable.m:109
    else:
        # call block-krylov subspace routine to get an orthogonal basis
    # of the controllable subspace.
        U,__,Ucols=krylov(a,b,n,tol,1,nargout=3)
# inst\is_controllable.m:113
        retval=(Ucols == n)
# inst\is_controllable.m:114
    
    return retval,U
    
if __name__ == '__main__':
    pass
    