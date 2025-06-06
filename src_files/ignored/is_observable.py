# Generated with SMOP  0.41-beta
from libsmop import *
# inst\is_observable.m

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
# @deftypefn {Function File} {[@var{retval}, @var{u}] =} is_observable (@var{a}, @var{c}, @var{tol})
# @deftypefnx {Function File} {[@var{retval}, @var{u}] =} is_observable (@var{sys}, @var{tol})
# Logical check for system observability.
    
    # Default: tol = 10*norm(a,'fro')*eps
    
    # Returns 1 if the system @var{sys} or the pair (@var{a},@var{c}) is
# observable, 0 if not.
    
    # @strong{See} @code{is_controllable} for detailed description of arguments
# and default values.
# @end deftypefn
# @seealso{size, rows, columns, length, ismatrix, isscalar, and isvector}
    
    # Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
# Created: August 1993
# Updated by John Ingram (ingraje@eng.auburn.edu) July 1996.
    
    
@function
def is_observable(a=None,c=None,tol=None,*args,**kwargs):
    varargin = is_observable.varargin
    nargin = is_observable.nargin

    if (nargin < 1):
        usage('[retval,U] = is_observable(a , c {, tol})')
    else:
        if (isstruct(a)):
            # system form
            if (nargin == 2):
                tol=copy(c)
# inst\is_observable.m:45
            else:
                if (nargin > 2):
                    usage('[retval,U] = is_observable(sys {, tol})')
            a,__,c=sys2ss(a,nargout=3)
# inst\is_observable.m:49
        else:
            if (nargin > 3):
                usage('[retval,U] = is_observable(a , c {, tol})')
    
    if (exist('tol')):
        retval,U=is_controllable(a.T,c.T,tol,nargout=2)
# inst\is_observable.m:54
    else:
        retval,U=is_controllable(a.T,c.T,nargout=2)
# inst\is_observable.m:56
    
    return retval,U
    
if __name__ == '__main__':
    pass
    