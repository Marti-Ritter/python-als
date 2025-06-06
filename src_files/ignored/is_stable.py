# Generated with SMOP  0.41-beta
from libsmop import *
# inst\is_stable.m

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
# @deftypefn {Function File} {} is_stable (@var{a}, @var{tol}, @var{dflg})
# @deftypefnx {Function File} {} is_stable (@var{sys}, @var{tol})
# Returns 1 if the matrix @var{a} or the system @var{sys}
# is stable, or 0 if not.
    
    # @strong{Inputs}
# @table @var
# @item  tol
# is a roundoff paramter, set to 200*@var{eps} if omitted.
# @item dflg
# Digital system flag (not required for system data structure):
# @table @code
# @item @var{dflg} != 0
# stable if eig(a) in unit circle
    
    # @item @var{dflg} == 0
# stable if eig(a) in open LHP (default)
# @end table
# @end table
# @end deftypefn
# @seealso{size, rows, columns, length, ismatrix, isscalar, isvector
# is_observable, is_stabilizable, is_detectable, krylov, and krylovb}
    
    # Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
# Created: August 1993
# Updated by John Ingram (ingraje@eng.auburn.edu) July, 1996 for systems
# Updated to simpler form by a.s.hodel 1998
    
    
@function
def is_stable(a=None,tol=None,disc=None,*args,**kwargs):
    varargin = is_stable.varargin
    nargin = is_stable.nargin

    if (nargin < 1) or (nargin > 3):
        usage('is_stable(a {,tol,disc})')
    else:
        if (isstruct(a)):
            # system was passed
            if nargin < 3:
                disc=is_digital(a)
# inst\is_stable.m:55
            else:
                if disc != is_digital(a):
                    warning('is_stable: disc =%d does not match system',disc)
            sys=sysupdate(a,'ss')
# inst\is_stable.m:59
            a=sys2ss(sys)
# inst\is_stable.m:60
        else:
            if (nargin < 3):
                disc=0
# inst\is_stable.m:63
            if (issquare(a) == 0):
                error('A(%dx%d) must be square',rows(A),columns(A))
    
    if (nargin < 2):
        tol=dot(200,eps)
# inst\is_stable.m:71
    else:
        if length(tol) != 1:
            error('is_stable: tol(%dx%d) must be a scalar',rows(tol),columns(tol))
    
    l=eig(a)
# inst\is_stable.m:76
    if disc:
        nbad=sum(dot(abs(l),(1 + tol)) > 1)
# inst\is_stable.m:78
    else:
        nbad=sum(real(l) + tol > 0)
# inst\is_stable.m:80
    
    retval=(nbad == 0)
# inst\is_stable.m:82
    return retval
    
if __name__ == '__main__':
    pass
    