# Generated with SMOP  0.41-beta
from libsmop import *
# inst\is_sample.m

    # Copyright (C) 1996 Auburn University.  All rights reserved.
    
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
# @deftypefn {Function File} {} is_sample (@var{ts})
# Return true if @var{ts} is a valid sampling time
# (real,scalar, > 0)
# @end deftypefn
    
    # Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
# Created: July 1995
    
    
@function
def is_sample(Ts=None,*args,**kwargs):
    varargin = is_sample.varargin
    nargin = is_sample.nargin

    out=(isscalar(Ts) and (Ts == abs(Ts)) and (Ts != 0))
# inst\is_sample.m:30
    return out
    
if __name__ == '__main__':
    pass
    