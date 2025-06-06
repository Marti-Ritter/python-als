# Generated with SMOP  0.41-beta
from libsmop import *
# inst\is_digital.m

    # Copyright (C) 1996, 1999 Auburn University.  All rights reserved.
    
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
# @deftypefn {Function File} {} is_digital (@var{sys})
# Return nonzero if system is digital;
# inputs:
# sys: system data structure
# eflg: 0 [default] exit with an error if system is mixed (continuous and
# discrete components)
#     : 1 print a warning if system is mixed (continuous and discrete)
#     : 2 silent operation
# outputs:
# DIGITAL:  0: system is purely continuous
#        :  1: system is purely discrete
#        : -1: system is mixed continuous and discrete
# Exits with an error of sys is a mixed (continuous and discrete) system
# @end deftypefn
    
    # Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
# Created: July 1996
    
    
@function
def is_digital(sys=None,eflg=None,*args,**kwargs):
    varargin = is_digital.varargin
    nargin = is_digital.nargin

    if (1) == (nargin):
        eflg=0
# inst\is_digital.m:41
    else:
        if (2) == (nargin):
            if (isempty(find(eflg == concat([0,1,2]),1))):
                error('invalid value of eflg=%d (%e)',eflg,eflg)
        else:
            usage('DIGITAL = is_digital(sys{,eflg})')
    
    # checked for sampled data system (mixed)
# discrete system
    sysyd=sysgetsignals(sys,'yd')
# inst\is_digital.m:52
    nn,nz=sysdimensions(sys,nargout=2)
# inst\is_digital.m:53
    cont=sum(sysyd == 0) + nn
# inst\is_digital.m:54
    tsam=sysgettsam(sys)
# inst\is_digital.m:55
    dig=sum(sysyd != 0) + nz + tsam
# inst\is_digital.m:56
    # check for mixed system
    if (dot(cont,dig) != 0):
        if (0) == (eflg):
            error('continuous/discrete system; use syscont, sysdisc, or c2d first')
        else:
            if (1) == (eflg):
                warning('is_digital: mixed continuous/discrete system')
        dig_sign=- 1
# inst\is_digital.m:66
    else:
        dig_sign=1
# inst\is_digital.m:68
    
    DIGITAL=dot(dig_sign,(tsam > 0))
# inst\is_digital.m:71
    return DIGITAL
    
if __name__ == '__main__':
    pass
    