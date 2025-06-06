# Generated with SMOP  0.41-beta
from libsmop import *
# inst\qp.m

    
    # Copyright (C) 2000-2021 The Octave Project Developers
    
    # See the file COPYRIGHT.md in the top-level directory of this
# distribution or <https://octave.org/copyright/>.
    
    # This file is part of Octave.
    
    # Octave is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
    
    # Octave is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
    
    # You should have received a copy of the GNU General Public License
# along with Octave; see the file COPYING.  If not, see
# <https://www.gnu.org/licenses/>.
    
    # -*- texinfo -*-
# @deftypefn  {} {[@var{x}, @var{obj}, @var{info}, @var{lambda}] =} qp (@var{x0}, @var{H})
# @deftypefnx {} {[@var{x}, @var{obj}, @var{info}, @var{lambda}] =} qp (@var{x0}, @var{H}, @var{q})
# @deftypefnx {} {[@var{x}, @var{obj}, @var{info}, @var{lambda}] =} qp (@var{x0}, @var{H}, @var{q}, @var{A}, @var{b})
# @deftypefnx {} {[@var{x}, @var{obj}, @var{info}, @var{lambda}] =} qp (@var{x0}, @var{H}, @var{q}, @var{A}, @var{b}, @var{lb}, @var{ub})
# @deftypefnx {} {[@var{x}, @var{obj}, @var{info}, @var{lambda}] =} qp (@var{x0}, @var{H}, @var{q}, @var{A}, @var{b}, @var{lb}, @var{ub}, @var{A_lb}, @var{A_in}, @var{A_ub})
# @deftypefnx {} {[@var{x}, @var{obj}, @var{info}, @var{lambda}] =} qp (@dots{}, @var{options})
# Solve a quadratic program (QP).
    
    # Solve the quadratic program defined by
# @tex
# $$
#  \min_x {1 \over 2} x^T H x + x^T q
# $$
# @end tex
# @ifnottex
    
    # @example
# @group
# min 0.5 x'*H*x + x'*q
#  x
# @end group
# @end example
    
    # @end ifnottex
# subject to
# @tex
# $$
#  A x = b \qquad lb \leq x \leq ub \qquad A_{lb} \leq A_{in} x \leq A_{ub}
# $$
# @end tex
# @ifnottex
    
    # @example
# @group
# A*x = b
# lb <= x <= ub
# A_lb <= A_in*x <= A_ub
# @end group
# @end example
    
    # @end ifnottex
# @noindent
# using a null-space active-set method.
    
    # Any bound (@var{A}, @var{b}, @var{lb}, @var{ub}, @var{A_in}, @var{A_lb},
# @var{A_ub}) may be set to the empty matrix (@code{[]}) if not present.  The
# constraints @var{A} and @var{A_in} are matrices with each row representing
# a single constraint.  The other bounds are scalars or vectors depending on
# the number of constraints.  The algorithm is faster if the initial guess is
# feasible.
    
    # @var{options} is a structure specifying additional parameters which
# control the algorithm.  Currently, @code{qp} recognizes these options:
# @qcode{"MaxIter"}, @qcode{"TolX"}.
    
    # @qcode{"MaxIter"} proscribes the maximum number of algorithm iterations
# before optimization is halted.  The default value is 200.
# The value must be a positive integer.
    
    # @qcode{"TolX"} specifies the termination tolerance for the unknown variables
# @var{x}.  The default is @code{sqrt (eps)} or approximately 1e-8.
    
    # On return, @var{x} is the location of the minimum and @var{fval} contains
# the value of the objective function at @var{x}.
    
    # @table @var
# @item info
# Structure containing run-time information about the algorithm.  The
# following fields are defined:
    
    # @table @code
# @item solveiter
# The number of iterations required to find the solution.
    
    # @item info
# An integer indicating the status of the solution.
    
    # @table @asis
# @item 0
# The problem is feasible and convex.  Global solution found.
    
    # @item 1
# The problem is not convex.  Local solution found.
    
    # @item 2
# The problem is not convex and unbounded.
    
    # @item 3
# Maximum number of iterations reached.
    
    # @item 6
# The problem is infeasible.
# @end table
# @end table
# @end table
# @seealso{sqp}
# @end deftypefn
    
    # PKG_ADD: ## Discard result to avoid polluting workspace with ans at startup.
# PKG_ADD: [~] = __all_opts__ ("qp");
    
    
@function
def qp(x0=None,H=None,varargin=None,*args,**kwargs):
    varargin = qp.varargin
    nargin = qp.nargin

    if (nargin == 1 and ischar(x0) and strcmp(x0,'defaults')):
        x=struct('MaxIter',200,'TolX',sqrt(eps))
# inst\qp.m:129
        return x,obj,INFO,lambda_
    
    nargs=copy(nargin)
# inst\qp.m:133
    if (nargs > 2 and isstruct(varargin[end()])):
        options=varargin[end()]
# inst\qp.m:135
        nargs -= 1
    else:
        options=struct()
# inst\qp.m:138
    
    if (nargs != 2 and nargs != 3 and nargs != 5 and nargs != 7 and nargs != 10):
        print_usage()
    
    if (nargs >= 3):
        q=varargin[1]
# inst\qp.m:146
    else:
        q=[]
# inst\qp.m:148
    
    if (nargs >= 5):
        A=varargin[2]
# inst\qp.m:152
        b=varargin[3]
# inst\qp.m:153
    else:
        A=[]
# inst\qp.m:155
        b=[]
# inst\qp.m:156
    
    if (nargs >= 7):
        lb=varargin[4]
# inst\qp.m:160
        ub=varargin[5]
# inst\qp.m:161
    else:
        lb=[]
# inst\qp.m:163
        ub=[]
# inst\qp.m:164
    
    if (nargs == 10):
        A_lb=varargin[6]
# inst\qp.m:168
        A_in=varargin[7]
# inst\qp.m:169
        A_ub=varargin[8]
# inst\qp.m:170
    else:
        A_lb=[]
# inst\qp.m:172
        A_in=[]
# inst\qp.m:173
        A_ub=[]
# inst\qp.m:174
    
    maxit=optimget(options,'MaxIter',200)
# inst\qp.m:177
    tol=optimget(options,'TolX',sqrt(eps))
# inst\qp.m:178
    
    if (logical_not(issquare(H))):
        error('qp: quadratic penalty matrix must be square')
    else:
        if (logical_not(ishermitian(H))):
            # warning ("qp: quadratic penalty matrix not hermitian");
            H=(H + H.T) / 2
# inst\qp.m:185
    
    n=rows(H)
# inst\qp.m:187
    
    # If empty it is resized to the right dimension and filled with 0.
    if (isempty(x0)):
        x0=zeros(n,1)
# inst\qp.m:192
    else:
        if (logical_not(isvector(x0))):
            error('qp: the initial guess X0 must be a vector')
        else:
            if (numel(x0) != n):
                error('qp: the initial guess X0 has incorrect length')
        x0=ravel(x0)
# inst\qp.m:199
    
    # Validate linear penalty.
    if (isempty(q)):
        q=zeros(n,1)
# inst\qp.m:204
    else:
        if (logical_not(isvector(q))):
            error('qp: Q must be a vector')
        else:
            if (numel(q) != n):
                error('qp: Q has incorrect length')
        q=ravel(q)
# inst\qp.m:211
    
    # Validate equality constraint matrices.
    if (isempty(A) or isempty(b)):
        A=zeros(0,n)
# inst\qp.m:216
        b=zeros(0,1)
# inst\qp.m:217
        n_eq=0
# inst\qp.m:218
    else:
        n_eq,n1=size(A,nargout=2)
# inst\qp.m:220
        if (n1 != n):
            error('qp: equality constraint matrix has incorrect column dimension')
        if (numel(b) != n_eq):
            error('qp: equality constraint matrix and vector have inconsistent dimensions')
    
    # Validate bound constraints.
    Ain=zeros(0,n)
# inst\qp.m:230
    bin=zeros(0,1)
# inst\qp.m:231
    n_in=0
# inst\qp.m:232
    if (nargs > 5):
        if (logical_not(isempty(lb))):
            if (numel(lb) != n):
                error('qp: lower bound LB has incorrect length')
            else:
                if (isempty(ub)):
                    Ain=concat([[Ain],[eye(n)]])
# inst\qp.m:238
                    bin=concat([[bin],[lb]])
# inst\qp.m:239
        if (logical_not(isempty(ub))):
            if (numel(ub) != n):
                error('qp: upper bound UB has incorrect length')
            else:
                if (isempty(lb)):
                    Ain=concat([[Ain],[- eye(n)]])
# inst\qp.m:247
                    bin=concat([[bin],[- ub]])
# inst\qp.m:248
        if (logical_not(isempty(lb)) and logical_not(isempty(ub))):
            rtol=copy(tol)
# inst\qp.m:253
            for i in arange(1,n).reshape(-1):
                if (abs(lb(i) - ub(i)) < dot(rtol,(1 + max(abs(lb(i) + ub(i)))))):
                    # These are actually an equality constraint
                    tmprow=zeros(1,n)
# inst\qp.m:257
                    tmprow[i]=1
# inst\qp.m:258
                    A=concat([[A],[tmprow]])
# inst\qp.m:259
                    b=concat([[b],[dot(0.5,(lb(i) + ub(i)))]])
# inst\qp.m:260
                    n_eq += 1
                else:
                    tmprow=zeros(1,n)
# inst\qp.m:263
                    tmprow[i]=1
# inst\qp.m:264
                    Ain=concat([[Ain],[tmprow],[- tmprow]])
# inst\qp.m:265
                    bin=concat([[bin],[lb(i)],[- ub(i)]])
# inst\qp.m:266
                    n_in += 2
    
    # Validate inequality constraints.
    if (nargs > 7 and isempty(A_in) and logical_not((isempty(A_lb) or isempty(A_ub)))):
        warning('qp: empty inequality constraint matrix but non-empty bound vectors')
    
    if (nargs > 7 and logical_not(isempty(A_in))):
        dimA_in,n1=size(A_in,nargout=2)
# inst\qp.m:279
        if (n1 != n):
            error('qp: inequality constraint matrix has incorrect column dimension, expected %i',n1)
        else:
            if (logical_not(isempty(A_lb))):
                if (numel(A_lb) != dimA_in):
                    error('qp: inequality constraint matrix and lower bound vector are inconsistent, %i != %i',dimA_in,numel(A_lb))
                else:
                    if (isempty(A_ub)):
                        Ain=concat([[Ain],[A_in]])
# inst\qp.m:287
                        bin=concat([[bin],[A_lb]])
# inst\qp.m:288
            if (logical_not(isempty(A_ub))):
                if (numel(A_ub) != dimA_in):
                    error('qp: inequality constraint matrix and upper bound vector are inconsistent, %i != %i',dimA_in,numel(A_ub))
                else:
                    if (isempty(A_lb)):
                        Ain=concat([[Ain],[- A_in]])
# inst\qp.m:295
                        bin=concat([[bin],[- A_ub]])
# inst\qp.m:296
            if (logical_not(isempty(A_lb)) and logical_not(isempty(A_ub))):
                rtol=copy(tol)
# inst\qp.m:301
                for i in arange(1,dimA_in).reshape(-1):
                    if (abs(A_lb(i) - A_ub(i)) < dot(rtol,(1 + max(abs(A_lb(i) + A_ub(i)))))):
                        # These are actually an equality constraint
                        tmprow=A_in(i,arange())
# inst\qp.m:306
                        A=concat([[A],[tmprow]])
# inst\qp.m:307
                        b=concat([[b],[dot(0.5,(A_lb(i) + A_ub(i)))]])
# inst\qp.m:308
                        n_eq += 1
                    else:
                        tmprow=A_in(i,arange())
# inst\qp.m:311
                        Ain=concat([[Ain],[tmprow],[- tmprow]])
# inst\qp.m:312
                        bin=concat([[bin],[A_lb(i)],[- A_ub(i)]])
# inst\qp.m:313
                        n_in += 2
    
    # Now we should have the following QP:
    
    #   min_x  0.5*x'*H*x + x'*q
  #   s.t.   A*x = b
  #          Ain*x >= bin
    
    # Discard inequality constraints that have -Inf bounds since those
  # will never be active.
    idx=(bin == - Inf)
# inst\qp.m:329
    bin[idx]=[]
# inst\qp.m:331
    Ain[idx,arange()]=[]
# inst\qp.m:332
    n_in=numel(bin)
# inst\qp.m:334
    
    if (isa(x0,'single') or isa(H,'single') or isa(q,'single') or isa(A,'single') or isa(b,'single')):
        rtol=sqrt(eps('single'))
# inst\qp.m:339
    else:
        rtol=copy(tol)
# inst\qp.m:341
    
    eq_infeasible=(n_eq > 0 and norm(dot(A,x0) - b) > dot(rtol,(1 + abs(b))))
# inst\qp.m:344
    in_infeasible=(n_in > 0 and any(dot(Ain,x0) - bin < dot(- rtol,(1 + abs(bin)))))
# inst\qp.m:345
    info=0
# inst\qp.m:347
    if (eq_infeasible or in_infeasible):
        # The initial guess is not feasible.
    # First, define an xbar that is feasible with respect to the
    # equality constraints.
        if (eq_infeasible):
            if (rank(A) < n_eq):
                error('qp: equality constraint matrix must be full row rank')
            xbar=dot(pinv(A),b)
# inst\qp.m:356
        else:
            xbar=copy(x0)
# inst\qp.m:358
        # Second, check that xbar is also feasible with respect to the
    # inequality constraints.
        if (n_in > 0):
            res=dot(Ain,xbar) - bin
# inst\qp.m:364
            if (any(res < dot(- rtol,(1 + abs(bin))))):
                # xbar is not feasible with respect to the inequality constraints.
        # Compute a step in the null space of the equality constraints,
        # by solving a QP.  If the slack is small, we have a feasible initial
        # guess.  Otherwise, the problem is infeasible.
                if (n_eq > 0):
                    Z=null(A)
# inst\qp.m:371
                    if (isempty(Z)):
                        # The problem is infeasible because A is square and full rank,
            # but xbar is not feasible.
                        info=6
# inst\qp.m:375
                if (info != 6):
                    # Solve an LP with additional slack variables
          # to find a feasible starting point.
                    gamma=eye(n_in)
# inst\qp.m:382
                    if (n_eq > 0):
                        Atmp=concat([dot(Ain,Z),gamma])
# inst\qp.m:384
                        btmp=- res
# inst\qp.m:385
                    else:
                        Atmp=concat([Ain,gamma])
# inst\qp.m:387
                        btmp=copy(bin)
# inst\qp.m:388
                    ctmp=concat([[zeros(n - n_eq,1)],[ones(n_in,1)]])
# inst\qp.m:390
                    lb=concat([[- Inf(n - n_eq,1)],[zeros(n_in,1)]])
# inst\qp.m:391
                    ub=[]
# inst\qp.m:392
                    ctype=repmat('L',n_in,1)
# inst\qp.m:393
                    P,dummy,status=glpk(ctmp,Atmp,btmp,lb,ub,ctype,nargout=3)
# inst\qp.m:394
                    if ((status == 0) and all(abs(P(arange(n - n_eq + 1,end()))) < dot(rtol,(1 + norm(btmp))))):
                        # We found a feasible starting point
                        if (n_eq > 0):
                            x0=xbar + dot(Z,P(arange(1,n - n_eq)))
# inst\qp.m:399
                        else:
                            x0=P(arange(1,n))
# inst\qp.m:401
                    else:
                        # The problem is infeasible
                        info=6
# inst\qp.m:405
            else:
                # xbar is feasible.  We use it a starting point.
                x0=copy(xbar)
# inst\qp.m:410
        else:
            # xbar is feasible.  We use it a starting point.
            x0=copy(xbar)
# inst\qp.m:414
    
    if (info == 0):
        # The initial (or computed) guess is feasible.  Call the solver.
        x,lambda_,info,iter=__qp__(x0,H,q,A,b,Ain,bin,maxit,rtol,nargout=4)
# inst\qp.m:420
    else:
        iter=0
# inst\qp.m:422
        x=copy(x0)
# inst\qp.m:423
        lambda_=[]
# inst\qp.m:424
    
    if (isargout(2)):
        obj=dot(dot(dot(0.5,x.T),H),x) + dot(q.T,x)
# inst\qp.m:427
    
    if (isargout(3)):
        INFO.solveiter = copy(iter)
# inst\qp.m:430
        INFO.info = copy(info)
# inst\qp.m:431
    
    return x,obj,INFO,lambda_
    
if __name__ == '__main__':
    pass
    
    # Test infeasible initial guess
    
    
    
    
    
    
    
    
    