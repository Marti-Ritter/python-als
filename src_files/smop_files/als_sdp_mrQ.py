# Generated with SMOP  0.41-beta
from libsmop import *


# als\inst\als_sdp_mrQ.m


@function
def als_sdp_mrQ(data=None, N=None, model=None, estimator=None, varargin=None, *args, **kwargs):
    varargin = als_sdp_mrQ.varargin
    nargin = als_sdp_mrQ.nargin

    ## A modified ALS-SDP function.  Uses Single Column ALS form and imposes 
    ## semidefinite constraints.
    ## Options allow semidefinite constraints to be removed.
    ##
    ## [Qest_cell,Rest_cell,trQ,Phi,Phi_tot,phi0,bhat_mtr,bhat_mtr_thry,cov_bound, ...
    ## Iter,Iter_maxed,timespent] = als_sdp_mrQ(data,N,model,estimator,varargin)
    ##
    ##
    ## Function Inputs:
    ##
    ## data.yk (measurements),
    ## data.uk (inputs),
    ## data.xhatk (state estimates- optional),
    ## data.datapts (number of data points considered, optional),
    ## data.start (data to be ignored in the beginning until initial condition is
    ## negligible - optional, default is 100),
    ## N (window size),
    ## model.A,
    ## model.B (optional),
    ## model.C,
    ## model.G (optional, default is identity matrix)
    ## estimator.Q (Qw initial guess, optional)
    ## estimator.R (Rv initial guess, optional),
    ## estimator.L (initial estimator gain - optional)
    ##
    ## Either (Q,R) or L must be supplied.  If L is supplied, (Q,R) are used as
    ## initial guess but not in calculating the estimator gain.  If L is not supplied
    ## then (Q,R) are used to calculate initial estimator
    ##
    ## Optional inputs for vargin:
    ## 'rho_values', [rho1 rho2 rho3 ...] (vector of penalties on trace(Q), default is
    ## [0])
    ## 'Rform', 'sym' or 'Rform', 'diag' (imposes symmetric or diagonal constraints on
    ## Rv, default is diagonal)
    ## 'Weight', 'I', or "Weight', 'data' (imposes identity matrix weighting or
    ## data-based weighting on least-squares term, default is data-based)
    ## 'Plot', 0 (turns off all plots) or 'plot', 1 (produces plots of L-innovations
    ## and plots of autocovariance estimates and fits), default is to produce plots
    ## 'sdp', 0 (solves unconstrained problem), default is 'sdp', 1 (enforces
    ## semidefinite constraints)
    ## 'tracestates', [1 1 ... 0 0] (vector indicating which states are to be
    ## penalized in the trace(Q) term, default is [1 1 ... 1], to penalize all states)
    ##
    ##
    ## Function outputs :
    ##
    ## Qest_cell - cell containing estimated Qw for each penalty on tr(Q)
    ## Rest_cell - cell containing estimated Rv for each penalty on tr(Q)
    ## trQ - vector containing trace(Q) for each penalty on tr(Q)
    ## Phi - vector containing least-squares portion of objective for each
    ## penalty on tr(Q); scaled so that phi = 1 when rho = 0
    ## Phi_tot - vector containing full objective function value for each
    ## penalty on tr(Q)
    ## phi0 - vector containing objective function value at rho = 0 (prior
    ## to scaling)
    ## bhat_mtr - matrix of autocovariances estimated from data, columns are
    ## <y_1(k),y_1(k-j)>, ..., <y_1(k),y_p(k-j)>, ..., <y_p(k),y_1(k-j)>, ...,
    ## <y_p(k),y_p(k-j)>
    ## bhat_mtr - matrix of theoretical autocovariances calculated from model
    ## and ALS results (for smallest penalty on trace(Q))
    ## cov_bound - estimated standard deviations for each auto-cross covariances,
    ## [sigma_11; sigma_12; ... sigma_pp], use +/- 2*sigma as 95# confidence intervals
    ## Iter - vector containing number of iterations for each penalty on tr(Q)
    ## Iter_maxed(i) = 1 if maximum number of iterations are reached for that
    ## penalty on trQ
    ## Timespent - vector of time spent in seconds of CPU time for each
    ## penalty on tr(Q)

    if logical_not(isfield(data, 'datapts')):
        data.datapts = copy(size(data.yk, 2))
    # als\inst\als_sdp_mrQ.m:71

    datapts = data.datapts
    # als\inst\als_sdp_mrQ.m:72
    Aa = model.A
    # als\inst\als_sdp_mrQ.m:74
    Ca = model.C
    # als\inst\als_sdp_mrQ.m:75
    if isfield(model, 'G'):
        Ga = model.G
    # als\inst\als_sdp_mrQ.m:77
    else:
        Ga = eye(size(Aa))
    # als\inst\als_sdp_mrQ.m:78

    n, g = size(Ga, nargout=2)
    # als\inst\als_sdp_mrQ.m:80
    p = size(Ca, 1)
    # als\inst\als_sdp_mrQ.m:81
    if isfield(model, 'B'):
        Ba = model.B
        # als\inst\als_sdp_mrQ.m:83
        m = size(Ba, 2)
    # als\inst\als_sdp_mrQ.m:83
    else:
        Ba = zeros(n)
        # als\inst\als_sdp_mrQ.m:84
        m = copy(n)
        # als\inst\als_sdp_mrQ.m:84
        data.uk = copy(zeros(m, datapts))
    # als\inst\als_sdp_mrQ.m:84

    if isfield(data, 'start'):
        start = data.start
    # als\inst\als_sdp_mrQ.m:87
    else:
        start = 100
    # als\inst\als_sdp_mrQ.m:88

    na = copy(n)
    # als\inst\als_sdp_mrQ.m:92
    ga = copy(g)
    # als\inst\als_sdp_mrQ.m:92
    pa = copy(p)
    # als\inst\als_sdp_mrQ.m:92
    ## Deal with variable arguments

    # Default values:
    lam_vec = 0
    # als\inst\als_sdp_mrQ.m:97

    trstates = ones(ga, 1)
    # als\inst\als_sdp_mrQ.m:98

    Rsym = 0
    # als\inst\als_sdp_mrQ.m:99

    dataweight = 1
    # als\inst\als_sdp_mrQ.m:100

    plot_flag = 1
    # als\inst\als_sdp_mrQ.m:101

    sdp = 1
    # als\inst\als_sdp_mrQ.m:102

    # Inputs
    okargs = cellarray(['rho_values', 'tracestates', 'rform', 'weight', 'plot', 'sdp'])
    # als\inst\als_sdp_mrQ.m:105
    for j in arange(1, (nargin - 4), 2).reshape(-1):
        pname = varargin[j]
        # als\inst\als_sdp_mrQ.m:107
        pval = varargin[j + 1]
        # als\inst\als_sdp_mrQ.m:108
        k = strmatch(lower(pname), okargs)
        # als\inst\als_sdp_mrQ.m:109
        if isempty(k):
            error('als_sdp_mrQ: Unknown parameter name: %s.', pname)
        else:
            if length(k) > 1:
                error('als_sdp_mrQ: Ambiguous parameter name: %s.', pname)
            else:
                if 1 == (k):
                    lam_vec = ravel(pval)
                # als\inst\als_sdp_mrQ.m:117
                else:
                    if 2 == (k):
                        trstates = ravel(pval)
                    # als\inst\als_sdp_mrQ.m:119
                    else:
                        if 3 == (k):
                            if strmatch(lower(pval), 'sym'):
                                Rsym = 1
                            # als\inst\als_sdp_mrQ.m:122
                            else:
                                if strmatch(lower(pval), 'diag'):
                                    Rsym = 0
                                # als\inst\als_sdp_mrQ.m:124
                                else:
                                    warning('als_sdp_mrQ: Unknown structure type for R; defaulting to diagonal')
                        else:
                            if 4 == (k):
                                if strmatch(upper(pval), 'I'):
                                    dataweight = 0
                                # als\inst\als_sdp_mrQ.m:130
                                else:
                                    if strmatch(lower(pval), 'data'):
                                        dataweight = 1
                                    # als\inst\als_sdp_mrQ.m:132
                                    else:
                                        warning('als_sdp_mrQ: Unknown weighting type; defaulting to data-based')
                            else:
                                if 5 == (k):
                                    plot_flag = copy(pval)
                                # als\inst\als_sdp_mrQ.m:137
                                else:
                                    if 6 == (k):
                                        sdp = copy(pval)
    # als\inst\als_sdp_mrQ.m:139

    ## Estimator Simulation
    y = data.yk
    # als\inst\als_sdp_mrQ.m:146
    u = data.uk
    # als\inst\als_sdp_mrQ.m:147
    # Calculate estimator gain if not supplied:
    if isfield(estimator, 'L'):
        L = estimator.L
    # als\inst\als_sdp_mrQ.m:150
    else:
        L = dlqe(Aa, Ga, Ca, estimator.Q, estimator.R)
    # als\inst\als_sdp_mrQ.m:151

    # Estimate states
    if (logical_not(isfield(data, 'xhatk'))):
        xhat = zeros(na, datapts)
        # als\inst\als_sdp_mrQ.m:156
        xhat_ = zeros(n, datapts)
        # als\inst\als_sdp_mrQ.m:157
        xhat_[arange(1, n), 1] = model.xhat0
        # als\inst\als_sdp_mrQ.m:158
        for i in arange(1, datapts).reshape(-1):
            xhat[arange(), i] = xhat_(arange(), i) + dot(L, (y(arange(), i) - dot(Ca, xhat_(arange(), i))))
            # als\inst\als_sdp_mrQ.m:160
            xhat_[arange(), i + 1] = dot(Aa, xhat(arange(), i)) + dot(Ba, u(arange(), i))
        # als\inst\als_sdp_mrQ.m:161
        xhat_ = xhat_(arange(), arange(1, end() - 1))
    # als\inst\als_sdp_mrQ.m:163
    else:
        xhat_ = data.xhatk
    # als\inst\als_sdp_mrQ.m:165

    # Estimate L-innovations
    inntrun = y(arange(), arange(start + 1, end())) - dot(Ca, xhat_(arange(), arange(start + 1, end())))
    # als\inst\als_sdp_mrQ.m:169
    # Optional - use to look at plots to check if error looks white:
    if plot_flag == 1:
        figure
        for i in arange(1, p, 1).reshape(-1):
            if p > 1:
                subplot(ceil(p / 2), 2, i)
            plot(inntrun(i, arange()), 'm', 'linewidth', 2)
            xlabel('Time')
            ylabel(sprintf('Innovation of y_%i', i))

    ## Calculation of Autocorrelations for one column ALS
    datatrun = datapts - start
    # als\inst\als_sdp_mrQ.m:185
    Eyy = []
    # als\inst\als_sdp_mrQ.m:186
    for i in arange(0, N - 1).reshape(-1):
        temp = dot(inntrun(arange(), arange(i + 1, end())), inntrun(arange(), arange(1, end() - i)).T)
        # als\inst\als_sdp_mrQ.m:188
        temp = temp / (datatrun - i)
        # als\inst\als_sdp_mrQ.m:189
        Eyy = concat([[Eyy], [temp]])
    # als\inst\als_sdp_mrQ.m:190

    Eyy = ravel(Eyy)
    # als\inst\als_sdp_mrQ.m:192
    if dataweight == 1:
        ## Calculation of Weighting matrix for one column ALS
        Nd = size(inntrun, 2)
        # als\inst\als_sdp_mrQ.m:196
        nt = Nd - N + 1
        # als\inst\als_sdp_mrQ.m:197
        trials = dot(2, N)
        # als\inst\als_sdp_mrQ.m:198
        datapts2 = floor(nt / trials)
        # als\inst\als_sdp_mrQ.m:199
        Eyy_store2 = []
        # als\inst\als_sdp_mrQ.m:200
        Yysmall = []
        # als\inst\als_sdp_mrQ.m:201
        covb = zeros(dot(pa ** 2, N), dot(pa ** 2, N))
        # als\inst\als_sdp_mrQ.m:202
        for i_trial in arange(1, trials, 1).reshape(-1):
            Yysmall = zeros(dot(N, p), datapts2)
            # als\inst\als_sdp_mrQ.m:204
            for i in arange(1, datapts2, 1).reshape(-1):
                yyst = inntrun(arange(),
                               arange(dot((i - 1), trials) + i_trial, dot((i - 1), trials) + (i_trial - 1) + N))
                # als\inst\als_sdp_mrQ.m:206
                Yysmall[arange(), i] = ravel(yyst)
            # als\inst\als_sdp_mrQ.m:207
            Py = cov(Yysmall.T)
            # als\inst\als_sdp_mrQ.m:209
            Px = Py(arange(1, p), arange(1, p))
            # als\inst\als_sdp_mrQ.m:210
            Pyx = Py(arange(), arange(1, p))
            # als\inst\als_sdp_mrQ.m:211
            covb = covb + dot(1 / datapts2, (kron(Px, Py) + dot(comm_mat(pa, dot(N, pa)), kron(Pyx, Pyx.T))))
        # als\inst\als_sdp_mrQ.m:212
        covb = covb / trials
        # als\inst\als_sdp_mrQ.m:214
        Wm = pinv(covb)
    # als\inst\als_sdp_mrQ.m:215
    else:
        ## Use identity-based weighting
        Wm = eye(length(Eyy))
    # als\inst\als_sdp_mrQ.m:218

    ##################################################
    ## Building the constant matrix for the LS problem
    ##################################################

    Ain = Aa - dot(dot(Aa, L), Ca)
    # als\inst\als_sdp_mrQ.m:225
    OO = []
    # als\inst\als_sdp_mrQ.m:227
    temp = eye(n)
    # als\inst\als_sdp_mrQ.m:228
    for i in arange(1, N).reshape(-1):
        OO = concat([[OO], [dot(Ca, temp)]])
        # als\inst\als_sdp_mrQ.m:230
        temp = dot(temp, Ain)
    # als\inst\als_sdp_mrQ.m:231

    ## temporary variables
    M1 = zeros(na ** 2, ga ** 2)
    # als\inst\als_sdp_mrQ.m:235
    i = 1
    # als\inst\als_sdp_mrQ.m:236
    for j in arange(1, ga).reshape(-1):
        for k in arange(1, ga).reshape(-1):
            II = zeros(ga)
            # als\inst\als_sdp_mrQ.m:239
            II[k, j] = 1
            # als\inst\als_sdp_mrQ.m:240
            t1 = dlyap(Ain, dot(dot(Ga, II), Ga.T))
            # als\inst\als_sdp_mrQ.m:241
            M1[arange(), i] = ravel(t1)
            # als\inst\als_sdp_mrQ.m:242
            i = i + 1
    # als\inst\als_sdp_mrQ.m:243

    M2 = zeros(na ** 2, pa ** 2)
    # als\inst\als_sdp_mrQ.m:246
    i = 1
    # als\inst\als_sdp_mrQ.m:247
    for j in arange(1, pa).reshape(-1):
        for k in arange(1, pa).reshape(-1):
            II = zeros(pa)
            # als\inst\als_sdp_mrQ.m:250
            II[k, j] = 1
            # als\inst\als_sdp_mrQ.m:251
            t2 = dlyap(Ain, dot(dot(dot(dot(Aa, L), II), L.T), Aa.T))
            # als\inst\als_sdp_mrQ.m:252
            M2[arange(), i] = ravel(t2)
            # als\inst\als_sdp_mrQ.m:253
            i = i + 1
    # als\inst\als_sdp_mrQ.m:254

    ##########################################
    ## Single column ALS method
    ##########################################

    PSI = eye(pa)
    # als\inst\als_sdp_mrQ.m:262
    for i in arange(1, N - 1).reshape(-1):
        PSI = concat([[PSI], [dot(dot(dot(- Ca, Ain ** (i - 1)), Aa), L)]])
    # als\inst\als_sdp_mrQ.m:264

    OOtemp = kron(Ca, OO)
    # als\inst\als_sdp_mrQ.m:267
    PSItemp = kron(eye(pa), PSI)
    # als\inst\als_sdp_mrQ.m:268
    LH1 = dot(OOtemp, M1)
    # als\inst\als_sdp_mrQ.m:270
    LH2 = dot(OOtemp, M2) + PSItemp
    # als\inst\als_sdp_mrQ.m:271
    if Rsym == 1:
        LHSsingc = concat([dot(LH1, symtran(ga)), dot(LH2, symtran(pa))])
    # als\inst\als_sdp_mrQ.m:275
    else:
        # Adds symmetric constraint to Q_w and diagonal constraint to R_v
        LHSsingc = concat([dot(LH1, symtran(ga)), LH2(arange(), arange(1, pa ** 2, pa + 1))])
    # als\inst\als_sdp_mrQ.m:278

    M = dot(dot(kron(Ca, eye(n)), M1), symtran(ga))
    # als\inst\als_sdp_mrQ.m:281
    Mrank = rank(M, 0.0001)
    # als\inst\als_sdp_mrQ.m:282
    nr, nc = size(M, nargout=2)
    # als\inst\als_sdp_mrQ.m:283
    if nc > Mrank:
        fprintf('The covariance estimates are not unique!\\n Use trade-off curve to find lowest rank solution\\n')

    if isfield(estimator, 'Q') and isfield(estimator, 'R'):
        Q0 = celldiag(cellarray([estimator.Q, estimator.R]))
    # als\inst\als_sdp_mrQ.m:288
    else:
        Q0 = eye(ga + pa)
    # als\inst\als_sdp_mrQ.m:290

    numax = 1
    # als\inst\als_sdp_mrQ.m:292
    if sdp == 0:
        # Solve least-squares problem without semi-definite constraints
        if lam_vec != 0:
            warning('als_sdp_mrQ: Cannot solve indefinite problem with penalty on trQ')
        warning('off', 'Octave:singular-matrix-div')
        # warning('off','MATLAB:nearlySingularMatrix') ##Matlab compatible
        time1 = copy(cputime)
        # als\inst\als_sdp_mrQ.m:301
        X = dot(dot(numpy.linalg.solve((dot(dot(LHSsingc.T, Wm), LHSsingc)), LHSsingc.T), Wm), Eyy)
        # als\inst\als_sdp_mrQ.m:302
        timespent = cputime - time1
        # als\inst\als_sdp_mrQ.m:303
        Qest_cell[1] = cellarray([reshape(dot(symtran(ga), X(arange(1, dot(ga, (ga + 1)) / 2))), ga, ga)])
        # als\inst\als_sdp_mrQ.m:304
        if Rsym == 1:
            Rest_cell[1] = cellarray([reshape(dot(symtran(pa), X(arange(dot(ga, (ga + 1)) / 2 + 1, end()))), pa, pa)])
        # als\inst\als_sdp_mrQ.m:306
        else:
            Rest_cell[1] = cellarray([diag(X(arange(end() - pa + 1, end())))])
        # als\inst\als_sdp_mrQ.m:308
        Phi = dot(dot((Eyy - dot(LHSsingc, X)).T, Wm), (Eyy - dot(LHSsingc, X)))
        # als\inst\als_sdp_mrQ.m:310
        phi0 = copy(Phi)
        # als\inst\als_sdp_mrQ.m:311
        Phi_tot = copy(Phi)
        # als\inst\als_sdp_mrQ.m:312
        Iter = 0
        # als\inst\als_sdp_mrQ.m:313
        Iter_maxed = 0
        # als\inst\als_sdp_mrQ.m:314
        trQ = trace(dot(Qest_cell[1], diag(trstates)))
    # als\inst\als_sdp_mrQ.m:315
    else:
        # Solve ALS with semi-definite constraints
        # Scale weighting matrix and initial objective function
        QR0, phi0 = sdp_QR_mrQ(LHSsingc, Eyy, Wm, Q0, 0, numax, pa, Rsym, trstates, nargout=2)
        # als\inst\als_sdp_mrQ.m:319
        Wm = Wm / phi0
        # als\inst\als_sdp_mrQ.m:320
        ## Loop through and solve for each value of lam_vec
        for i in arange(1, length(lam_vec), 1).reshape(-1):
            lam = lam_vec(i)
            # als\inst\als_sdp_mrQ.m:323
            time1 = copy(cputime)
            # als\inst\als_sdp_mrQ.m:324
            QR1, phi, phi_tot, iter, nsq, iter_maxed = sdp_QR_mrQ(LHSsingc, Eyy, Wm, Q0, lam, numax, pa, Rsym, trstates,
                                                                  nargout=6)
            # als\inst\als_sdp_mrQ.m:325
            timespent[i, 1] = cputime - time1
            # als\inst\als_sdp_mrQ.m:326
            Qest = QR1(arange(1, ga), arange(1, ga))
            # als\inst\als_sdp_mrQ.m:327
            Rest = QR1(arange(ga + 1, end()), arange(ga + 1, end()))
            # als\inst\als_sdp_mrQ.m:328
            Qest_cell[i] = cellarray([(Qest + Qest.T) / 2])
            # als\inst\als_sdp_mrQ.m:329
            Rest_cell[i] = cellarray([(Rest + Rest.T) / 2])
            # als\inst\als_sdp_mrQ.m:330
            Phi[i, 1] = phi
            # als\inst\als_sdp_mrQ.m:331
            Phi_tot[i, 1] = phi_tot
            # als\inst\als_sdp_mrQ.m:332
            Iter[i, 1] = iter
            # als\inst\als_sdp_mrQ.m:333
            Iter_maxed[i, 1] = iter_maxed
            # als\inst\als_sdp_mrQ.m:334
            trQ[i, 1] = trace(dot(Qest_cell[i], diag(trstates)))
    # als\inst\als_sdp_mrQ.m:335

    ## Form matrices of autocovariances from data and theory
    if Rsym == 1:
        # Eyy_thry = LHSsingc*[symtranT(ga)*vec(Qest_cell{1}); symtranT(pa)*vec(Rest_cell{1})];
        Eyy_thry = dot(LHSsingc,
                       concat([[dot(symtranT(ga), ravel(Qest_cell[1]))], [dot(symtranT(pa), vec(Rest_cell[1]))]]))
    # als\inst\als_sdp_mrQ.m:343
    else:
        # Eyy_thry = LHSsingc*[symtranT(ga)*vec(Qest_cell{1}); diag(Rest_cell{1})];
        Eyy_thry = dot(LHSsingc, concat([[dot(symtranT(ga), ravel(Qest_cell[1]))], [diag(Rest_cell[1])]]))
    # als\inst\als_sdp_mrQ.m:346

    Eyy = reshape(Eyy, dot(N, pa), pa)
    # als\inst\als_sdp_mrQ.m:348
    Eyy_thry = reshape(Eyy_thry, dot(N, pa), pa)
    # als\inst\als_sdp_mrQ.m:349
    for i in arange(1, pa, 1).reshape(-1):
        for j in arange(1, pa, 1).reshape(-1):
            for k in arange(1, N, 1).reshape(-1):
                bhat_mtr[k, dot(pa, (i - 1)) + j] = Eyy(dot(pa, (k - 1)) + i, j)
                # als\inst\als_sdp_mrQ.m:353
                bhat_mtr_thry[k, dot(pa, (i - 1)) + j] = Eyy_thry(dot(pa, (k - 1)) + i, j)
    # als\inst\als_sdp_mrQ.m:354

    npts = datapts - start - N
    # als\inst\als_sdp_mrQ.m:359
    for i in arange(1, pa, 1).reshape(-1):
        for j in arange(1, pa, 1).reshape(-1):
            cov_bound[i, j] = sqrt(dot(bhat_mtr(1, dot(pa, (i - 1)) + i), bhat_mtr(1, dot(pa, (j - 1)) + j)) / (npts))
    # als\inst\als_sdp_mrQ.m:361

    cov_bound = ravel(cov_bound)
    # als\inst\als_sdp_mrQ.m:363
    if plot_flag == 1:
        figure
        for i in arange(1, pa ** 2, 1).reshape(-1):
            subplot(pa, pa, i)
            plot(arange(0, N - 1, 1), bhat_mtr(arange(), i), '-*', arange(0, N - 1, 1), bhat_mtr_thry(arange(), i),
                 '-o', 'linewidth', 2)
            # xlabel('Lag')
            hold('on')
            plot(concat([[0], [N]]), dot(2, kron(concat([[1, - 1], [1, - 1]]), concat([cov_bound(i)]))), 'k')
            if i == 1:
                legend('Data', 'Fit')
            if i <= pa:
                title('Autocorrelation')
            if i > pa ** 2 - pa:
                xlabel('Lag')
            hold('off')
