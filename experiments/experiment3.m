%% EXPERIMENT3
% Experiment 3 from:
%
%   A. H. Al-Mohy,
%   "Computing Linear Combinations of phi-Function Actions
%    for Exponential Integrators"
%
% This script compares PHIMV, BAMPHI, and KIOPS on three large
% matrix-free low-rank operators
%
%       A = U*V',
%
% where U has orthonormal DCT-II columns and V = U*M'.
% Hence V'*U = M, and the small core M controls the spectrum
% and nonnormality of A.
%
% The script reproduces the numerical data reported in the three
% tables of Experiment 3.
%
% Expected repository structure:
%
%   phimv/
%   |
%   |-- phimv.m
%   |-- find_optimal_param.m
%   |-- setup.m
%   |
%   |-- experiments/
%   |   |-- experiment1.m
%   |   |-- experiment2.m
%   |   |-- experiment3.m
%   |   `-- results/
%   |
%   `-- external/
%       |-- phi_funm/
%       |-- bamphi/
%       `-- kiops/
%
% The external open-source dependencies are installed by setup.m.
% The small-core reference values are computed with phi_funm, as
% prescribed in the paper.
%
% MATLAB version used for the paper: R2022b.
%
% -------------------------------------------------------------------------

clear
clc
close all
format compact

rng(1,'twister');


%% ------------------------------------------------------------------------
%  Repository paths
% -------------------------------------------------------------------------

expdir      = fileparts(mfilename('fullpath'));
rootdir     = fileparts(expdir);
externaldir = fullfile(rootdir,'external');
resultsdir  = fullfile(expdir,'results');

addpath(rootdir);

% Add only the intended external package directories.  Do not use
% addpath(genpath(externaldir)), because BAMPHI contains a bundled copy
% of KIOPS that may otherwise shadow the standalone KIOPS installation.
add_external_package(externaldir,'phi_funm','phi_funm.m');
add_external_package(externaldir,'bamphi','bamphi.m');
add_external_package(externaldir,'kiops','kiops.m');

rehash;

if ~exist(resultsdir,'dir')
    mkdir(resultsdir);
end


%% ------------------------------------------------------------------------
%  Check dependencies
% -------------------------------------------------------------------------

required_functions = { ...
    'phimv', ...
    'phi_funm', ...
    'bamphi', ...
    'kiops'};

fprintf('Checking dependencies:\n');

for k = 1:numel(required_functions)

    fname = required_functions{k};
    fpath = which(fname);

    if isempty(fpath)
        error(['Required function "%s" was not found.\n' ...
               'Run setup.m or check the repository structure.'],fname);
    end

    fprintf('  %-18s %s\n',fname,fpath);
end

% Ensure that MATLAB is using the standalone KIOPS installation rather
% than the copy bundled with BAMPHI.
kiops_path = which('kiops');
if contains(lower(kiops_path), [filesep 'bamphi' filesep])
    error(['MATLAB is resolving KIOPS through the BAMPHI directory:\n%s\n' ...
           'Run setup.m again and check the MATLAB path.'],kiops_path);
end

fprintf('\n');


%% ------------------------------------------------------------------------
%  Timeout configuration
% -------------------------------------------------------------------------

timeout_duration = 10*60;       % 10 minutes

% A hard timeout requires Parallel Computing Toolbox because MATLAB must
% execute the solver on a worker that can be cancelled.  If PCT is not
% available, the script falls back to a serial "soft" timeout: the call
% cannot be interrupted, but it is marked as a timeout if its elapsed time
% exceeds timeout_duration after returning.

have_pct = license('test','Distrib_Computing_Toolbox') && ...
           ~isempty(ver('parallel'));

use_hard_timeout = have_pct;

pool_created_here = false;

if use_hard_timeout

    pool = gcp('nocreate');

    if isempty(pool)
        fprintf('Starting a one-worker parallel pool for hard timeouts...\n');
        parpool('local',1);
        pool_created_here = true;
    elseif pool.NumWorkers ~= 1
        warning(['An existing parallel pool with %d workers was found. ' ...
                 'The experiment submits one task at a time.'], ...
                 pool.NumWorkers);
    end

    fprintf('Hard timeout enabled: %d seconds per solver call.\n\n', ...
            timeout_duration);

else

    warning(['Parallel Computing Toolbox was not found. ' ...
             'The experiment will use serial calls with a soft timeout.']);

    fprintf('Soft timeout: %d seconds per solver call.\n\n', ...
            timeout_duration);
end


%% ------------------------------------------------------------------------
%  Test cases
% -------------------------------------------------------------------------

case_names = {'imgEig','nonNorm','Moler'};
num_cases  = numel(case_names);

method_names = {'phimv','bamphi','kiops'};
num_methods  = numel(method_names);

max_tvals = 5;

times   = NaN(num_cases,max_tvals,num_methods);
errors  = NaN(num_cases,max_tvals,num_methods);
status  = strings(num_cases,max_tvals,num_methods);

tgrid  = cell(num_cases,1);
n_list = zeros(num_cases,1);
p_list = zeros(num_cases,1);

fprintf('Starting Experiment 3.\n');
fprintf('Cases: %s\n',strjoin(case_names,', '));
fprintf('%s\n',repmat('=',1,72));


%% ------------------------------------------------------------------------
%  Main benchmark
% -------------------------------------------------------------------------

for icase = 1:num_cases

    %% Select small core, problem size, and time grid

    switch case_names{icase}

        case 'imgEig'

            Mcore = [ ...
                 0,  1e1; ...
                -1e1, 0];

            tvals = [1e-1,1,10,50,100];

            n = 2e5;
            p = 3;


        case 'nonNorm'

            Mcore = [ ...
                -1, 1e5; ...
                 0, -10];

            tvals = [1e-1,1,10,50,100];

            n = 4e5;
            p = 4;


        case 'Moler'

            a = 2e10;
            b = 4e8/6;
            c = 200/3;
            d = 3;
            e = 1e-8;

            Mcore = [ ...
                 0,      e, 0; ...
                -(a+b), -d, a; ...
                 c,      0, -c];

            tvals = [1e-5,1e-3,1e-1,1,10];

            n = 5e5;
            p = 2;


        otherwise

            error('Unknown test case "%s".',case_names{icase});
    end


    %% Record case metadata

    tgrid{icase}  = tvals;
    n_list(icase) = n;
    p_list(icase) = p;


    %% Random vectors b_0,...,b_p
    %
    % This call occurs at the same point in the random-number sequence as
    % in the original experiment.

    B = randn(n,p+1);


    %% Construct the low-rank matrix-free operator
    %
    % U contains the first r orthonormal DCT-II basis vectors.
    % Vfac = U*Mcore', so that
    %
    %       Vfac'*U = Mcore,
    %
    % and
    %
    %       A = U*Vfac'.

    r = size(Mcore,1);

    [U,Vfac,Afun] = make_dct_low_rank_operator(n,Mcore);
    
    fprintf('\nCase %d/%d: %s\n',icase,num_cases,case_names{icase});
    fprintf('  n = %d, rank = %d, p = %d\n',n,r,p);
    fprintf('  t = %s\n',mat2str(tvals));


    %% Loop over step sizes

    for it = 1:numel(tvals)

        t = tvals(it);

        fprintf('\n  t = %.3e\n',t);


        % ================================================================
        % Reference solution
        % ================================================================
        %
        % Compute
        %
        %   yex = sum_{j=0}^p phi_j(tA)b_j
        %
        % from the small core Mcore using
        %
        %   A^j = U*Mcore^(j-1)*Vfac',  j >= 1.
        %

        try

            yex = sum_phi_low_rank(U,Vfac,Mcore,B,t);

            fprintf('    reference : completed\n');

        catch ME

            error(['Reference computation failed for %s at t = %.3e:\n%s'], ...
                  case_names{icase},t,ME.message);
        end


        % ================================================================
        % PHIMV
        % ================================================================

        [times(icase,it,1), ...
         errors(icase,it,1), ...
         status(icase,it,1)] = run_solver( ...
             @() phimv_error(t,Afun,B,yex), ...
             timeout_duration, ...
             use_hard_timeout, ...
             case_names{icase},t,'phimv');


        % ================================================================
        % BAMPHI
        % ================================================================

        [times(icase,it,2), ...
         errors(icase,it,2), ...
         status(icase,it,2)] = run_solver( ...
             @() bamphi_error(t,Afun,B,p,yex), ...
             timeout_duration, ...
             use_hard_timeout, ...
             case_names{icase},t,'bamphi');


        % ================================================================
        % KIOPS
        % ================================================================

        [times(icase,it,3), ...
         errors(icase,it,3), ...
         status(icase,it,3)] = run_solver( ...
             @() kiops_error(t,Afun,B,p,yex), ...
             timeout_duration, ...
             use_hard_timeout, ...
             case_names{icase},t,'kiops');


        %% Print one-line summary

        fprintf(['    summary   : phimv  %s | ' ...
                 'bamphi %s | kiops %s\n'], ...
                 result_string(times(icase,it,1), ...
                               errors(icase,it,1), ...
                               status(icase,it,1)), ...
                 result_string(times(icase,it,2), ...
                               errors(icase,it,2), ...
                               status(icase,it,2)), ...
                 result_string(times(icase,it,3), ...
                               errors(icase,it,3), ...
                               status(icase,it,3)));

    end

    fprintf('\n%s\n',repmat('-',1,72));

end


%% ------------------------------------------------------------------------
%  Timeout / failure summary
% -------------------------------------------------------------------------

fprintf('\n=== FAILURE / TIMEOUT SUMMARY ===\n');

any_failure = false;

for icase = 1:num_cases

    for it = 1:numel(tgrid{icase})

        for imethod = 1:num_methods

            st = status(icase,it,imethod);

            if st ~= "ok"

                any_failure = true;

                fprintf('%-8s : %-8s, t = %.3e, status = %s\n', ...
                    method_names{imethod}, ...
                    case_names{icase}, ...
                    tgrid{icase}(it), ...
                    st);
            end
        end
    end
end

if ~any_failure
    fprintf('No failures or timeouts.\n');
end


%% ------------------------------------------------------------------------
%  Performance summary
% -------------------------------------------------------------------------

fprintf('\n=== DETAILED PERFORMANCE REPORT ===\n');

for imethod = 1:num_methods

    configured = false(num_cases,max_tvals);

    for icase = 1:num_cases
        configured(icase,1:numel(tgrid{icase})) = true;
    end

    st = status(:,:,imethod);

    ok_mask      = configured & (st == "ok");
    timeout_mask = configured & (st == "timeout");
    bd_mask      = configured & (st == "breakdown");

    total_calls   = nnz(configured);
    success_calls = nnz(ok_mask);
    timeout_calls = nnz(timeout_mask);
    bd_calls      = nnz(bd_mask);

    valid_times = times(:,:,imethod);
    valid_times(~ok_mask) = NaN;

    avg_time = mean(valid_times(:),'omitnan');

    fprintf('\n--- %s ---\n',method_names{imethod});

    fprintf('Successful completions : %d/%d (%.1f%%)\n', ...
        success_calls,total_calls,100*success_calls/total_calls);

    fprintf('Timeouts               : %d/%d (%.1f%%)\n', ...
        timeout_calls,total_calls,100*timeout_calls/total_calls);

    fprintf('Breakdowns              : %d/%d (%.1f%%)\n', ...
        bd_calls,total_calls,100*bd_calls/total_calls);

    fprintf('Average successful time : %.3f s\n',avg_time);
end


%% ------------------------------------------------------------------------
%  Save MAT file with full metadata
% -------------------------------------------------------------------------

results = struct();

results.errors           = errors;
results.times            = times;
results.status           = status;
results.case_names       = case_names;
results.method_names     = method_names;
results.tgrid            = tgrid;
results.n_list           = n_list;
results.p_list           = p_list;
results.timeout_duration = timeout_duration;
results.hard_timeout     = use_hard_timeout;
results.timestamp        = datetime('now');
results.machine_info     = computer;
results.matlab_version   = version;
results.phimv_file        = which('phimv');
results.phi_funm_file     = which('phi_funm');
results.bamphi_file       = which('bamphi');
results.kiops_file        = which('kiops');

mat_file = fullfile(resultsdir,'experiment3_results.mat');

save(mat_file,'-struct','results');

fprintf('\nFull results saved to:\n  %s\n',mat_file);


%% ------------------------------------------------------------------------
%  Save one CSV table per test case
% -------------------------------------------------------------------------

for icase = 1:num_cases

    nt = numel(tgrid{icase});

    Tout = table( ...
        tgrid{icase}(:), ...
        squeeze(times(icase,1:nt,1)).', ...
        squeeze(errors(icase,1:nt,1)).', ...
        squeeze(status(icase,1:nt,1)).', ...
        squeeze(times(icase,1:nt,2)).', ...
        squeeze(errors(icase,1:nt,2)).', ...
        squeeze(status(icase,1:nt,2)).', ...
        squeeze(times(icase,1:nt,3)).', ...
        squeeze(errors(icase,1:nt,3)).', ...
        squeeze(status(icase,1:nt,3)).', ...
        'VariableNames',{ ...
            't', ...
            'phimv_time','phimv_error','phimv_status', ...
            'bamphi_time','bamphi_error','bamphi_status', ...
            'kiops_time','kiops_error','kiops_status'});

    csv_file = fullfile( ...
        resultsdir, ...
        sprintf('experiment3_%s.csv',case_names{icase}));

    writetable(Tout,csv_file);

    fprintf('CSV results saved to:\n  %s\n',csv_file);
end


%% ------------------------------------------------------------------------
%  Print complete LaTeX tables in the manuscript format
% -------------------------------------------------------------------------

fprintf('\n\n=== LATEX-READY TABLES ===\n');

for icase = 1:num_cases

    print_latex_table( ...
        icase, ...
        tgrid{icase}, ...
        squeeze(times(icase,:,:)), ...
        squeeze(errors(icase,:,:)), ...
        squeeze(status(icase,:,:)));

end


%% ------------------------------------------------------------------------
%  Clean up pool created by this script
% -------------------------------------------------------------------------

if pool_created_here

    pool = gcp('nocreate');

    if ~isempty(pool)
        delete(pool);
    end
end

fprintf('\nExperiment 3 completed successfully.\n');


%% ========================================================================
%  Local functions
% ========================================================================

function [U,Vfac,Afun] = make_dct_low_rank_operator(n,M)
%MAKE_DCT_LOW_RANK_OPERATOR Construct A = U*Vfac' from a DCT-II basis.
%
% U contains the first r orthonormal DCT-II basis vectors, where
% r = size(M,1), and
%
%       Vfac = U*M'.
%
% Hence
%
%       Vfac'*U = M.
%
% Afun applies A without forming the n-by-n matrix.

    r = size(M,1);

    i = (0:n-1)';
    k = 0:r-1;

    U = sqrt(2/n) * cos(pi*(i + 0.5) * (k/n));
    U(:,1) = U(:,1)/sqrt(2);

    Vfac = U*M.';

    Afun = @(X) U*(Vfac.'*X);

end


% -------------------------------------------------------------------------

function y = sum_phi_low_rank(U,Vfac,M,B,t)
%SUM_PHI_LOW_RANK Reference for
%
%       y = sum_{i=0}^p phi_i(tA)b_i,
%
% where A = U*Vfac' and Vfac'*U = M.
%
% For j >= 1,
%
%       A^j = U*M^(j-1)*Vfac',
%
% which gives
%
%   sum_i phi_i(tA)b_i
%      =
%   sum_i b_i/i!
%      +
%   U * [ t * sum_i phi_{i+1}(tM)*(Vfac'*b_i) ].

    p = size(B,2)-1;

    Phi = phi_funm(t*M,1:(p+1));

    C = Vfac.'*B;

    invfact = [1,cumprod(1./(1:p))];

    y = B*invfact(:);

    small_sum = zeros(size(C,1),1,'like',C);

    for i = 0:p
        small_sum = small_sum ...
            + t*Phi{i+1}*C(:,i+1);
    end

    y = y + U*small_sum;

end


% -------------------------------------------------------------------------

function e = phimv_error(t,Afun,B,yex)
%PHIMV_ERROR Relative 1-norm error for PHIMV.

    y = phimv(t,1,Afun,B(:,1),B(:,2:end));

    e = norm(yex-y,1)/norm(yex,1);

end


% -------------------------------------------------------------------------

function e = bamphi_error(t,Afun,B,p,yex)
%BAMPHI_ERROR Relative 1-norm error for BAMPHI.
%
% BAMPHI evaluates the conventional weighted combination involving t^j.
% Scaling b_j by t^{-j} makes its target equal to
%
%       sum_j phi_j(tA)b_j.

    Bscaled = B*diag(t.^(0:-1:-p));

    y = bamphi(t,Afun,[],Bscaled);

    e = norm(yex-y,1)/norm(yex,1);

end


% -------------------------------------------------------------------------

function e = kiops_error(t,Afun,B,p,yex)
%KIOPS_ERROR Relative 1-norm error for KIOPS.
%
% KIOPS uses its own default internal Krylov parameters.  Its default
% task1=true applies a final factor t^(-p), so scale the input columns by
% t^(p-j):
%
%   t^(-p) * sum_{j=0}^p t^j*phi_j(tA)*t^(p-j)b_j
%      = sum_{j=0}^p phi_j(tA)b_j.
%
% The convergence tolerance eps/2 is prescribed by this experiment; all
% other optional KIOPS parameters are left at their internal defaults.

    Bscaled = B*diag(t.^(p:-1:0));

    y = kiops(t,Afun,Bscaled,eps/2);

    e = norm(yex-y,1)/norm(yex,1);

end


% -------------------------------------------------------------------------

function [execution_time,error_val,status] = run_solver( ...
    method_func,timeout_duration,use_hard_timeout, ...
    case_name,t_val,method_name)
%RUN_SOLVER Execute one solver call and record timing/status.

    execution_time = NaN;
    error_val       = NaN;
    status          = "breakdown";

    start_time = tic;

    if use_hard_timeout

        f = [];

        try

            % Pass the anonymous function in a cell, matching the mechanism
            % used in the original experiment.
            f = parfeval(@parallel_worker,1,{method_func});

            [completed,result] = fetchNext(f,timeout_duration);

            execution_time = toc(start_time);

            if completed

                if isscalar(result) && isfinite(result)
                    error_val = result;
                    status    = "ok";
                else
                    status = "breakdown";
                end

            else

                cancel(f);

                status = "timeout";

                fprintf(['    %-9s: TIMEOUT for %s, t = %.3e ' ...
                         '(>%d s)\n'], ...
                         method_name,case_name,t_val,timeout_duration);
            end

        catch ME

            execution_time = toc(start_time);

            if ~isempty(f)
                try
                    cancel(f);
                catch
                end
            end

            msg = lower(ME.message);

            if contains(msg,'timeout') || ...
               contains(msg,'time out') || ...
               contains(msg,'cancel')

                status = "timeout";

                fprintf(['    %-9s: TIMEOUT for %s, t = %.3e ' ...
                         '(>%d s)\n'], ...
                         method_name,case_name,t_val,timeout_duration);

            else

                status = "breakdown";

                fprintf(['    %-9s: BREAKDOWN for %s, t = %.3e: %s\n'], ...
                         method_name,case_name,t_val,ME.message);
            end
        end

    else

        try

            result = method_func();

            execution_time = toc(start_time);

            if execution_time > timeout_duration
                status = "timeout";
            elseif isscalar(result) && isfinite(result)
                error_val = result;
                status    = "ok";
            else
                status = "breakdown";
            end

        catch ME

            execution_time = toc(start_time);
            status = "breakdown";

            fprintf(['    %-9s: BREAKDOWN for %s, t = %.3e: %s\n'], ...
                    method_name,case_name,t_val,ME.message);
        end
    end

end


% -------------------------------------------------------------------------

function result = parallel_worker(func_cell)
%PARALLEL_WORKER Execute the supplied zero-input function.

    result = func_cell{1}();

end


% -------------------------------------------------------------------------

function s = result_string(t,e,status)
%RESULT_STRING Compact status string for command-window reporting.

    switch status

        case "ok"
            s = sprintf('%.2fs / %.2e',t,e);

        case "timeout"
            s = 'TO';

        otherwise
            s = sprintf('BD (%.2fs)',t);
    end

end


% -------------------------------------------------------------------------

function print_latex_table(icase,tvals,times,errors,status)
%PRINT_LATEX_TABLE Print one complete table in the manuscript's syntax.

    switch icase
        case 1
            caption = ['$A=UW^{T}$ with \eqref{M1}. Relative forward ' ...
                       'errors and run times (s) over varying $t$.'];
            label = 'tab:imgeig-compact';

        case 2
            caption = ['$A=UW^{T}$ with \eqref{M2}. Relative forward ' ...
                       'errors and run times (s) over varying $t$. ' ...
                       '\textsc{NA} indicates no answer; \textsc{TO} ' ...
                       'indicates a timeout (10 min).'];
            label = 'tab:nonnorm-compact';

        case 3
            caption = ['$A=UW^{T}$ with \eqref{M3}. Relative forward ' ...
                       'errors and run times (s) over varying $t$. ' ...
                       '\textsc{BD}, \textsc{NA}, and \textsc{TO} indicate ' ...
                       'breakdown, no answer, and a timeout (10 min), ' ...
                       'respectively.'];
            label = 'tab:moler-compact';

        otherwise
            error('Unexpected case index.');
    end

    fprintf('\n\\begin{table}[t]\n');
    fprintf('\\centering\n');
    fprintf('\\setlength{\\tabcolsep}{4pt}\n');
    fprintf('\\caption{%s}\n',caption);
    fprintf('\\label{%s}\n',label);
    fprintf('\\begin{tabular}{lcc cc cc}\n');
    fprintf('\\toprule\n');
    fprintf(['& \\multicolumn{2}{c}{\\phimv} ' ...
             '& \\multicolumn{2}{c}{\\bamphi} ' ...
             '& \\multicolumn{2}{c}{\\kiops} \\\\\n']);
    fprintf(['\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}' ...
             '\\cmidrule(lr){6-7}\n']);
    fprintf(['$t$ & time (s) & rel.\\ err & time (s) & rel.\\ err ' ...
             '& time (s) & rel.\\ err \\\\\n']);
    fprintf('\\midrule\n');

    for i = 1:numel(tvals)

        fprintf('$%s$',format_t_for_table(tvals(i),icase));

        for m = 1:size(times,2)

            st = status(i,m);

            if st == "ok"

                fprintf(' & %.2f & $%s$', ...
                    times(i,m),format_scientific(errors(i,m)));

            elseif st == "timeout"

                fprintf(' & \\textsc{TO} & \\textsc{NA}');

            else

                if isfinite(times(i,m))
                    fprintf(' & %.2f & \\textsc{BD}',times(i,m));
                else
                    fprintf(' & \\textsc{BD} & \\textsc{NA}');
                end
            end
        end

        fprintf(' \\\\\n');
    end

    fprintf('\\bottomrule\n');
    fprintf('\\end{tabular}\n');
    fprintf('\\end{table}\n');

end


% -------------------------------------------------------------------------

function s = format_t_for_table(t,icase)
%FORMAT_T_FOR_TABLE Match the t syntax used in the manuscript.
%
% Cases 1 and 2:
%    0.1, 1, 10, 50, 100
%
% Case 3:
%    10^{-5}, 10^{-3}, 10^{-1}, 1, 10

    if icase <= 2
        s = sprintf('%g',t);
    else
        if t == 1
            s = '1';
        elseif t == 10
            s = '10';
        else
            exponent = round(log10(t));
            s = sprintf('10^{%d}',exponent);
        end
    end

end


% -------------------------------------------------------------------------

function s = format_scientific(x)
%FORMAT_SCIENTIFIC Format x as 1.65\mathrm{e}{-16}.

    if x == 0
        s = '0';
        return
    end

    exponent = floor(log10(abs(x)));
    mantissa = x/10^exponent;

    % Protect against rounding a mantissa such as 9.995... to 10.00.
    mantissa_rounded = round(mantissa*100)/100;
    if abs(mantissa_rounded) >= 10
        mantissa = mantissa/10;
        exponent = exponent+1;
    end

    if exponent >= 0
        estr = sprintf('+%d',exponent);
    else
        estr = sprintf('%d',exponent);
    end

    s = sprintf('%.2f\\mathrm{e}{%s}',mantissa,estr);

end


% -------------------------------------------------------------------------

function add_external_package(externaldir,pkgname,marker_name)
%ADD_EXTERNAL_PACKAGE Add only the directory containing the main routine.

    pkgroot = fullfile(externaldir,pkgname);

    if ~exist(pkgroot,'dir')
        return
    end

    marker = dir(fullfile(pkgroot,'**',marker_name));

    if isempty(marker)
        return
    end

    % Choose the shallowest occurrence within this package.
    depths = zeros(numel(marker),1);

    for k = 1:numel(marker)
        rel = erase(marker(k).folder,pkgroot);
        depths(k) = numel(strfind(rel,filesep));
    end

    [~,idx] = min(depths);

    addpath(marker(idx).folder,'-begin');

end
