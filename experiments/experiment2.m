%% EXPERIMENT2
% Experiment 2 from:
%
%   A. H. Al-Mohy,
%   "Computing Linear Combinations of phi-Function Actions
%    for Exponential Integrators"
%
% This script compares
%
%       phimv, bamphi, kiops, and phipm
%
% for the computation
%
%   w(t) = sum_{j=0}^p t^j phi_j(tA) v_j
%
% where A is obtained from a Chebyshev spectral collocation
% discretization of the second derivative with homogeneous
% Dirichlet boundary conditions.
%
% The experiment generates the numerical data reported in Table 1
% of the manuscript.
%
% Expected repository structure:
%
%   phimv/
%   |
%   |-- phimv.m
%   |-- find_optimal_param.m
%   |
%   |-- experiments/
%   |   |-- experiment1.m
%   |   |-- experiment2.m
%   |   |-- common/
%   |   |   `-- phi_func_ex.m
%   |   |
%   |   `-- results/
%   |
%   `-- external/
%       |-- bamphi/
%       |-- kiops/
%       `-- phipm/
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

% Directory containing this script.
expdir = fileparts(mfilename('fullpath'));

% Root directory of the phimv repository.
rootdir = fileparts(expdir);

% Auxiliary routines.
commondir = fullfile(expdir,'common');

% External comparison codes.
externaldir = fullfile(rootdir,'external');

% Output directory.
resultsdir = fullfile(expdir,'results');

% Add repository paths.
addpath(rootdir);
addpath(commondir);

if exist(externaldir,'dir')
    addpath(genpath(externaldir));
end

if ~exist(resultsdir,'dir')
    mkdir(resultsdir);
end


%% ------------------------------------------------------------------------
%  Check dependencies
% -------------------------------------------------------------------------

required_functions = { ...
    'phimv', ...
    'phi_func_ex', ...
    'bamphi', ...
    'kiops', ...
    'phipm'};

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

fprintf('\n');


%% ------------------------------------------------------------------------
%  Restrict MATLAB numerical computation to one computational thread
% -------------------------------------------------------------------------

% This was used in the original experiment to improve fairness of
% run-time comparisons.
%
% maxNumCompThreads is supported in the MATLAB version used for the paper,
% although newer MATLAB releases may issue a warning.

try
    maxNumCompThreads(1);
catch
    warning(['Could not restrict MATLAB to one computational thread. ' ...
             'Timing results may therefore differ from those in the paper.']);
end


%% ------------------------------------------------------------------------
%  Chebyshev spectral discretization
% -------------------------------------------------------------------------

N = 100;
L = 2;

j = (0:N)';

x = (cos(pi*j/N)+1)*(L/2);

c = [2;ones(N-1,1);2].*(-1).^j;

X  = repmat(x,1,N+1);
dX = X-X';

D = (c*(1./c)')./(dX+eye(N+1));
D = D-diag(sum(D,2));

% Second derivative on [0,L].
A = (2/L)^2*(D^2);

% Remove boundary rows and columns to enforce homogeneous
% Dirichlet boundary conditions.
A = A(2:N,2:N);

n = size(A,1);


%% ------------------------------------------------------------------------
%  Experiment parameters
% -------------------------------------------------------------------------

% Highest phi-function index.
p = 5;

% Random vectors:
%
%   u = v_0,
%   V = [v_1,...,v_p].
%
u = randn(n,1);
V = randn(n,p);

% Matrix-free operator used by bamphi.
Afun = @(Y) A*Y;

% Step sizes.
tvals = 10.^(-4:0);

% Tolerance explicitly supplied to kiops and phipm in the
% original experiment.
tol = 1e-14;

% Soft timeout.
%
% NOTE: this does not interrupt a running MATLAB routine. It only marks a
% computation as a timeout if the routine eventually returns after TMAX.
TMAX = 5*60;


%% ------------------------------------------------------------------------
%  Precompute high-accuracy reference solutions
% -------------------------------------------------------------------------

fprintf('Computing reference solutions...\n');

References = cell(numel(tvals),1);

for it = 1:numel(tvals)

    tt = tvals(it);

    % High-accuracy phi-function matrices.
    Phi_exact = phi_func_ex(tt*A,0:p);

    % Form
    %
    %   R = sum_{j=0}^p t^j phi_j(tA) v_j.
    %
    R = Phi_exact{1}*u;

    for jpow = 1:p
        R = R ...
            + tt^jpow*Phi_exact{jpow+1}*V(:,jpow);
    end

    References{it} = R;

    fprintf('  t = %.1e completed\n',tt);
end

fprintf('\n');


%% ------------------------------------------------------------------------
%  Solver definitions
% -------------------------------------------------------------------------

solver_names = { ...
    'phimv', ...
    'bamphi', ...
    'kiops', ...
    'phipm'};

ns = numel(solver_names);

Ttime = NaN(numel(tvals),ns);
Terr  = NaN(numel(tvals),ns);
Tstat = strings(numel(tvals),ns);


%% ------------------------------------------------------------------------
%  Benchmark
% -------------------------------------------------------------------------

fprintf('Running benchmark...\n\n');

for it = 1:numel(tvals)

    tt = tvals(it);

    R  = References{it};
    nR = norm(R,1);


    %% phimv

    fprintf('t = %.1e | phimv  : ',tt);

    t0 = tic;

    try

        % This reproduces the call used in the original experiment.
        % Hence phimv uses its default tolerance.
        w = phimv(tt,tt,A,u,V);

        elapsed = toc(t0);

        Ttime(it,1) = elapsed;

        if elapsed > TMAX

            Tstat(it,1) = "timeout";
            fprintf('timeout (%.3f s)\n',elapsed);

        else

            Terr(it,1)  = norm(R-w,1)/nR;
            Tstat(it,1) = "ok";

            fprintf('time = %.3f s, error = %.3e\n', ...
                elapsed,Terr(it,1));
        end

    catch ME

        elapsed = toc(t0);

        Ttime(it,1) = elapsed;
        Tstat(it,1) = "error";

        fprintf('ERROR: %s\n',ME.message);

    end


    %% bamphi

    fprintf('t = %.1e | bamphi : ',tt);

    t0 = tic;

    try

        % Single-output call to bamphi.
        w = bamphi_first_output( ...
            tt,Afun,[],[u,V]);

        elapsed = toc(t0);

        Ttime(it,2) = elapsed;

        if elapsed > TMAX

            Tstat(it,2) = "timeout";
            fprintf('timeout (%.3f s)\n',elapsed);

        else

            Terr(it,2)  = norm(R-w,1)/nR;
            Tstat(it,2) = "ok";

            fprintf('time = %.3f s, error = %.3e\n', ...
                elapsed,Terr(it,2));
        end

    catch ME

        elapsed = toc(t0);

        Ttime(it,2) = elapsed;
        Tstat(it,2) = "error";

        fprintf('ERROR: %s\n',ME.message);

    end


    %% kiops

    fprintf('t = %.1e | kiops   : ',tt);

    t0 = tic;

    try

        w = kiops( ...
            tt,A,[u,V],tol,[],[],[],false);

        elapsed = toc(t0);

        Ttime(it,3) = elapsed;

        if elapsed > TMAX

            Tstat(it,3) = "timeout";
            fprintf('timeout (%.3f s)\n',elapsed);

        else

            Terr(it,3)  = norm(R-w,1)/nR;
            Tstat(it,3) = "ok";

            fprintf('time = %.3f s, error = %.3e\n', ...
                elapsed,Terr(it,3));
        end

    catch ME

        elapsed = toc(t0);

        Ttime(it,3) = elapsed;
        Tstat(it,3) = "error";

        fprintf('ERROR: %s\n',ME.message);

    end


    %% phipm

    fprintf('t = %.1e | phipm   : ',tt);

    t0 = tic;

    try

        w = phipm(tt,A,[u,V],tol);

        elapsed = toc(t0);

        Ttime(it,4) = elapsed;

        if elapsed > TMAX

            Tstat(it,4) = "timeout";
            fprintf('timeout (%.3f s)\n',elapsed);

        else

            Terr(it,4)  = norm(R-w,1)/nR;
            Tstat(it,4) = "ok";

            fprintf('time = %.3f s, error = %.3e\n', ...
                elapsed,Terr(it,4));
        end

    catch ME

        elapsed = toc(t0);

        Ttime(it,4) = elapsed;
        Tstat(it,4) = "error";

        fprintf('ERROR: %s\n',ME.message);

    end

    fprintf('\n');

end


%% ------------------------------------------------------------------------
%  Display summary
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('===============================================================\n');
fprintf('Experiment 2 summary\n');
fprintf('Soft per-call timeout: %d s\n',TMAX);
fprintf('===============================================================\n\n');

fprintf(['%-9s  ' ...
         '%-22s %-22s %-22s %-22s\n'], ...
         't',solver_names{:});

for it = 1:numel(tvals)

    fprintf('%-9.1e ',tvals(it));

    for s = 1:ns

        if Tstat(it,s) == "ok"

            fprintf(' %7.3fs / %8.1e     ', ...
                Ttime(it,s),Terr(it,s));

        else

            fprintf(' %-20s   ',char(Tstat(it,s)));

        end

    end

    fprintf('\n');

end


%% ------------------------------------------------------------------------
%  Save results
% -------------------------------------------------------------------------

results_table = table( ...
    tvals(:), ...
    Ttime(:,1),Terr(:,1),Tstat(:,1), ...
    Ttime(:,2),Terr(:,2),Tstat(:,2), ...
    Ttime(:,3),Terr(:,3),Tstat(:,3), ...
    Ttime(:,4),Terr(:,4),Tstat(:,4), ...
    'VariableNames',{ ...
        't', ...
        'phimv_time','phimv_error','phimv_status', ...
        'bamphi_time','bamphi_error','bamphi_status', ...
        'kiops_time','kiops_error','kiops_status', ...
        'phipm_time','phipm_error','phipm_status'});

results_file = fullfile( ...
    resultsdir,'experiment2_results.csv');

writetable(results_table,results_file);

fprintf('\nResults saved to:\n  %s\n',results_file);


%% ------------------------------------------------------------------------
%  Finish
% -------------------------------------------------------------------------

fprintf('\nExperiment 2 completed successfully.\n');


%% ========================================================================
%  Local helper
% ========================================================================

function f = bamphi_first_output(t,A,At,uV)
%BAMPHI_FIRST_OUTPUT Return only the first output of bamphi.
%
% The call below reproduces the default bamphi configuration used
% in the original experiment.

    f = bamphi(t,A,At,uV);

end