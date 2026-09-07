%% EXPERIMENT1
% First numerical experiment from:
%
% A. H. Al-Mohy,
% "Computing Linear Combinations of phi-Function Actions
%  for Exponential Integrators"
%
% This script compares PHIMV and PHI_FUNM and generates:
%
%   Figure 1: Relative forward errors
%   Figure 2: Accuracy performance profile
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
%   |   |-- common/
%   |   |   |-- testmats.m
%   |   |   |-- phi_func_ex.m
%   |   |   |-- expm_mp.m
%   |   |   |-- expm_frechet.m
%   |   |   |-- expm_frechet_cond.m
%   |   |   |-- fd_condest.m
%   |   |   |-- expm_frechet_cn.m
%   |   |   `-- myperfprof.m
%   |   |
%   |   |-- results/
%   |   `-- figures/
%   |
%   `-- external/
%       |-- phi_funm/
%       |   `-- phi_funm.m
%       |
%       `-- anymatrix/
%           |-- anymatrix.m
%           `-- ...
%
% MATLAB version used for the paper: R2022b.
%
% -------------------------------------------------------------------------

clear
clc
close all
format compact


%% ------------------------------------------------------------------------
% Repository paths
% -------------------------------------------------------------------------

% Directory containing this script.
expdir = fileparts(mfilename('fullpath'));

% Root directory of the phimv repository.
rootdir = fileparts(expdir);

% Auxiliary routines used by the experiments.
commondir = fullfile(expdir,'common');

% External comparison/dependency codes.
externaldir = fullfile(rootdir,'external');

% Output directories.
resultsdir = fullfile(expdir,'results');
figdir     = fullfile(expdir,'figures');

% Add repository paths.
addpath(rootdir);
addpath(commondir);

if exist(externaldir,'dir')
    addpath(genpath(externaldir));
end

% Create output directories if needed.
if ~exist(resultsdir,'dir')
    mkdir(resultsdir);
end

if ~exist(figdir,'dir')
    mkdir(figdir);
end


%% ------------------------------------------------------------------------
% Check dependencies
% -------------------------------------------------------------------------

required_functions = { ...
    'phimv', ...
    'find_optimal_param',...
    'testmats', ...
    'anymatrix', ...
    'phi_func_ex', ...
    'phi_funm', ...
    'fd_condest', ...
    'expm_frechet_cn', ...
    'myperfprof',...
    'expm_frechet',...
    'expm_frechet_cond',...
    'expm_mp'};

fprintf('Checking dependencies:\n');

for j = 1:numel(required_functions)

    fname = required_functions{j};
    fpath = which(fname);

    if isempty(fpath)
        error(['Required function "%s" was not found.\n' ...
               'Please run setup.m or check the repository structure.'], ...
               fname);
    end

    fprintf('  %-20s %s\n',fname,fpath);
end

fprintf('\n');


%% ------------------------------------------------------------------------
% Experiment parameters
% -------------------------------------------------------------------------

ids_min = 1;

[~,ids_max] = testmats();

num_mats  = ids_max - ids_min + 1;
n_default = 20;

% Highest phi-function index.
pdeg = 5;

% Linear combination:
%
%   w = sum_{j=0}^p alpha^j phi_j(tA)v_j.
%
t     = 1;
alpha = 1;

% Unit roundoff / prescribed tolerance.
tol = 2^(-53);

% Generate figure files.
generate_figures = true;


%% ------------------------------------------------------------------------
% Storage
% -------------------------------------------------------------------------

err_phimv   = NaN(num_mats,1);
err_phifunm = NaN(num_mats,1);
kappa       = NaN(num_mats,1);


%% ------------------------------------------------------------------------
% Main experiment
% -------------------------------------------------------------------------

fprintf('Running Experiment 1 on %d test matrices.\n\n',num_mats);

for k = ids_min:ids_max

    idx = k - ids_min + 1;

    % Reproducible random vectors for each test matrix.
    rng(1234+k,'twister');


    %% Generate test matrix

    A = testmats(k,n_default);

    % Scalings used in the original experiment.
    if k == 14 || k == 46
        A = A/10;
    elseif k == 80
        A = A/20;
    end

    n = size(A,1);


    %% Generate v_0,...,v_p

    V = randn(n,pdeg+1);

    % This reproduces the experiment used in the paper:    
    % V(:,1) = zeros(n,1);


    %% High-accuracy reference solution

    phi_exact = phi_func_ex(t*A,0:pdeg);

    w_exact = zeros(n,1);

    for j = 1:pdeg+1
        w_exact = w_exact ...
            + alpha^(j-1)*phi_exact{j}*V(:,j);
    end

    norm_exact = norm(w_exact,2);


    %% phi_funm

    phi_approx = phi_funm(t*A,0:pdeg);

    w_phifunm = zeros(n,1);

    for j = 1:pdeg+1
        w_phifunm = w_phifunm ...
            + alpha^(j-1)*phi_approx{j}*V(:,j);
    end


    %% phimv

    w_phimv = phimv( ...
        t, ...
        alpha, ...
        A, ...
        V(:,1), ...
        V(:,2:end), ...
        tol);


    %% Relative forward errors

    err_phimv(idx) = ...
        norm(w_exact-w_phimv,2)/norm_exact;

    err_phifunm(idx) = ...
        norm(w_exact-w_phifunm,2)/norm_exact;


    %% Condition-number estimate

    Vr = V(:,2:end);

    J = diag(ones(pdeg-1,1),1);

    Eblock = [ ...
        zeros(size(A)),  Vr; ...
        zeros(size(Vr')), zeros(pdeg)];

    kappa(idx) = fd_condest( ...
        @expm_frechet_cn, ...
        t*blkdiag(A,J), ...
        Eblock);


    %% Progress information

    fprintf(['Matrix %3d: kappa = %.3e | ' ...
             'err_phimv = %.3e | err_phi_funm = %.3e\n'], ...
             k, ...
             kappa(idx), ...
             err_phimv(idx), ...
             err_phifunm(idx));

end


%% ------------------------------------------------------------------------
% Save numerical results
% -------------------------------------------------------------------------

matrix_id = (ids_min:ids_max).';

Tresults = table( ...
    matrix_id, ...
    err_phimv, ...
    err_phifunm, ...
    kappa, ...
    kappa.*tol, ...
    'VariableNames', ...
    {'matrix_id', ...
     'err_phimv', ...
     'err_phi_funm', ...
     'kappa', ...
     'kappa_tol'});

results_file = fullfile( ...
    resultsdir,'experiment1_results.csv');

writetable(Tresults,results_file);

fprintf('\nNumerical results saved to:\n');
fprintf('  %s\n',results_file);


%% ------------------------------------------------------------------------
% Sort by decreasing condition number
% -------------------------------------------------------------------------

% IMPORTANT:
% Sorting uses the actual computed condition numbers.
[~,order] = sort(kappa,'descend');

kappa_sorted       = kappa(order);
err_phimv_sorted   = err_phimv(order);
err_phifunm_sorted = err_phifunm(order);

rank_idx = 1:num_mats;

% Actual reference levels.
kappa_tol = kappa_sorted*tol;

% Figure 1 has y <= 1.
% Clip ONLY the plotted reference curve to the displayed range.
% The true condition numbers and kappa*tol values remain unchanged
% in the saved numerical results.
kappa_tol_plot = min(kappa_tol,1);


%% ------------------------------------------------------------------------
% Plotting parameters
% -------------------------------------------------------------------------

lg_linewidth  = 1.2;
lg_markersize = 4;
lg_fontsize   = 12;

ax_linewidth = 1.0;
ax_fontsize  = 10;

color_cond    = [0,0,0];
color_phifunm = [0.23,0.48,0.34];
color_phimv   = [0.635,0.078,0.184];


%% ------------------------------------------------------------------------
% Figure 1: relative forward errors
% -------------------------------------------------------------------------

figure(1)
clf

% Publication-friendly figure dimensions.
set(gcf,'Position',[100,100,720,520]);

ga = gobjects(3,1);

ga(1) = semilogy( ...
    rank_idx, ...
    kappa_tol_plot, ...
    '-', ...
    'Color',color_cond, ...
    'LineWidth',lg_linewidth, ...
    'DisplayName','$\kappa\cdot\mathrm{tol}$');

hold on

ga(2) = semilogy( ...
    rank_idx, ...
    err_phifunm_sorted, ...
    'v', ...
    'Color',color_phifunm, ...
    'LineWidth',lg_linewidth, ...
    'MarkerSize',lg_markersize, ...
    'DisplayName','\texttt{phi\_funm}');

ga(3) = semilogy( ...
    rank_idx, ...
    err_phimv_sorted, ...
    'o', ...
    'Color',color_phimv, ...
    'LineWidth',lg_linewidth, ...
    'MarkerSize',lg_markersize, ...
    'DisplayName','\texttt{phimv}');

grid on
box on

% Horizontal axis.
xlim([1,num_mats]);

xt = unique([1,20:20:100,num_mats]);
xt = xt(xt >= 1 & xt <= num_mats);
xticks(xt);

% Logarithmic vertical axis.
ylim([1e-18,1]);

yticks(10.^(-18:3:0));

yticklabels({ ...
    '$10^{-18}$', ...
    '$10^{-15}$', ...
    '$10^{-12}$', ...
    '$10^{-9}$', ...
    '$10^{-6}$', ...
    '$10^{-3}$', ...
    '$10^{0}$'});

set(gca, ...
    'LineWidth',ax_linewidth, ...
    'FontSize',ax_fontsize, ...
    'TickLabelInterpreter','latex');

legend( ...
    ga, ...
    'NumColumns',1, ...
    'FontSize',lg_fontsize, ...
    'Interpreter','latex', ...
    'Location','northeast');

hold off


%% Save Figure 1

if generate_figures

    eps_file = fullfile(figdir,'fig_relerr.eps');
    png_file = fullfile(figdir,'fig_relerr.png');

    print(gcf,eps_file,'-depsc2','-painters');
    print(gcf,png_file,'-dpng','-r300');

    fprintf('\nFigure 1 saved to:\n');
    fprintf('  %s\n',eps_file);
    fprintf('  %s\n',png_file);
end


%% ------------------------------------------------------------------------
% Figure 2: accuracy performance profile
% -------------------------------------------------------------------------

figure(2)
clf

set(gcf,'Position',[100,100,720,520]);

% Retain the layout used in the original figure.
subplot(3,1,[1,2]);

ratio_max = 12;

error_matrix = [ ...
    err_phifunm_sorted, ...
    err_phimv_sorted];

method_names = { ...
    '\texttt{phi\_funm}', ...
    '\texttt{phimv}'};

markers = {'-v','-^'};

mycolors = [ ...
    color_phifunm; ...
    color_phimv];

myperfprof( ...
    error_matrix, ...
    ratio_max, ...
    mycolors, ...
    markers, ...
    lg_markersize, ...
    lg_linewidth);

grid on
box on

ylim([0.2,1]);

set(gca, ...
    'LineWidth',ax_linewidth, ...
    'FontSize',ax_fontsize);

legend( ...
    method_names, ...
    'Interpreter','latex', ...
    'FontSize',lg_fontsize, ...
    'Location','southeast');


%% Save Figure 2

if generate_figures

    eps_file = fullfile(figdir,'fig_perf.eps');
    png_file = fullfile(figdir,'fig_perf.png');

    print(gcf,eps_file,'-depsc2','-painters');
    print(gcf,png_file,'-dpng','-r300');

    fprintf('\nFigure 2 saved to:\n');
    fprintf('  %s\n',eps_file);
    fprintf('  %s\n',png_file);
end


%% ------------------------------------------------------------------------
% Finish
% -------------------------------------------------------------------------

fprintf('\nExperiment 1 completed successfully.\n');