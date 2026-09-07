%% EXPERIMENT4
% Experiment 4 from:
%
%   A. H. Al-Mohy,
%   "Computing Linear Combinations of phi-Function Actions
%    for Exponential Integrators"
%
% Advection-Diffusion-Reaction (ADR) test with EXPRK4S6.
%
% The script compares PHIMV, BAMPHI, and KIOPS on
%
%   u_t = epsilon*Delta*u - alpha*(u_x+u_y)
%         + gamma*u*(u-1/2)*(1-u),
%
% on Omega = [0,1]^2, t in [0,1/2], with homogeneous Neumann
% conditions on x=0 and y=0 and homogeneous Dirichlet conditions
% on x=1 and y=1.
%
% Spatial discretization:
%   second-order centered finite differences,
%   Nx = Ny = 100 unknowns, hence n = 10000.
%
% Time integration:
%   EXPRK4S6, implemented with four grouped phi-combination calls
%   per time step.  The PHIMV implementation follows the scheme used
%   in the paper.
%
% Step sizes:
%   h_i = 2^(-8)*(tend-t0)/(i+1), i=1,...,7,
% equivalently
%   Nt = 512, 768, 1024, 1280, 1536, 1792, 2048.
%
% Solver tolerance:
%   tol = 1e-7.
%
% Reference:
%   ODE15S with RelTol = AbsTol = 1e-12.
%   An exact sparse Jacobian is supplied to ODE15S.  Only three requested
%   output states are stored, and only the final state is retained.
%
% The external BAMPHI and KIOPS packages are installed by setup.m.
%
% -------------------------------------------------------------------------

clear
clc
close all
format compact


%% ------------------------------------------------------------------------
% Repository paths
% -------------------------------------------------------------------------

expdir      = fileparts(mfilename('fullpath'));
rootdir     = fileparts(expdir);
externaldir = fullfile(rootdir,'external');
resultsdir  = fullfile(expdir,'results');
figuresdir  = fullfile(expdir,'figures');

addpath(rootdir);

add_external_package(externaldir,'bamphi','bamphi.m');
add_external_package(externaldir,'kiops','kiops.m');

rehash;

if ~exist(resultsdir,'dir')
    mkdir(resultsdir);
end

if ~exist(figuresdir,'dir')
    mkdir(figuresdir);
end


%% ------------------------------------------------------------------------
% Check dependencies
% -------------------------------------------------------------------------

required_functions = {'phimv','find_optimal_param','bamphi','kiops'};

fprintf('Checking dependencies:\n');

for k = 1:numel(required_functions)

    fname = required_functions{k};
    fpath = which(fname);

    if isempty(fpath)
        error(['Required function "%s" was not found.\n' ...
               'Run setup.m from the repository root first.'],fname);
    end

    fprintf('  %-20s %s\n',fname,fpath);
end

% Guard against accidentally using BAMPHI's bundled KIOPS.
kiops_path = which('kiops');

if contains(lower(kiops_path),[filesep 'bamphi' filesep])
    error(['MATLAB is resolving KIOPS through the BAMPHI tree:\n%s\n' ...
           'Run setup.m again and correct the MATLAB path.'],kiops_path);
end

fprintf('\n');


%% ------------------------------------------------------------------------
% Problem parameters
% -------------------------------------------------------------------------

t0   = 0;
tend = 1/2;

epsilon = 1e-3;
alpha   = -1/2;
gamma   = 1000;

Nx = 100;
Ny = 100;
n  = Nx*Ny;

tol = 1e-7;

% h_i = 2^(-8)*(tend-t0)/(i+1), i=1,...,7.
ii   = 1:7;
hvec = 2^(-8)*(tend-t0)./(ii+1);
Nt   = round((tend-t0)./hvec);

fprintf('ADR problem:\n');
fprintf('  Nx = Ny = %d, n = %d\n',Nx,n);
fprintf('  epsilon = %.1e, alpha = %.2f, gamma = %.0f\n', ...
        epsilon,alpha,gamma);
fprintf('  t interval = [%.1f, %.1f]\n',t0,tend);
fprintf('  tolerance = %.1e\n',tol);
fprintf('  time-step counts = %s\n\n',mat2str(Nt));


%% ------------------------------------------------------------------------
% Spatial discretization
% -------------------------------------------------------------------------
%
% For each coordinate we retain the left Neumann boundary node and
% exclude the right Dirichlet boundary node:
%
%     x_j = j/Nx,  j=0,...,Nx-1,
%     y_j = j/Ny,  j=0,...,Ny-1.
%
% At the left Neumann boundary, centered ghost-point elimination gives
%
%     u_x(0)  = 0,
%     u_xx(0) ~ 2(u_1-u_0)/h^2.
%
% At the right Dirichlet boundary, the missing boundary value is zero.

hx = 1/Nx;
hy = 1/Ny;

x = (0:Nx-1)'*hx;
y = (0:Ny-1)'*hy;

% Build the one-dimensional ADR operators first.  This avoids forming and
% storing the four n-by-n sparse matrices Dx, Dy, Dxx, and Dyy separately.
D1x = first_derivative_mixed_bc(Nx,hx);
D2x = second_derivative_mixed_bc(Nx,hx);
Lx  = epsilon*D2x - alpha*D1x;

D1y = first_derivative_mixed_bc(Ny,hy);
D2y = second_derivative_mixed_bc(Ny,hy);
Ly  = epsilon*D2y - alpha*D1y;

% Kronecker-sum representation of
%   epsilon*Delta - alpha*(d/dx + d/dy).
A = kron(speye(Ny),Lx) + kron(Ly,speye(Nx));

Afun = @(X) A*X;

% The initial condition is separable:
%
%   u_0(x,y) = 256 [x(1-x)]^2 [y(1-y)]^2.
%
% Construct it directly in vector form and avoid storing the Nx-by-Ny
% coordinate arrays X and Y.
ux = (x.*(1-x)).^2;
uy = (y.*(1-y)).^2;
u0 = 256*kron(uy,ux);

g  = @(t,u) gamma*u.*(u-1/2).*(1-u); %#ok<INUSD>
dg = @(u) -3*gamma*(u.^2-u+1/6);

fprintf('Sparse ADR operator constructed: nnz(A) = %d\n\n',nnz(A));


%% ------------------------------------------------------------------------
% Reference solution with ODE15S
% -------------------------------------------------------------------------

reference_file = fullfile(resultsdir,'experiment4_reference.mat');
force_reference = false;

if exist(reference_file,'file') && ~force_reference

    tmp = load(reference_file,'uref','Nx','Ny','epsilon','alpha','gamma', ...
               't0','tend');

    reference_matches = ...
        isfield(tmp,'uref')    && ...
        isfield(tmp,'Nx')      && tmp.Nx == Nx && ...
        isfield(tmp,'Ny')      && tmp.Ny == Ny && ...
        isfield(tmp,'epsilon') && tmp.epsilon == epsilon && ...
        isfield(tmp,'alpha')   && tmp.alpha == alpha && ...
        isfield(tmp,'gamma')   && tmp.gamma == gamma && ...
        isfield(tmp,'t0')      && tmp.t0 == t0 && ...
        isfield(tmp,'tend')    && tmp.tend == tend;

else

    reference_matches = false;
end

if reference_matches

    uref = tmp.uref;
    fprintf('Loaded ODE15S reference from:\n  %s\n\n',reference_file);

else

    fprintf('Computing ODE15S reference (RelTol = AbsTol = 1e-12)...\n');

    rhs = @(t,u) A*u + g(t,u); %#ok<NASGU>
    jac = @(t,u) A + spdiags(dg(u),0,n,n); %#ok<NASGU>

    ode_opts = odeset( ...
        'RelTol',1e-12, ...
        'AbsTol',1e-12, ...
        'Jacobian',jac);

    % We need only the final reference state.  Supplying more than two
    % requested output times makes ODE15S return the solution only at
    % those requested times instead of at its complete internal mesh.
    % Hence only three n-vectors are stored here.
    tspan_ref = [t0, (t0+tend)/2, tend];

    tref = tic;

    [~,Uref] = ode15s(rhs,tspan_ref,u0,ode_opts);
    uref = Uref(end,:).';

    reference_time = toc(tref);

    clear Uref

    save(reference_file, ...
        'uref','Nx','Ny','epsilon','alpha','gamma','t0','tend', ...
        'reference_time');

    fprintf('Reference completed in %.2f s.\n',reference_time);
    fprintf('Saved to:\n  %s\n\n',reference_file);
end


%% ------------------------------------------------------------------------
% PHIMV parameters
% -------------------------------------------------------------------------
%
% The shift and scaling parameter depend only on A, so they are computed
% once and reused throughout the time integration, as described in the
% paper.  This one-time setup cost is not included in the time-marching
% timings below.

m = 61;

fprintf('Computing PHIMV shift/scaling parameters once...\n');

tsetup = tic;
[s_phi,shift_phi] = find_optimal_param(Afun,n,m,tol);
setup_time_phi = toc(tsetup);

fprintf('  s = %.16g\n',s_phi);
fprintf('  shift = %.16g\n',shift_phi);
fprintf('  setup time = %.3f s\n\n',setup_time_phi);


%% ------------------------------------------------------------------------
% Optional untimed JIT warm-up
% -------------------------------------------------------------------------
%
% The warm-up is performed on one step only and all method-specific state
% is discarded afterwards.  Thus the measured runs still start with fresh
% BAMPHI spectral information and fresh KIOPS Krylov-size estimates.

do_warmup = true;

if do_warmup

    fprintf('Performing one untimed warm-up step for each method...\n');

    hw = hvec(1);

    exprk4s6_phimv_step(u0,hw,t0,Afun,g,tol,s_phi,shift_phi);

    info_warm = cell(1,4);
    exprk4s6_bamphi_step(u0,hw,t0,Afun,g,tol,info_warm);

    mopt_warm = 10*ones(1,4);
    exprk4s6_kiops_step(u0,hw,t0,Afun,g,tol,mopt_warm);

    fprintf('Warm-up complete.\n\n');
end


%% ------------------------------------------------------------------------
% Main benchmark
% -------------------------------------------------------------------------

method_names = {'phimv','bamphi','kiops'};
num_methods  = numel(method_names);
num_runs     = numel(Nt);

times  = NaN(num_runs,num_methods);
errors = NaN(num_runs,num_methods);

fprintf('Starting Experiment 4.\n');
fprintf('%s\n',repmat('=',1,72));


for irun = 1:num_runs

    h     = hvec(irun);
    nstep = Nt(irun);

    fprintf('\nRun %d/%d: Nt = %d, h = %.16e\n', ...
            irun,num_runs,nstep,h);


    % ================================================================
    % PHIMV
    % ================================================================

    u = u0;

    tt = tic;

    for istep = 1:nstep

        tn = t0 + (istep-1)*h;

        u = exprk4s6_phimv_step( ...
            u,h,tn,Afun,g,tol,s_phi,shift_phi);
    end

    times(irun,1) = toc(tt);
    errors(irun,1) = norm(u-uref,inf)/norm(uref,inf);

    fprintf('  phimv  : %8.3f s, rel_inf = %.6e\n', ...
            times(irun,1),errors(irun,1));


    % ================================================================
    % BAMPHI
    % ================================================================
    %
    % We follow the recommended EXPRK4S6 strategy from the public
    % BAMPHI implementation: one information object is associated with
    % each of the four phi-combination calls and is reused across time
    % steps.  All four stores start empty for each independent run.

    u = u0;
    info_store = cell(1,4);

    tt = tic;

    for istep = 1:nstep

        tn = t0 + (istep-1)*h;

        [u,info_store] = exprk4s6_bamphi_step( ...
            u,h,tn,Afun,g,tol,info_store);
    end

    times(irun,2) = toc(tt);
    errors(irun,2) = norm(u-uref,inf)/norm(uref,inf);

    fprintf('  bamphi : %8.3f s, rel_inf = %.6e\n', ...
            times(irun,2),errors(irun,2));


    % ================================================================
    % KIOPS
    % ================================================================
    %
    % The returned Krylov dimension from each of the four calls is reused
    % as the starting dimension for the corresponding call at the next
    % time step, as recommended in the public EXPRK4S6/KIOPS code.
    %
    % Because grouped EXPRK4S6 requires task1=false, the current standalone
    % KIOPS interface requires m_init, mmin, and mmax to be supplied
    % explicitly.  We use 10, 10, and 128, exactly its built-in defaults.

    u = u0;
    m_opt = 10*ones(1,4);

    tt = tic;

    for istep = 1:nstep

        tn = t0 + (istep-1)*h;

        [u,m_opt] = exprk4s6_kiops_step( ...
            u,h,tn,Afun,g,tol,m_opt);
    end

    times(irun,3) = toc(tt);
    errors(irun,3) = norm(u-uref,inf)/norm(uref,inf);

    fprintf('  kiops  : %8.3f s, rel_inf = %.6e\n', ...
            times(irun,3),errors(irun,3));

end

fprintf('\n%s\n',repmat('=',1,72));


%% ------------------------------------------------------------------------
% Save numerical results
% -------------------------------------------------------------------------

Tout = table( ...
    (1:num_runs).', ...
    Nt(:), ...
    hvec(:), ...
    times(:,1),errors(:,1), ...
    times(:,2),errors(:,2), ...
    times(:,3),errors(:,3), ...
    'VariableNames',{ ...
        'i','Nt','h', ...
        'phimv_time','phimv_error', ...
        'bamphi_time','bamphi_error', ...
        'kiops_time','kiops_error'});

csv_file = fullfile(resultsdir,'experiment4_results.csv');
writetable(Tout,csv_file);

mat_file = fullfile(resultsdir,'experiment4_results.mat');

save(mat_file, ...
    'times','errors','Nt','hvec','uref', ...
    'Nx','Ny','epsilon','alpha','gamma','t0','tend','tol', ...
    's_phi','shift_phi','setup_time_phi','method_names');

fprintf('\nResults saved to:\n');
fprintf('  %s\n',csv_file);
fprintf('  %s\n',mat_file);


%% ------------------------------------------------------------------------
% Time-accuracy figure
% -------------------------------------------------------------------------

fig = figure;

semilogy(times(:,2),errors(:,2),'-o','LineWidth',1.2);
hold on
semilogy(times(:,3),errors(:,3),'-x','LineWidth',1.2);
semilogy(times(:,1),errors(:,1),'-s','LineWidth',1.2);

grid on
box on

xlabel('running time (s)');
ylabel('||\cdot||_\infty-relative error');

legend( ...
    'exprk4s6 with BAMPHI', ...
    'exprk4s6 with KIOPS', ...
    'exprk4s6 with PHIMV', ...
    'Location','best');

figure_png = fullfile(figuresdir,'experiment4_adr.png');
figure_pdf = fullfile(figuresdir,'experiment4_adr.pdf');

exportgraphics(fig,figure_png,'Resolution',300);
exportgraphics(fig,figure_pdf,'ContentType','vector');

fprintf('\nFigure saved to:\n');
fprintf('  %s\n',figure_png);
fprintf('  %s\n',figure_pdf);


%% ------------------------------------------------------------------------
% Command-window summary
% -------------------------------------------------------------------------

fprintf('\n=== EXPERIMENT 4 RESULTS ===\n\n');

disp(Tout);

fprintf('Experiment 4 completed successfully.\n');


%% ========================================================================
% Local functions
% ========================================================================

function D1 = first_derivative_mixed_bc(N,h)
%FIRST_DERIVATIVE_MIXED_BC
% Second-order centered first derivative with
% homogeneous Neumann at the left boundary and
% homogeneous Dirichlet at the right boundary.

    v = ones(N,1)/(2*h);

    D1 = spdiags(v*[-1,1],[-1,1],N,N);

    % Ghost-point elimination for u_x(0)=0.
    D1(1,2) = 0;

end


% -------------------------------------------------------------------------

function D2 = second_derivative_mixed_bc(N,h)
%SECOND_DERIVATIVE_MIXED_BC
% Second-order centered second derivative with
% homogeneous Neumann at the left boundary and
% homogeneous Dirichlet at the right boundary.

    v = ones(N,1)/(h^2);

    D2 = spdiags(v*[1,-2,1],[-1,0,1],N,N);

    % Ghost-point elimination for u_x(0)=0:
    % u_{-1}=u_1.
    D2(1,2) = 2/(h^2);

end


% -------------------------------------------------------------------------

function u = exprk4s6_phimv_step(u,k,t,A,g,tol,s,shift)
%EXPRK4S6_PHIMV_STEP One EXPRK4S6 step using four PHIMV calls.
%
% The scaling and shift parameters are computed once outside the time
% integration and supplied explicitly to every PHIMV call.
%
% In the final call, t=k and alpha=1 are kept separate, exploiting the
% decoupling of the phi-function argument and the polynomial weights.

    gn = g(t,u);
    Fn = A(u) + gn;

    U = zeros(length(u),2);


    % Stage 2.
    U(:,1) = phimv( ...
        (1/2)*k,(1/2)*k,A,[],Fn,tol,s,shift);

    U(:,1) = u + U(:,1);


    % Stages 4 and 3, evaluated together.
    U = phimv( ...
        k*[1/3,1/2], ...
        k*[1/3,1/2], ...
        A,[], ...
        [Fn, ...
         (2/k)*(g(t+(1/2)*k,U(:,1))-gn)], ...
        tol,s,shift);

    % U(:,2) corresponds to stage 3 (c_3=1/2).
    % U(:,1) corresponds to stage 4 (c_4=1/3).
    U(:,2) = g(t+(1/2)*k,u+U(:,2)) - gn;   % D_3
    U(:,1) = g(t+(1/3)*k,u+U(:,1)) - gn;   % D_4


    % Stages 6 and 5, evaluated together.
    U = phimv( ...
        k*[1/3,5/6], ...
        k*[1/3,5/6], ...
        A,[], ...
        [Fn, ...
         (6/k)*(-(2/3)*U(:,2)+(3/2)*U(:,1)), ...
         (12/k^2)*(2*U(:,2)-3*U(:,1))], ...
        tol,s,shift);

    % U(:,2) corresponds to stage 5 (c_5=5/6).
    % U(:,1) corresponds to stage 6 (c_6=1/3).
    U(:,2) = g(t+(5/6)*k,u+U(:,2)) - gn;   % D_5
    U(:,1) = g(t+(1/3)*k,u+U(:,1)) - gn;   % D_6


    % Final update.
    z = phimv( ...
        k,1,A,[], ...
        [Fn, ...
         2*(-(2/5)*U(:,2)+(5/2)*U(:,1)), ...
         4*((6/5)*U(:,2)-3*U(:,1))], ...
        tol,s,shift);

    u = u + k*z;

end

% -------------------------------------------------------------------------

function [u,info_store] = exprk4s6_bamphi_step( ...
    u,h,t,A,g,tol,info_store)
%EXPRK4S6_BAMPHI_STEP One EXPRK4S6 step using four BAMPHI calls.
%
% The four information structures are recycled independently across time
% steps, matching the recommended strategy in BAMPHI's public EXPRK4S6
% implementation.

    gn = g(t,u);
    Fn = A(u) + gn;

    z = zeros(size(u));

    opts = [];
    opts.tol = tol;


    % Stage 2.
    [Z,info_store{1}] = bamphi( ...
        (1/2)*h,A,[],[z,Fn],opts,info_store{1});

    U2 = u + Z;

    D2 = g(t+(1/2)*h,U2) - gn;


    % Stages 4 and 3, in this output order.
    [Z,info_store{2}] = bamphi( ...
        h*[1/3,1/2],A,[], ...
        [z,Fn,(2/h)*D2], ...
        opts,info_store{2});

    U4 = u + Z(:,1);
    U3 = u + Z(:,2);

    D3 = g(t+(1/2)*h,U3) - gn;
    D4 = g(t+(1/3)*h,U4) - gn;


    % Stages 6 and 5, in this output order.
    v2 = (6/h)*((3/2)*D4-(2/3)*D3);
    v3 = (12/h^2)*(2*D3-3*D4);

    [Z,info_store{3}] = bamphi( ...
        h*[1/3,5/6],A,[], ...
        [z,Fn,v2,v3], ...
        opts,info_store{3});

    U6 = u + Z(:,1);
    U5 = u + Z(:,2);

    D5 = g(t+(5/6)*h,U5) - gn;
    D6 = g(t+(1/3)*h,U6) - gn;


    % Final update.
    v2 = (2/h)*((5/2)*D6-(2/5)*D5);
    v3 = (4/h^2)*((6/5)*D5-3*D6);

    [Z,info_store{4}] = bamphi( ...
        h,A,[],[z,Fn,v2,v3], ...
        opts,info_store{4});

    u = u + Z;

end


% -------------------------------------------------------------------------

function [u,m_opt] = exprk4s6_kiops_step(u,h,t,A,g,tol,m_opt)
%EXPRK4S6_KIOPS_STEP One EXPRK4S6 step using four KIOPS calls.
%
% The current standalone KIOPS implementation uses defaults
%
%   m_init = 10,  mmin = 10,  mmax = 128.
%
% Grouped EXPRK4S6 evaluation requires task1=false.  Since task1 is the
% eighth argument, KIOPS requires the Krylov-size arguments to be supplied
% explicitly here; the bounds below are exactly its built-in defaults.
% The returned m_opt from each call is recycled at the next time step.

    mmin = 10;
    mmax = 128;

    gn = g(t,u);
    Fn = A(u) + gn;

    z = zeros(size(u));


    % Stage 2.
    [Z,m_opt(1)] = kiops( ...
        (1/2)*h,A,[z,Fn],tol,m_opt(1),mmin,mmax,false);

    U2 = u + Z;

    D2 = g(t+(1/2)*h,U2) - gn;


    % Stages 4 and 3, in this output order.
    [Z,m_opt(2)] = kiops( ...
        h*[1/3,1/2],A, ...
        [z,Fn,(2/h)*D2], ...
        tol,m_opt(2),mmin,mmax,false);

    U4 = u + Z(:,1);
    U3 = u + Z(:,2);

    D3 = g(t+(1/2)*h,U3) - gn;
    D4 = g(t+(1/3)*h,U4) - gn;


    % Stages 6 and 5, in this output order.
    v2 = (6/h)*((3/2)*D4-(2/3)*D3);
    v3 = (12/h^2)*(2*D3-3*D4);

    [Z,m_opt(3)] = kiops( ...
        h*[1/3,5/6],A, ...
        [z,Fn,v2,v3], ...
        tol,m_opt(3),mmin,mmax,false);

    U6 = u + Z(:,1);
    U5 = u + Z(:,2);

    D5 = g(t+(5/6)*h,U5) - gn;
    D6 = g(t+(1/3)*h,U6) - gn;


    % Final update.
    v2 = (2/h)*((5/2)*D6-(2/5)*D5);
    v3 = (4/h^2)*((6/5)*D5-3*D6);

    [Z,m_opt(4)] = kiops( ...
        h,A,[z,Fn,v2,v3], ...
        tol,m_opt(4),mmin,mmax,false);

    u = u + Z;

end


% -------------------------------------------------------------------------

function add_external_package(externaldir,pkgname,marker_name)
%ADD_EXTERNAL_PACKAGE Add only the directory containing the package's main
%routine.  This avoids recursively adding bundled copies of other solvers.

    pkgroot = fullfile(externaldir,pkgname);

    if ~exist(pkgroot,'dir')
        return
    end

    marker = dir(fullfile(pkgroot,'**',marker_name));

    if isempty(marker)
        return
    end

    depths = zeros(numel(marker),1);

    for k = 1:numel(marker)
        rel = erase(marker(k).folder,pkgroot);
        depths(k) = numel(strfind(rel,filesep));
    end

    [~,idx] = min(depths);

    addpath(marker(idx).folder,'-begin');

end
