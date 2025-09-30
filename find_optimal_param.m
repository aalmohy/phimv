function [s, shift, mv] = find_optimal_param(A, sizeA, m, tol)
% FIND_OPTIMAL_PARAM  Pick scaling S and spectral shift SHIFT for PHIMV.
%
%   [S, SHIFT, MV] = FIND_OPTIMAL_PARAM(A, SIZEA, M, TOL)
%   returns a scaling factor S and a real shift SHIFT that are used by
%   PHIMVes "scaling-and-recovering method with Taylor series" routine to
%   efficiently and stably evaluate linear combinations of phi-functions
%   exp/phi_j(tA)v_j. The routine also returns MV, the number of A*x
%   products spent during the parameter search (MV=0 for cache hits).
%
% --- Inputs
%   A       : Matrix or function handle. If a matrix (n×n), it is wrapped as
%             a matvec handle A(X)=A*X. If a handle, it must accept both
%             n×1 and n×k inputs and return A*X of the same size.
%
%   SIZEA   : Problem size n (i.e., the number of rows/cols of A). Needed
%             when A is a function handle so the routine can allocate and
%             probe correctly.
%
%   M       : Intended Taylor polynomial degree used later by PHIMV.
%             This routine uses M both to decide how far to probe powers of
%             A and to normalize the objective it minimizes for SHIFT.
%
%   TOL     : Target accuracy for the subsequent PHIMV call. Influences the
%             recommended scaling S via the asymptotic formula
%                 S ≈ s0 * fval / (TOL * M!)^(1/M).
%
% --- Outputs
%   S       : Positive scaling parameter (integer number of recovery steps
%             in PHIMV is ceil(max(|t|)*S)). Larger S → more, smaller
%             time-steps inside PHIMV, stabilizing the Taylor remainder.
%
%   SHIFT   : Real spectral shift used inside PHIMV to build (A − SHIFT*I).
%             The shift is chosen to minimize an upper bound on the tail of
%             the truncated Taylor/φ expansion for the probed direction.
%
%   MV      : Count of matrix-vector applications used here (roughly M).
%             If the same (A, M, TOL) triplet is queried again, the result
%             is served from an internal cache and MV=0.
%
% --- What the routine does (high level)
%   1) If A is a matrix, it is wrapped as a matvec handle for uniformity.
%   2) Draw a random unit vector v and build the Krylov column stack
%      [v, Av, A^2 v, …, A^r v] up to M, stopping early if under/overflow
%      is anticipated (log-safe guards).
%   3) Estimate a “natural scale” s0 from the growth rate of ||A^k v||.
%   4) Rescale the stack by s0 so subsequent operations are well-scaled.
%   5) Minimize (via fminbnd on a small interval around 0) the 1D objective
%         f(ξ) = || ∑_{k=0}^{M} C(M,k) (−ξ/s0)^{M−k} A^k v ||^(1/M),
%      i.e., the (normalized) size of a weighted linear combination of the
%      powers, which models the Taylor/φ tail sensitivity to a real shift.
%      The minimizer defines SHIFT, and the minimized value adjusts S.
%
% --- Caching
%   The function keeps a small MRU cache (default 16 entries) keyed by the
%   tuple {A, M, TOL}. On a cache hit, it returns the previously computed
%   (S, SHIFT) and sets MV=0. This is helpful when the same configuration is
%   used across many experiments. (If A is a handle, the identity of the
%   handle, not the underlying operator, forms part of the key.)
%
% --- When to use
%   Call FIND_OPTIMAL_PARAM once per (A, M, TOL) combination, then pass S
%   and SHIFT into PHIMV. If you sweep over TOL or M, reuse the cache to
%   avoid repeated probing when values match.
%
% --- Example: matrix A
%   n = 200; A = gallery('poisson',n);  % SPD 2D Laplacian
%   [S, SHIFT, MV] = find_optimal_param(A, size(A,1), 60, 1e-8);
%   % -> pass S, SHIFT to PHIMV
%
% --- Example: function-handle A
%   n = 1000; L = delsq(numgrid('S', 34));    % sparse
%   Aop = @(X) L * X;                         % matvec
%   [S, SHIFT, MV] = find_optimal_param(Aop, 2^10, 80, 1e-10);
%
% See also: PHIMV, FMINBND. 
%   Reference:
%
%   Awad H. Al-Mohy, Computing Linear Combinations of φ-function Actions
%            for Exponential Integrators, 2025.

% ---------- Caching ----------
persistent cache
if isempty(cache)
    cache = struct('key', {}, 's', {}, 'shift', {}, 'mv', {});
end

% ---- build cache key BEFORE matvec unification ----
key = {A, m, tol};

% ---- cache lookup ----
for i = 1:numel(cache)
    if isequal(cache(i).key, key)
        % Cache hit: return cached values, set mv=0
        s     = cache(i).s;
        shift = cache(i).shift;
        mv    = 0;
       % fprintf('Cache hit: Using cached results for identical inputs (mv set to 0)\n');
        return
    end
end
% Cache miss: proceed with normal computation
% fprintf('Cache miss: Computing new results\n');

% ---------- matvec unification ----------
if ~isa(A, 'function_handle')
    Amat = A;                 % preserve user input
    A    = @(X) Amat * X;     % matvec
end
n = sizeA;

% ---------- parameters ----------
nrm        = 2;               % vector norm
clamp      = realmin;         % clamp
logRmax    = log(realmax);
logRmin    = log(realmin);
frac       = 0.8;             % headroom
logSafeMax = frac*logRmax;    % realmax^frac guard
logSafeMin = frac*logRmin;
mv         = m;               % matrix-vector product

% ---------- start vector ----------
rng(0,'twister');
v  = randn(n,1);
v  = v / norm(v, nrm);

% ---------- Binomial coefficients ----------
kvec        = 0:m;
logBinom    = gammaln(m+1) - gammaln(kvec+1) - gammaln(m-kvec+1);
logBinomMax = max(logBinom);

% ---------- grow r adaptively (over/undrflow-safe) -------
V        = zeros(n, m+1);
V(:,1)   = v;
logNorm  = -inf(m+1, 1);    % log ||A^k v||; start with ||v||=1 -> log(1)=0
logNorm(1) = 0;

for k = 1:m
    V(:,k+1)   = A(V(:,k));   % A^k v
    logNorm(k+1) = log(norm(V(:,k+1), nrm));
    lognw = logNorm(k+1);
    if lognw + log1p(k)/2 + logBinomMax > logSafeMax || lognw < logSafeMin
        break
    end
end

r = k;  % Now V has r+1 columns: [v, Av, A^2v, ..., A^r v]

% ---------- s0 from last few safe ratios (geometric mean) ----------
j_hi  = r;                    % corresponds to A^{r-1} v. A^r excluded
j_lo  = max(2, j_hi-5);       % use up to 5 ratios
dlogs = diff(logNorm(j_lo:j_hi));
s0    = exp(mean(dlogs));

% ---------- rescale accepted columns by s0^(0:-1:-r) ----------
V(:,1:r+1) = bsxfun(@times, V(:,1:r+1), s0.^(0:-1:-r));

% --------- extend to m with 1/s0 scaling & early gate on s_k ---------
for k = r+1:m
    V(:,k+1) = A(V(:,k)) / s0;
end

% ---------- binomial coefficients (fast & stable) ----------
bicoef = exp(logBinom).';

% ---------- shift search interval ----------
range = sqrt(n) * s0;

% ---------- fminbnd over objective ----------
opts = optimset('Display','off', 'FunValCheck','off', 'TolX', 1e-3);
try
    opts = optimset(opts, 'FunctionTolerance', 1e-3, 'StepTolerance', 1e-3);
catch
end

[shift, fval, exitflag] = fminbnd(@objective, -range, range, opts);
if exitflag ~= 1
    error('convergence is not guaranteed');
end
fval = max(fval, clamp);  % clamp > 0

% ---------- final s ----------
% s = s0 * fval / (tol * m!)^(1/m)
log_s = log(s0) + log(fval) - ( log(tol) + gammaln(m+1) ) / m;
s     = exp(log_s);
 
% ================ STORE RESULTS IN CACHE (with MRU cap) ==================
new_entry.key   = key;
new_entry.s     = s;
new_entry.shift = shift;
new_entry.mv    = mv; % store original mv value for reference
cache(end+1)    = new_entry;
% Cap cache growth: keep most recent MAXCACHE entries
maxCache = 16;
if numel(cache) > maxCache
    cache = cache(end - maxCache + 1:end);
end

%fprintf('Cache updated: Stored new results (cache size: %d)\n', numel(cache));

% ---------- nested objective function ----------
    function f = objective(xi)
        z      = -xi / s0;
        p      = cumprod([1; repmat(z, m, 1)]);  % [1; z; ...; z^m]
        xi_pwr = p(end:-1:1);                    % [z^m; ...; 1]
      % xi_pwr = z.^(m:-1:0).';
        alpha  = xi_pwr .* bicoef;               % weights
        y      = V * alpha;                      % linear combination
        f      = norm(y, nrm)^(1/m);
        f      = max(f, clamp);                  % clamp to avoid zeros
    end
end
