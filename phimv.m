function [W, st, mv] = phimv(t, alpha, A, v0, V, tol, s, shift)
% PHIMV  Action of multiple φ-functions via scaling & recovery (truncated Taylor).
%
%   W           = PHIMV(t, alpha, A, v0, V)
%   [W, st, mv] = PHIMV(t, alpha, A, v0, V, tol)
%   [W, st, mv] = PHIMV(t, alpha, A, v0, V, tol, s, shift)
%
% Purpose
% -------
% Efficiently evaluate *multiple linear combinations* of φ-functions
%   w_i = sum_{j=0}^p alpha_i^j * φ_j(t_i A) v_j,    i = 1..r
% without forming φ_j(tA) explicitly. The routine advances with a scaling
% parameter and a truncated Taylor “remainder” expansion (a.k.a. scaling
% & recovery), and couples the φ-ladders through a small right-side block
% recurrence. It returns the n×r result block W whose columns are w_i.
%
% Inputs
% ------
% t      : length-r vector of times {t_i}.
% alpha  : scalar or length-r vector. Internally forms Δ = diag(alpha(:)).
% A      : n×n matrix *or* function handle applying A(X)=A*X for X (n×k).
% v0     : n×1 vector (may be []) — starting vector for φ_0 (i.e., exp).
% V      : n×p matrix whose columns are [v1,...,vp] (may be []).
%          Here p is the highest φ-index used (φ_0 ≡ exp, φ_1, ... , φ_p).
% tol    : (optional) convergence tolerance for the truncated Taylor series.
%          Default: eps/2.
% s,shift: (optional) scaling S and spectral shift SHIFT. If omitted, they
%          are estimated by FIND_OPTIMAL_PARAM(A, n, m, tol).
%
% Outputs
% 
% Let r = numel(t), p = size(V,2) (0 if V = []).
% Internally the algorithm builds F = [ first n×r block | last n×r block ] and
% returns  W = F(:,1:r) + F(:,r+1:2*r) * Δ,  with Δ = diag(alpha(:)).
%
% Case A: V is unavailable/empty  (V = [], p = 0)
%   • Only the “exp” is requested.
%   • Output:
%         W(:,i) = exp(t(i) * A) * v0,    i = 1..r.
%     (No φ_j terms appear because no v_j were provided.)
%
% Case B: v0 is unavailable/empty  (isempty(v0) == true)
%   • Only the “last φ-block” is requested—from V = [v1, ..., vp].
%   • Output:
%         W(:,i) = sum_{j=1}^p alpha(i)^(j) * φ_j(t(i) * A) * v_j,   i = 1..r.
%     (No exp(tA)v0 term appears because v0 was not provided.)
%
% Case C: both v0 and V are provided  (general case)
%   • Full linear combination is returned:
%         W(:,i) = exp(t(i) * A) * v0  +  sum_{j=1}^p alpha(i)^(j) * φ_j(t(i) * A) * v_j,
%       for i = 1..r.
%
% Shapes and requirements
% -----------------------
%   • W is always n×r (one column per time t(i)).
%   • If alpha is scalar, Δ = alpha * I_r; if vector, Δ = diag(alpha(:)) and
%     must have length r.
%   • If both v0 and V are empty, the call is ill-posed (nothing to compute).
%     In that case, error out.
%
% See also: FIND_OPTIMAL_PARAM, EXPM, PHI_FUNM functions. 
%
%   Reference:
%
%   Awad H. Al-Mohy, Computing Linear Combinations of φ-function Actions
%            for Exponential Integrators, 2025.
%=======================================================================

    if nargin < 5, V = []; end
    if nargin < 6 || isempty(tol), tol = eps/2; end

    n = max(size(v0,1), size(V,1));  
    mv = 0;
    if nargin < 7 || isempty(s) || isempty(shift)
        m = 61;
        [s, shift, mv] = find_optimal_param(A, n , m, tol);
    end    

    % --- Build unshifted and shifted apply-ops without changing user A ---
    isHandle = isa(A,'function_handle');
    if isHandle
        % keep handle as-is (unshifted) and build shifted on-the-fly
        A1_apply = A;                             % unshifted AX
        A_apply  = @(X) A(X) - shift*X;           % shifted (A - shift I)X
    else
        % keep matrix; build shifted matrix once
        In     = speye(n);
        A_shift = A - shift*In;                   % shifted once
        A1_apply = @(X) A * X;                    % unshifted AX
        A_apply  = @(X) A_shift * X;              % shifted (A - shift I)X
    end
    
    emptyv0 = isempty(v0);       % if true, exp(t_iA)v0 is not wanted
    emptyV  = isempty(V);        % if true, exp(t_iA)v0 is the only wanted
    r     = numel(t);
    T     = t(:).'; 
    diagT = diag(T);    
    Alpha = eye(r).*alpha(:);
    tmax  = max(abs(T));
    st    = ceil(tmax * s);
    nrm   = 1;  % chosen norm 
    if tmax == 0, st = 1; end

    mu = exp(T .* shift / st);
    p  = size(V,2);

    J = tril(triu(ones(p), 1), 1); % robust choice to all p >= 0
    if emptyv0
        J1 = J;    p0 = 0;
    else
        J1 = [zeros(1, p+1); [zeros(p, 1), J]];  p0 = 1;
    end

    % Kronecker structures 
    J  = kron(J1, Alpha) - shift * kron(eye(p+p0), diagT);   
    V  = [v0, V(:, end:-1:1)] * kron(eye(p+p0), mu/st);

    % ---------------- Taylor series initialization ----------------
    J    = J / st;    
    krT  = kron(eye(p+p0), diagT); 
    S = V; D = V; 
    k = 1; 
    factk = 1;
    c1 = inf;    
    c2 = norm(D, nrm);   
    sizeD2 = size(D,2);
    while c1 + c2 > tol * norm(S, nrm)
        c1 = c2; k = k + 1; factk = factk * k;
        V  = V * J;
        D  = A_apply(D * (krT / st)) + V;
        mv = mv + sizeD2;
        c2 = norm(D, nrm) / factk;
        S  = S + D / factk;
    end    
    expJs = expm(kron(J1/st,Alpha));     
    % --------------- Extract requested columns -------------------
    idx_last = (p-1)*r + 1 : p*r;
    idx_both = [1:r, p*r + 1 : (p+1)*r];

    if emptyV
        F = S;
    elseif emptyv0
        F = S(:, idx_last);
    else
        F = S(:, idx_both);
    end
    %---------------------------------------------------------------
    if ~emptyv0
        % F(:, 1:r) ≈ exp(tA/s)v0  at the initial block
        F(:, 1:r) = A1_apply(F(:, 1:r)) .* T + v0 * ones(1, r);
    end
    mv = mv + r;
    r2 = 2*r;
    % --------------- Scaling and recovering loop -------------------
    if r == 1 || emptyv0 || emptyV
        krT  = T; 
        krmu = mu;
    else
        krT  = [T T];
        krmu = [mu mu];
    end
    krT = krT / st;
  
    for j = 1:(st - 1)
        EXP = F;
        k   = 0;
        c1  = inf;
        c2  = norm(F, nrm);

        while c1 + c2 > tol * norm(EXP, nrm)
            c1 = c2; k = k + 1;
            F  = A_apply(F);
            F  = F .* (krT / k);
            mv = mv + r2;
            EXP = EXP + F;
            c2  = norm(F, nrm);
        end

        EXP = EXP .* krmu;

        if ~emptyV
            S = S * expJs;
            if emptyv0
                idx = idx_last;
            else
                idx = idx_both;
            end
            F = EXP + S(:, idx);
        end

        if ~emptyv0
            F(:, 1:r) = EXP(:, 1:r);
        end
    end
    %------------ Output shape -------------
    if emptyV  
        W = F;
    elseif emptyv0
        W = F*Alpha;
    else
        W = F(:,1:r) + F(:,r+1:2*r) * Alpha;
    end
end

