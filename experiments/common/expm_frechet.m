function [F,L] = expm_frechet(A,E)
%EXPM_FRECHET   Matrix exponential and its Frechet derivative.
%   [F,L] = EXPM(X,E) is the matrix exponential and its Frechet
%   derivative L at X in the direction E.

%   Based on EXPM: $Revision: 5.10.4.6 $  $Date: 2005/11/18 14:15:53 $

wantL = (nargin == 2);

[m_vals, theta, ell, classA] = expmchk; % Initialization
if wantL, theta = ell; end % Frechet parameters determine scaling.
normA = norm(A,1);

if normA <= theta(end)
    % no scaling and squaring is required.
    for i = 1:length(m_vals)
        if normA <= theta(i)
            if wantL
               [F,L] = PadeApproximantOfDegree(m_vals(i));
%               normA, i, theta(i), mym = m_vals(i)
            else
               F = PadeApproximantOfDegree(m_vals(i));
            end
            break;
        end
    end
else
    [t s] = log2(normA/theta(end));
    s = s - (t == 0.5);    % Adjust s if normA/theta(end) is a power of 2.
    A = A/2^s; if wantL, E = E/2^s; end % Scaling
    if wantL
       [F,L] = PadeApproximantOfDegree(m_vals(end));
    else
       F = PadeApproximantOfDegree(m_vals(end));
    end
    for i = 1:s
        if wantL, L = F*L + L*F; end
        F = F*F;           % Squaring
    end
end
% End of expm

%%%%Nested Functions%%%%
    function [m_vals, theta, ell, classA] = expmchk
        %EXPMCHK Check the class of input A and
        %    initialize M_VALS and THETA accordingly.
        classA = class(A);
        switch classA
            case 'double'
                m_vals = [3 5 7 9 13];
                % theta_m for m=1:13.  Values for exp computation.
                theta = [%3.650024139523051e-008
                         %5.317232856892575e-004
                          1.495585217958292e-002  % m_vals = 3
                         %8.536352760102745e-002
                          2.539398330063230e-001  % m_vals = 5
                         %5.414660951208968e-001
                          9.504178996162932e-001  % m_vals = 7
                         %1.473163964234804e+000
                          2.097847961257068e+000  % m_vals = 9
                         %2.811644121620263e+000
                         %3.602330066265032e+000
                         %4.458935413036850e+000
                          5.371920351148152e+000];% m_vals = 13
                % ell_m for m=1:13.  Values for Frechet deriv. computation.
                 ell =  [%2.107342425540526e-008
                         %3.555847927942764e-004
                          1.081338577784837e-002  % m_vals = 3
                         %6.486276265575182e-002
                          1.998063206978949e-001  % m_vals = 5
                         %4.373149570114990e-001
                          7.834608472962045e-001  % m_vals = 7
                         %1.234629393622150e+000
                          1.782448623969279e+000  % m_vals = 9
                         %2.416804616122239e+000
                         %3.127459108657023e+000
                         %3.904829673913844e+000
                          4.740307543766806e+000];  % m_vals = 13
            case 'single'
                m_vals = [3 5 7];
                % theta_m for m=1:7.
                theta = [%8.457278879935396e-004
                         %8.093024012430565e-002
                          4.258730016922831e-001  % m_vals = 3
                         %1.049003250386875e+000
                          1.880152677804762e+000  % m_vals = 5
                         %2.854332750593825e+000
                          3.925724783138660e+000];% m_vals = 7
                % ell vector not yet computed for 'single'.
            otherwise
                error('MATLAB:expm:inputType','Input must be single or double.')
        end
    end

    function [F,L] = PadeApproximantOfDegree(m)
        %PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
        %   F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
        %   Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
        %   Series are evaluated in decreasing order of powers, which is
        %   in approx. increasing order of maximum norms of the terms.

        n = length(A);
        c = getPadeCoefficients;

        % Evaluate Pade approximant.
        switch m

            case {3, 5, 7, 9}

                % Pade for exp.
                m2 = (m+1)/2;
                Apowers = cell(m2,1);
                Apowers{1} = eye(n,classA);
                Apowers{2} = A*A;
                for j = 3:m2
                    Apowers{j} = Apowers{j-1}*Apowers{2};
                end
                U = zeros(n,classA); V = zeros(n,classA);

                for j = m+1:-2:2
                    U = U + c(j)*Apowers{j/2};
                end
                Usave = U;
                U = A*U;
                for j = m:-2:1
                    V = V + c(j)*Apowers{(j+1)/2};
                end

                if wantL
                   % NB: I'm not happy with the indexing below.
                   %     Better way to write it?
                   % Pade for Frechet derivative.
                   m3 = (m-1)/2;
                   M = cell(m3,1);
                   M{1} = E*A + A*E;
                   for j = 2:m3
                       M{j} = M{1}*Apowers{j} + Apowers{2}*M{j-1};
                   end
                   Lu = zeros(n,classA); Lv = zeros(n,classA);
                   for j = m+1:-2:4
%                       Lu = Lu + c(j)*M{j/2};
                       Lu = Lu + c(j)*M{j/2-1};
                   end
                   Lu = A*Lu + E*Usave;
                   for j = m:-2:3
%                       Lv = Lv + c(j)*M{(j+1)/2};
                       Lv = Lv + c(j)*M{(j-1)/2};
                   end
                end

            case 13

                % For optimal evaluation need different formula for m >= 12.

%                 A2 = A*A; A4 = A2*A2; A6 = A2*A4;
%                 UU = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) ...
%                     + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n,classA) );
%                 VV = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) ...
%                     + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n,classA);

                % Pade for exp.
                A2 = A*A; A4 = A2*A2; A6 = A2*A4;
                W1 = c(14)*A6 + c(12)*A4 + c(10)*A2;
                W2 = c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n,classA);
                W = A6*W1 + W2;
                U = A*W;
                Z1 = c(13)*A6 + c(11)*A4 + c(9)*A2;
                Z2 = c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n,classA); 
                V = A6*Z1 + Z2;

%                m,UV_errs = [norm(U-UU,1)/norm(UU)  norm(V-VV,1)/norm(VV)]

                if wantL
                   % Pade for Frechet derivative.
                   M2 = E*A + A*E;
                   M4 = M2*A2 + A2*M2;
                   M6 = M4*A2 + A4*M2;
                   Lw1 = c(14)*M6 + c(12)*M4 + c(10)*M2;
                   Lw2 = c(8)*M6 + c(6)*M4 + c(4)*M2;
                   Lw = A6*Lw1 + M6*W1 + Lw2;
                   Lu = A*Lw + E*W;
                   Lz1 = c(13)*M6 + c(11)*M4 + c(9)*M2;
                   Lz2 = c(7)*M6 + c(5)*M4 + c(3)*M2;
                   Lv = A6*Lz1 + M6*Z1 + Lz2;
                end

        end

        F = (-U+V)\(U+V);
        if wantL
           [XL,XU] = lu(-U+V);
           F = XU\(XL\(U+V));
           L = XU\(XL\(Lu+Lv + (Lu-Lv)*F)); 
        end

        function c = getPadeCoefficients
            % GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
            %    C = GETPADECOEFFICIENTS returns coefficients of numerator
            %    of [M/M] Pade approximant, where M = 3,5,7,9,13.
            switch m
                case 3
                    c = [120, 60, 12, 1];
                case 5
                    c = [30240, 15120, 3360, 420, 30, 1];
                case 7
                    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
                case 9
                    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
                         2162160, 110880, 3960, 90, 1];
                case 13
                    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
                         1187353796428800,  129060195264000,   10559470521600, ...
                         670442572800,      33522128640,       1323241920,...
                         40840800,          960960,            16380,  182,  1];
            end
        end
    end
%%%%Nested Functions%%%%
end
