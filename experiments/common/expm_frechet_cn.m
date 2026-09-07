function [F,L,cn] = expm_frechet_cn(A,E)

[F,L] = expm_frechet(A,E);
[~,cn] = expm_frechet_cond(A);
end