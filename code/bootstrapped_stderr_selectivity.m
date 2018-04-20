function [E, med] = bootstrapped_stderr_selectivity(X, n_smps)

% Calculates standard errors for a selectivity index using matlab's bootci function.
% 
% X = [randn(100,1)+2, randn(100,1)+1];
% [E, med] = bootstrapped_stderr_selectivity(X, n_smps)

selectivityfn = @(X)selectivity(X);
[E, bootstat] = bootci(n_smps, {selectivityfn, X}, 'alpha', normcdf(-1)*2);
med = median(bootstat);

function s = selectivity(X)

m = mean(X);
s = (m(1) - m(2)) / (m(1) + m(2));