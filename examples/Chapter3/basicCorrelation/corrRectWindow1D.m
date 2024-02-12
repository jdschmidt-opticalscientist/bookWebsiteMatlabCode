% corrRectWindow1D.m

clear variables; clc;

N = 10; % number of samples

w = [0 0 0 1 1 1 1 0 0 0]; % window function

wCorrNneg = zeros(1,N); % allocate space for correlation, just nonnegative lags
% explicit sum for correlation:
for m = 0 : N-1
    wCorrNneg(m+1) = sum(w(m+1:N) .* w(1:N-m));
end
% use symmetry to fill in negative lags:
wCorrLoop = cat(2, wCorrNneg(N:-1:2), wCorrNneg);
% use xcorr to compute the same thing:
wXCorr = xcorr(w, 'none');
% display correlation from loop:
fprintf('wCorrLoop = \n')
fprintf('%1.0f ', wCorrLoop)
fprintf('\n')
% display correlation from xcorr:
fprintf('wXCorr = \n')
fprintf('%1.0f ', wXCorr)
fprintf('\n')