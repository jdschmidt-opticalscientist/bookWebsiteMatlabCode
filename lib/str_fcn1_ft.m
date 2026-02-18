function [D, lags] = str_fcn1_ft(x)
% STR_FCN1_FT Computes the unbiased 1-D structure function using FFT
%
% Usage: [D, lags] = str_func_unbiased(x)
%
% This function uses the identity:
% D(tau) = mean((x(t+tau) - x(t))^2)
% Expanded: D(tau) = UnbiasedCorr(x^2, 1) + UnbiasedCorr(1, x^2)...
%  - 2*UnbiasedCorr(x, x)
%
% Input:
%   x    - 1D signal
%
% Outputs:
%   D    - Unbiased structure function
%   lags - Vector of time lags
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

N = length(x);
x = x(:);

% The minimum padding to avoid circular wrap-around
L = 2*N - 1;

% Component A: Cross-term sum(x_i * x_{i+tau})
X = fft(x, L);
corr_raw = real(ifft(X .* conj(X)));

% Component B: Sum of squares in the overlapping windows
x2 = x.^2;
X2 = fft(x2, L);
Ones = fft(ones(N, 1), L);

sum_x_sq_left  = real(ifft(X2 .* conj(Ones)));
sum_x_sq_right = real(ifft(Ones .* conj(X2)));

% Truncate to positive lags
m = (N:-1:1)';
s_left  = sum_x_sq_left(1:N);
s_right = sum_x_sq_right(1:N);
c_raw   = corr_raw(1:N);

% D(tau) = [Sum(x_left^2) + Sum(x_right^2) - 2*Sum(x_left*x_right)] / (N-tau)
D = (s_left + s_right - 2*c_raw) ./ m;
lags = (0:N-1).';
end