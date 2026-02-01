function g = ftGaussianProc1(N, delta, PSDFcn)
% FTGAUSSIANPROC1 Synthesize a 1D Gaussian random process using the FT method
%
% Usage: g = ftGaussianProc1(N, delta, PSDFcn)
%
% Inputs:
%   N      - Number of samples
%   delta  - Spatial sampling interval [m]
%   PSDFcn - Handle to a function that accepts a frequency vector [1/m]
%            and returns the Power Spectral Density
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ift.m) copyright SPIE.

df = 1/(N*delta);   % frequency grid spacing [1/m]
f = (-N/2 : N/2-1) * df;

PSD = PSDFcn(f);
PSD = PSD(:); % ensure that PSD is a column
del_f = 1/(N*delta);   % frequency grid spacing [1/m]
% random draws of Fourier coefficients:
cn = complex(randn(N, 1), randn(N, 1)) .* sqrt(PSD*del_f);
% synthesize the process:
g = real(ift2(cn, 1));
end