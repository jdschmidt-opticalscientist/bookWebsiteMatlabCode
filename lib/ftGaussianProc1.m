function g = ftGaussianProc1(N, delta, PSDFcn, zeroDC)
% FTGAUSSIANPROC1 Synthesize a 1D Gaussian random process using the FT method
%
% Usage: g = ftGaussianProc1(N, delta, PSDFcn, zeroDC)
%
% Inputs:
%   N      - Number of samples
%   delta  - Spatial sampling interval [m]
%   PSDFcn - Handle to a function that accepts a frequency vector [1/m]
%            and returns the Power Spectral Density
%   zeroDC - Flag for removing mean from Fourier coefficients%

% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ift.m) copyright SPIE.

df = 1/(N*delta);   
f = (-N/2 : N/2-1) * df;
PSD = PSDFcn(f);
PSD = PSD(:); 

% Random draws of Fourier coefficients:
cn = complex(randn(N, 1), randn(N, 1)) .* sqrt(PSD*df);

% Spectral Partitioning: Zero out DC if used with subharmonics
if zeroDC
    cn(N/2+1) = 0;
end

% Synthesize the process:
g = real(ift(cn, 1)); % Note: ift, not ift2 for 1-D
end