function g = ftGaussianProc2(N, delta, PSDFcn, zeroDC)
% FTGAUSSIANPROC2 Synthesize a 2D Gaussian random process using the FT method
%
% Usage: g = ftGaussianProc2(N, delta, PSDFcn, zeroDC)
%
% Inputs:
%   N      - Number of samples per side (square grid)
%   delta  - Spatial sampling interval [m]
%   PSDFcn - Handle to a function that accepts two frequency vectors [1/m]
%            and returns the 2D Power Spectral Density
%   zeroDC - Flag for removing mean from Fourier coefficients
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ift2.m) copyright SPIE.

df = 1/(N*delta);   % frequency grid spacing [1/m]
f = (-N/2 : N/2-1) * df;
[fx, fy] = meshgrid(f);

% Evaluate PSD on the 2D grid
PSD = PSDFcn(fx, fy);

% Random draws of Fourier coefficients:
% We multiply by df (sqrt(PSD * df_x * df_y))
cn = complex(randn(N), randn(N)) .* sqrt(PSD) * df;

% make the DC component zero, like if being used with a subharmonic method
% to avoid double-counting the DC component
if zeroDC
    cn(N/2+1, N/2+1) = 0;
end

% Synthesize the process:
% The factor of N^2 or similar is handled inside ift2 if using book definitions
g = real(ift2(cn, 1)); 
end