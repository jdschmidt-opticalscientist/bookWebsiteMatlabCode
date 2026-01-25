function [phz_lo, phz_hi] = ftShGaussianProc1(M, dt, psdThFcn)
% FTSHGAUSSIANPROC1 Synthesize a 1D Gaussian random process with subharmonics
%
% Usage: [phz_lo, phz_hi] = ftShGaussianProc1(M, dt, psdThFcn)
%
% This function uses the FFT-based method for high frequencies and the
% subharmonic method to compensate for under-sampled low frequencies.
%
% Inputs:
%   M        - Number of samples
%   dt       - Spatial sampling interval [m]
%   psdThFcn - Handle to a function for the PSD (Power Spectral Density)
%
% Outputs:
%   phz_lo   - Low-frequency component (subharmonics)
%   phz_hi   - High-frequency component (FFT method)
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ftGaussianProc1.m) copyright SPIE.

% high-frequency screen from FFT method:
phz_hi = ftGaussianProc1(M, dt, psdThFcn);

% subharmonics:
x = (-M/2 : M/2-1).' * dt;
D = M*dt;
% initialize low-freq screen:
phz_lo = zeros(M, 1);
% loop over frequency grids with spacing 1/(3^p*D)
NP = 5;
for p = 1:NP
    % setup the PSD:
    del_f = 1 / (3^p*D); % frequency grid spacing [1/m]
    fx = linspace(-0.5, 0.5, 3) * del_f;
    PSD_phi = psdThFcn(fx);
    %PSD_phi(2) = 0;
    % random draws of Fourier coefficients:
    cn = complex(randn(1,3), randn(1,3)) .* sqrt(PSD_phi*del_f);
    SH = zeros(M, 1);
    % loop over frequencies on this grid:
    for ii = 1:3
        SH = SH + cn(ii) * exp(1i*2*pi*fx(ii)*x);
    end
    phz_lo = phz_lo + SH;   % accumulate subharmonics
end
phz_lo = real(phz_lo) - mean(real(phz_lo(:)));
end