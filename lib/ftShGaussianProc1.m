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

% 1. High-frequency component (FFT method)
% Zero the DC bin to allow subharmonics to handle the low-frequency limit
zeroDC = true;
phz_hi = ftGaussianProc1(M, dt, psdThFcn, zeroDC);

% 2. Low-frequency component (Subharmonics)
x = (-M/2 : M/2-1).' * dt;
D = M*dt;
phz_lo = zeros(M, 1);
NP = 3; % Typically 3 levels are sufficient for Gaussian processes

for p = 1:NP
    df_p = 1 / (3^p * D); 
    fx = [-1, 0, 1] * df_p;
    PSD_p = psdThFcn(fx);
    
    cn = complex(randn(1, 3), randn(1, 3)) .* sqrt(PSD_p * df_p);
    
    SH = zeros(M, 1);
    for ii = 1:3
        % SPECTRAL PARTITIONING:
        % Only include the center bin (ii=2) at the very last (smallest) 
        % level to represent the ensemble mean variance exactly once.
        if (ii == 2) && (p < NP)
            continue; 
        end
        SH = SH + cn(ii) * exp(1i * 2 * pi * fx(ii) * x);
    end
    phz_lo = phz_lo + SH;
end
phz_lo = real(phz_lo);
% DO NOT subtract mean here to preserve absolute variance
end