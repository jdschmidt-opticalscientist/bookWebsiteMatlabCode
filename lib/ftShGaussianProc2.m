function [phz_lo, phz_hi] = ftShGaussianProc2(N, delta, PSDFcn)
% FTSHGAUSSIANPROC2 Synthesize a 2D Gaussian random process with subharmonics
%
% Usage: [phz_lo, phz_hi] = ftShGaussianProc2(N, delta, PSDFcn)
%
% This function uses spectral partitioning to ensure absolute variance 
% matches analytical theory. The high-frequency component (FFT) has its 
% DC component zeroed out, and the true DC power is sampled only at the 
% deepest subharmonic level to prevent double-counting.
%
% Inputs:
%   N      - Number of samples per side
%   delta  - Spatial sampling interval [m]
%   PSDFcn - Handle to a function for the 2D PSD
%
% Outputs:
%   phz_lo - Low-frequency component (subharmonics)
%   phz_hi - High-frequency component (FFT method)
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

% 1. High-frequency component from FFT method
% Pass 'true' to zeroDC to prevent overlap with subharmonic DC power
zeroDC = true;
phz_hi = ftGaussianProc2(N, delta, PSDFcn, zeroDC);

% 2. Low-frequency component using subharmonics
x = (-N/2 : N/2-1) * delta;
[X, Y] = meshgrid(x);
D = N * delta;
phz_lo = zeros(N);

% Number of subharmonic levels
NP = 3; 

for p = 1:NP
    % Frequency grid spacing at this level
    df_p = 1 / (3^p * D); 
    f_p = [-1, 0, 1] * df_p;
    [fx, fy] = meshgrid(f_p);
    
    % Evaluate PSD on the 3x3 subharmonic grid
    PSD_p = PSDFcn(fx, fy);
    
    % Random draws of 3x3 Fourier coefficients scaled by frequency spacing
    cn = complex(randn(3), randn(3)) .* sqrt(PSD_p) * df_p;
    
    SH = zeros(N);
    % Sum over the 3x3 frequency grid
    for ii = 1:3
        for jj = 1:3
            % SPECTRAL PARTITIONING:
            % Only include the (0,0) center bin at the very last (smallest) 
            % level to represent the ensemble mean variance exactly once.
            if (ii == 2 && jj == 2) && (p < NP)
                continue; 
            end
            
            % Add this frequency component to the spatial screen
            SH = SH + cn(ii, jj) * exp(1i * 2 * pi * (fx(ii, jj) * X + fy(ii, jj) * Y));
        end
    end
    phz_lo = phz_lo + SH;
end

% Ensure the result is real-valued
phz_lo = real(phz_lo);

% Note: We avoid subtracting the ensemble mean here. While individual 
% realizations have non-zero means, forcing a zero-mean on finite grids 
% biases the variance and autocorrelation at large lags.
end