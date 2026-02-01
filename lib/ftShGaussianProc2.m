function [phz_lo, phz_hi] = ftShGaussianProc2(N, delta, PSDFcn)
% FTSHGAUSSIANPROC2 Synthesize a 2D Gaussian random process with subharmonics
%
% Usage: [phz_lo, phz_hi] = ftShGaussianProc2(N, delta, PSDFcn)
%
% Inputs:
%   N      - Number of samples per side
%   delta  - Spatial sampling interval [m]
%   PSDFcn - Handle to a function for the 2D PSD
%
% Outputs:
%   phz_lo - Low-frequency component (subharmonics)
%   phz_hi - High-frequency component (FFT method)

% High-frequency component from FFT method:
phz_hi = ftGaussianProc2(N, delta, PSDFcn);

% Coordinates for the low-frequency screen:
x = (-N/2 : N/2-1) * delta;
[X, Y] = meshgrid(x);
D = N * delta;

% Initialize low-freq component:
phz_lo = zeros(N);

% Subharmonic levels
NP = 3; 
for p = 1:NP
    df_p = 1 / (3^p * D); % frequency grid spacing at this level
    f_p = [-1, 0, 1] * df_p;
    [fx, fy] = meshgrid(f_p);
    
    % Evaluate PSD on the 3x3 grid
    PSD_p = PSDFcn(fx, fy);
    
    % Random draws of 3x3 Fourier coefficients:
    cn = complex(randn(3), randn(3)) .* sqrt(PSD_p) * df_p;
    
    SH = zeros(N);
    % Sum over the 3x3 frequency grid (excluding the center DC component)
    for ii = 1:3
        for jj = 1:3
            if ii == 2 && jj == 2
                continue; % skip DC
            end
            SH = SH + cn(ii, jj) * exp(1i * 2 * pi * (fx(ii, jj) * X + fy(ii, jj) * Y));
        end
    end
    phz_lo = phz_lo + SH;
end
phz_lo = real(phz_lo);
% Subtract mean of the total signal to center it
total_mean = mean(phz_lo(:) + phz_hi(:));
phz_lo = phz_lo - total_mean;
end