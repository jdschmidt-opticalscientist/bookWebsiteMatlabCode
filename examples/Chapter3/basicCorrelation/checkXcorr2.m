% checkXcorr2.m Supplemental code for "Numerical Simulation of Optical 
% Wave Propagation with Examples in MATLAB"
% 
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ft2.m) copyright SPIE.
clear variables; close all; clc;

N = 64;           % number of grid points per side
L = 5.0;          % grid size [m]
dx = L/N;         % grid spacing [m]
x = (-N/2 : N/2-1) * dx;
[X, Y] = meshgrid(x);

% set up theoretical covariance:
w_sig = 8*dx;     % width parameter for Gaussian covariance [m]
varTh = 1.5;      % theoretical variance
lags = (-(N-1) : (N-1)) * dx;
[tauX, tauY] = meshgrid(lags);
% Theoretical 2-D Covariance:
corrTh2D = varTh * exp(-pi * (tauX.^2 + tauY.^2) / w_sig^2);

% set up theoretical PSD function handle:
psdThFcn = @(fx, fy) varTh * w_sig^2 * exp(-pi * (fx.^2 + fy.^2) * w_sig^2);

% Ensemble parameters:
NR = 500; % number of random draws
Npad = 2*N;
avgCorrUnbiased = zeros(Npad, Npad); 

% Use a simple square mask covering the whole grid for this check
mask = ones(N); 
areamask = sum(mask(:)) * dx^2;

% For FT-based unbiasing:
df = 1/(Npad * dx);

% Pre-calculate mask autocorrelation (the unbiasing denominator)
maskPad = zeros(Npad, Npad);
maskPad(1:N, 1:N) = mask; 
W = ft2(maskPad, dx);
maskCorr = ift2(abs(W).^2, df);

% Threshold to avoid division by zero at edges
idxValid = abs(maskCorr) >= dx^2 / areamask;

fprintf('Running ensemble of %i realizations...\n', NR);

for idx = 1 : NR
    % 1. generate 2-D random process realization:
    [phz_lo, phz_hi] = ftShGaussianProc2(N, dx, psdThFcn);
    g = phz_lo + phz_hi;
    
    % 2. remove mean of the realization to compute covariance:
    g = g - mean(g(:));
    
    % 3. Apply mask and manually zero-pad:
    gPad = zeros(Npad, Npad);
    gPad(1:N, 1:N) = g .* mask;
    
    % 4. FT-based correlation (Correlation Theorem):
    G = ft2(gPad, dx);
    rawCorr = ift2(abs(G).^2, df);
    
    % 5. unbias using the mask autocorrelation:
    unbiasedRealization = zeros(Npad, Npad);
    unbiasedRealization(idxValid) = rawCorr(idxValid) ./ maskCorr(idxValid);
    
    % accumulate ensemble average:
    avgCorrUnbiased = avgCorrUnbiased + unbiasedRealization / NR;
end

% --- DATA ANALYSIS & PLOTTING ---

% extract 1-D slice for plotting comparison (along tau_y = 0)
centerIdx = N + 1;
sliceRange = (centerIdx - (N-1)) : (centerIdx + (N-1));
corrSlice = real(avgCorrUnbiased(centerIdx, sliceRange));

% Define grids for radial averaging
lagsPad = (-Npad/2 : Npad/2-1) * dx;
[tauXPad, tauYPad] = meshgrid(lagsPad);

% Perform azimuthal average
[rvals, radialAvg] = azimuthal_average(tauXPad, tauYPad, real(avgCorrUnbiased));

% Theoretical radial function for comparison
corrThRadial = varTh * exp(-pi * rvals.^2 / w_sig^2);

% FIGURE 1: Realization and 1-D Slice
f1 = figure(1); clf;
set(f1, 'OuterPosition', [200 200 1000 450]);
tlayout1 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(x, x, g);
axis image xy; colorbar;
xlabel('x [m]'); ylabel('y [m]');
title('Final Signal Realization');

nexttile;
plot(lags, corrTh2D(N, :), 'k-', 'LineWidth', 2); hold on;
plot(lags, corrSlice, 'r--', 'LineWidth', 1.5);
grid on; xlim([-L L]/2);
xlabel('Lag \Delta x [m]');
ylabel('Autocovariance');
legend('Theory', '2-D Unbiased (Slice)', 'Location', 'NorthEast');
title(sprintf('Ensemble Average (%i draws)', NR));

% FIGURE 2: Azimuthal Average Comparison
f2 = figure(2); clf;
plot(rvals, corrThRadial, 'k-', 'LineWidth', 2); hold on;
plot(rvals, radialAvg, 'r--', 'LineWidth', 1.5);
grid on; xlim([0 L]/2);
xlabel('Radial Lag r [m]');
ylabel('Autocovariance');
legend('Theory', 'Azimuthal Average');
title('Radial Statistics Verification');

% export figures for use in the HTML article:
exportgraphics(f1, 'checkXcorr2D_Results.png');
exportgraphics(f2, 'checkXcorr2D_Radial.png');