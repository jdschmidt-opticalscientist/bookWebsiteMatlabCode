% checkXcorr2.m Supplemental code for "Numerical Simulation of Optical 
% Wave Propagation with Examples in MATLAB"
% 
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ft2.m) copyright SPIE.

clear variables; close all; clc;

% --- Configuration ---
N = 128;          % number of grid points per side
L = 5.0;          % grid size [m]
dx = L/N;         % grid spacing [m]
x = (-N/2 : N/2-1) * dx;
[xx, yy] = meshgrid(x);

% --- Theoretical Covariance & PSD ---
w_sig = 1.6;      % width parameter for Gaussian covariance [m]
varTh = 1.5;      % theoretical variance
psdThFcn = @(fx, fy) varTh * w_sig^2 * exp(-pi * (fx.^2 + fy.^2) * w_sig^2);

% --- Simulation Setup ---
NR = 2000;        % Ensemble size
Npad = 2*N;       
avgCorrNone = zeros(Npad, Npad); 
avgCorrUnbiased = zeros(Npad, Npad); 

% Circular aperture mask
D = 3.0;          
mask = (sqrt(xx.^2 + yy.^2) <= D/2);
areamask = sum(mask(:)) * dx^2;

% --- Unbiasing Denominator (Aperture Autocorrelation) ---
df = 1/(Npad * dx); % frequency spacing [cyc/m]
maskPad = zeros(Npad, Npad);
idxCenter = (1:N) + N/2;
maskPad(idxCenter,idxCenter) = mask; 
W = ft2(maskPad, dx);
maskCorr = ift2(abs(W).^2, df);

% Threshold to avoid noise amplification at aperture edges:
idxValid = abs(maskCorr) >= (dx^2 / areamask);

fprintf('Running ensemble of %i realizations...\n', NR);

for idx = 1 : NR
    % 1. Generate 2-D random process with subharmonics
    [phz_lo, phz_hi] = ftShGaussianProc2(N, dx, psdThFcn);
    g = phz_lo + phz_hi;
    
    % 2. Prepare the padded and masked arrays
    % Note: mPad is defined outside the loop as zeros(Npad, Npad)
    % with the mask in the top-left (or center) NxN block.
    gPad = zeros(Npad, Npad);
    gPad(idxCenter,idxCenter) = g .* mask;
    
    % 3. Use the new function for the unbiased correlation
    % We pass the padded signal and padded mask.
    [unbiasedRealization, rawCorr] = corr2_ft(gPad, gPad, maskPad, dx);
    
    % 4. Accumulate results
    % avgCorrNone stores the "raw" (biased) correlation for comparison
    avgCorrNone = avgCorrNone + rawCorr / NR;
    
    % avgCorrUnbiased stores the unmasked, scaled correlation
    avgCorrUnbiased = avgCorrUnbiased + unbiasedRealization / NR;
end

% --- ANALYSIS ---

% Prep coordinates for radial averaging
lagsPad = (-Npad/2 : Npad/2-1) * dx;
[tauXPad, tauYPad] = meshgrid(lagsPad);

% Radial average of the BIASED data (normalized by mask area)
[rvals, radialBiased] = azimuthal_average(tauXPad, tauYPad, real(avgCorrNone)/areamask);

% Radial average of the UNBIASED data
[~, radialUnbiased] = azimuthal_average(tauXPad, tauYPad, real(avgCorrUnbiased));

% Analytical theory for radial plot
corrThRadial = varTh * exp(-pi * rvals.^2 / w_sig^2);

% --- PLOT ---
f1 = figure(1); clf;
set(f1, 'OuterPosition', [200 200 1000 800]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Signal realization
nexttile([1 2]);
imagesc(x, x, g .* mask); axis image xy; colorbar;
xlabel('x [m]'); ylabel('y [m]');
title('Final Signal Realization within Circular Aperture (a)');

% (b) Biased Radial Average vs Theory
nexttile;
plot(rvals, corrThRadial, 'k-', 'LineWidth', 2); hold on;
plot(rvals, radialBiased, 'r--', 'LineWidth', 1.5);
grid on; xlim([0 L/2]); ylim([0 varTh*1.1]);
xlabel('Radial Lag r [m]'); ylabel('Autocovariance');
title('Biased Radial Average (b)');
legend('Theory', 'Biased (Raw/Area)', 'Location', 'NorthEast');

% (c) Unbiased Radial Average vs Theory
nexttile;
plot(rvals, corrThRadial, 'k-', 'LineWidth', 2); hold on;
plot(rvals, radialUnbiased, 'r--', 'LineWidth', 1.5);
grid on; xlim([0 L/2]); ylim([0 varTh*1.1]);
xlabel('Radial Lag r [m]'); ylabel('Autocovariance');
title('Unbiased Radial Average (c)');
legend('Theory', 'Unbiased', 'Location', 'NorthEast');

exportgraphics(f1, 'checkXcorr2D_Results.png');