% checkStrFcn2.m
% Verification of 2-D unbiased structure function calculation
%
% This script compares:
% 1. Direct unbiased structure function from str_fcn2_ft.m
% 2. Structure function derived from covariance: 2*sigma^2 - 2*gamma(r)
% 3. Theoretical Gaussian structure function
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

clear; close %clc;

% --- Parameters ---
N = 128;            % Grid size
dx = 1/N;           % Grid spacing [m]
w_sig = 0.15;       % Covariance width [m]
var_sig = 2.0;      % Process variance
NR = 1000;          % Number of realizations
D = 0.8;            % Aperture diameter [m]

% --- Setup Grids ---
Npad = 2 * N;
idxCenter = (1:N) + N/2; 
lagsPad = (-Npad/2 : Npad/2-1) * dx;
[tauXPad, tauYPad] = meshgrid(lagsPad);

% --- Aperture Mask ---
x = (-N/2 : N/2-1) * dx;
[xx, yy] = meshgrid(x);
mask = double(sqrt(xx.^2 + yy.^2) < D/2);
maskPad = zeros(Npad, Npad);
maskPad(idxCenter, idxCenter) = mask;

% --- Theory ---
psdThFcn = @(fx, fy) var_sig * (2*pi*w_sig^2) * exp(-2*pi^2 * w_sig^2 * (fx.^2 + fy.^2));
gamma_th = var_sig * exp(-(tauXPad.^2 + tauYPad.^2) / (2 * w_sig^2));
D_th = 2 * var_sig * (1 - exp(-(tauXPad.^2 + tauYPad.^2) / (2 * w_sig^2)));

% --- Initialize Accumulators ---
avgD = zeros(Npad, Npad);
avgDBiased = zeros(Npad, Npad); % New accumulator
avgCorr = zeros(Npad, Npad);

fprintf('Running realizations...\n');
for idx = 1:NR
    % 1. Generate Gaussian Process
    [phz_lo, phz_hi] = ftShGaussianProc2(N, dx, psdThFcn);
    g = phz_lo + phz_hi;
    g = g - mean(g(mask==1)); 
    
    % 2. Padding
    gPad = zeros(Npad, Npad);
    gPad(idxCenter, idxCenter) = g .* mask;
    
    % 3. Compute Stats
    % Get both the unbiased D and the raw maskCorr (denominator)
    [D_sim, maskCorr] = str_fcn2_ft(gPad, maskPad, dx);
    
    % Calculate the biased (raw) version
    % This is the raw numerator from the structure function logic
    % Consistent with: D_biased = D_unbiased .* (maskCorr / maskCorr(origin))
    D_biased_sim = D_sim .* (maskCorr ./ max(maskCorr(:)));
    
    [unbiasedCorr, ~] = corr2_ft(gPad, gPad, maskPad, dx);
    
    % 4. Accumulate
    avgD = avgD + D_sim / NR;
    avgDBiased = avgDBiased + D_biased_sim / NR;
    avgCorr = avgCorr + unbiasedCorr / NR;
end

% --- Radial Averaging ---
[rvals, D_radial] = azimuthal_average(tauXPad, tauYPad, avgD);
[~, D_biased_radial] = azimuthal_average(tauXPad, tauYPad, avgDBiased);
[~, G_radial] = azimuthal_average(tauXPad, tauYPad, avgCorr);
[~, D_th_radial] = azimuthal_average(tauXPad, tauYPad, D_th);

% Derived D from Gamma: 2*Gamma(0) - 2*Gamma(tau)
D_from_G = 2 * (G_radial(1) - G_radial);

% --- Plotting ---
f1 = figure(1);
plot(rvals, D_th_radial, 'k-', 'LineWidth', 2); hold on;
plot(rvals, D_radial, 'r--', 'LineWidth', 1.5);
plot(rvals, D_from_G, 'b:', 'LineWidth', 1.5);
plot(rvals, D_biased_radial, 'g-.', 'LineWidth', 1.5); % Biased curve
grid on;
xlim([0, D/2]); 
xlabel('Radial Lag r [m]');
ylabel('Structure Function D_\phi(r)');
legend('Theory', 'Direct Unbiased D', '2[\Gamma(0) - \Gamma(r)]', 'Biased D', ...
    'Location', 'SouthEast');
title(sprintf('2-D Structure Function Verification (N_R=%i)', NR));

exportgraphics(f1, 'checkStrFcn2D_Results.png');