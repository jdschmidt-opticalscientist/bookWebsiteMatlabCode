% checkFTSHScreens.m
% Comparison of FT-only (High-Freq) vs. Subharmonic phase screens
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

clear; close all;

% --- Parameters ---
N = 128;            % Grid size
D_aper = 0.8;       % Aperture diameter [m]       
r0 = 0.15;          % Fried parameter [m]
L0 = 100;           % Outer scale [m]
l0 = 0.01;          % Inner scale [m]
dx = 2 * D_aper / N; % Grid spacing [m]
NR = 500;           % Number of realizations

% --- Setup Grids ---
Npad = 2 * N;
idxCenter = (1:N) + N/2; 
lagsPad = (-Npad/2 : Npad/2-1) * dx;
[tauXPad, tauYPad] = meshgrid(lagsPad);

% --- Aperture Mask ---
x = (-N/2 : N/2-1) * dx;
[xx, yy] = meshgrid(x);
mask = double(sqrt(xx.^2 + yy.^2) < D_aper/2);
maskPad = zeros(Npad, Npad);
maskPad(idxCenter, idxCenter) = mask;

% --- Atmospheric PSD ---
fm = 5.92 / l0 / (2*pi); f0 = 1 / L0;
psdAtm = @(fx, fy) 0.023 * r0^(-5/3) * exp(-(sqrt(fx.^2+fy.^2)/fm).^2) ...
    ./ (fx.^2 + fy.^2 + f0^2).^(11/6);

% --- Theory ---
[~, r] = cart2pol(tauXPad, tauYPad);
D_th = 6.88 * (r / r0).^(5/3); 

% --- Initialize Accumulators ---
avgD_SH = zeros(Npad, Npad);
avgD_FT = zeros(Npad, Npad);

fprintf('Running %i realizations...\n', NR);
for idx = 1:NR
    % 1. Generate Components
    [phz_lo, phz_hi] = ftShGaussianProc2(N, dx, psdAtm);
    phzSH = phz_lo + phz_hi;
    phzFT = phz_hi; 
    
    % Remove piston over aperture
    phzSH = phzSH - mean(phzSH(mask==1));
    phzFT = phzFT - mean(phzFT(mask==1));
    
    % 2. Padding
    shPad = zeros(Npad, Npad);
    shPad(idxCenter, idxCenter) = phzSH .* mask;
    ftPad = zeros(Npad, Npad);
    ftPad(idxCenter, idxCenter) = phzFT .* mask;
    
    % 3. Compute Unbiased Stats
    [D_sh_sim, ~] = str_fcn2_ft(shPad, maskPad, dx);
    [D_ft_sim, ~] = str_fcn2_ft(ftPad, maskPad, dx);
    
    % 4. Accumulate
    avgD_SH = avgD_SH + D_sh_sim / NR;
    avgD_FT = avgD_FT + D_ft_sim / NR;
end

% --- Radial Averaging ---
[rvals, D_sh_radial] = azimuthal_average(tauXPad, tauYPad, avgD_SH);
[~, D_ft_radial] = azimuthal_average(tauXPad, tauYPad, avgD_FT);
[~, D_th_radial] = azimuthal_average(tauXPad, tauYPad, D_th);

% --- Plotting Figure 1: Statistics ---
f1 = figure(1);
set(f1, 'Position', [100, 500, 1000, 450]);
t1 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(rvals, D_th_radial, 'k-', 'LineWidth', 2); hold on;
plot(rvals, D_sh_radial, 'r--', 'LineWidth', 1.5);
plot(rvals, D_ft_radial, 'b:', 'LineWidth', 1.5);
grid on; xlim([0, D_aper]);
xlabel('Radial Lag r [m]');
ylabel('D_\phi(r) [rad^2]');
title('Linear Scale');
legend('Theory', 'Subharmonics (lo+hi)', 'FT-Only (hi)', 'Location', 'NorthWest');

nexttile;
loglog(rvals, D_th_radial, 'k-', 'LineWidth', 2); hold on;
loglog(rvals, D_sh_radial, 'r--', 'LineWidth', 1.5);
loglog(rvals, D_ft_radial, 'b:', 'LineWidth', 1.5);
grid('on');
xlim([dx, D_aper]);
ylim([1e-1, 1e2]);
xlabel('Radial Lag r [m]');
ylabel('D_\phi(r) [rad^2]');
title('Log-Log Scale');

title(t1, 'Statistical Impact of Subharmonics');
exportgraphics(f1, 'subharmonicStrFcn.png');

% --- Plotting Figure 2: Visual Realization ---
f2 = figure(2);
set(f2, 'Position', [100, 100, 1000, 450]);
t2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(x, x, phzFT);
xlabel('x [m]');
ylabel('y [m]');
axis('image', 'xy');
colorbar;
title('FT-Only Phase Screen (phz\_hi)');

nexttile;
imagesc(x, x, phzSH);
xlabel('x [m]');
ylabel('y [m]');
axis('image', 'xy');
colorbar;
title('Subharmonic Phase Screen (phz\_lo + phz\_hi)');

title(t2, 'Single Realization Comparison');
exportgraphics(f2, 'subharmonicRealization.png');