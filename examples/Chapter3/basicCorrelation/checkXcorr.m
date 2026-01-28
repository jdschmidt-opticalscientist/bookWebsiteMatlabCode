% checkXcorr.m Supplemental code for "Numerical Simulation of Optical 
% Wave Propagation with Examples in MATLAB"
% 
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ft.m) copyright SPIE.

clear variables; close all; clc;

M = 32; % number of grid points
T = 2.7;  % grid size [m]
dt = T/M; % grid spacing [m]
t = (-M/2 : M/2-1) * dt;
mLags = (-(M-1) : (M-1)); % lag index for xcorr
t2 = mLags * dt; % sample times [s]

% set up theoretical covariance:
w = 10*dt; % width parameter for Gaussian covariance [m]
varTh = 2.3; % variance
corrTh = varTh * exp(-pi*t2.^2/w^2); % covariance

% set up theoretical PSD:
df = 1/T;   % frequency grid spacing [1/m]
f = (-M/2 : M/2-1) * df;
psdThFcn = @(F) varTh * w*exp(-pi*F.^2*w^2);
psdTh = psdThFcn(f);
vThPSD = trapz(f, psdTh); % check PSD's variance

NR = 5000; % number of random draws

% allocate space for correlation variables:
rNone = zeros(1, 2*M-1);
rBiased = zeros(1, 2*M-1);
rUnbiased = zeros(1, 2*M-1);
gMean = 0;
gMeanSqr = 0;
% for FT-based calculations:
gPad = zeros(1, 2*M); % zero-padded array for g
dfBig = 1/(2*T);   % frequency grid spacing for double-size g [1/m]
idxFill = (-M/2 : M/2-1) + M+1; % indices of grBig to fill
mLagsFT = (-M : M-1); % lag index for FT-based correlation
corrGBig = zeros(1, 2*M);
for idx = 1 : NR
    % generate random process:
    [phz_lo, phz_hi] = ftShGaussianProc1(2*M, dt, psdThFcn);
    g = phz_lo + phz_hi;
    g = g(1:M).';
    
    % compute auto-correlation with xcorr & various options:
    rNone = rNone + xcorr(g, 'none')/NR; % default
    rBiased = rBiased + xcorr(g, 'biased')/NR;
    rUnbiased = rUnbiased + xcorr(g, 'unbiased')/NR;
    
    % use correlation theorem; be sure to pad with zeros:
    gPad(idxFill) = g; % fill center of gPad
    ftGBig = ft(gPad, dt); % g in frequency domain
    % compute auto-correlation:
    corrGBig = corrGBig + ift(abs(ftGBig).^2, dfBig)/NR;
end

% plots

% different normalizations:
f1 = figure(1); clf;
set(f1, 'OuterPosition', [672 128 920 835]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile([1 2]); % span two columns
% plot the random draw:
plot(t, g, 'k', 'LineWidth', 1.5);
grid on;
xlabel({'Time [sec]'; '(a)'});
ylabel('Signal Realization');
nexttile;
plot(t2, corrTh, 'k', t2, rNone/M, 'r--', t2, rBiased, 'g:', ...
    'LineWidth', 1.5);
grid on;
xlabel({'Time Lag [s]'; '(b)'});
ylabel('xcorr biased');
legend('Theory', 'xcorr none / M', 'xcorr biased', ...
    'location', 'NorthEast');
% unbiased correlation:
nexttile;
plot(t2, corrTh, 'k', t2, rUnbiased, 'r--', ...
    t2, rNone./(M-abs(mLags)), 'g:', 'LineWidth', 1.5);
grid on;
xlabel({'Time Lag [s]'; '(c)'});
ylabel('xcorr Unbiased');
legend('Theory', 'xcorr unbiased', 'xcorr none / (M - |m|)', ...
    'location', 'NorthEast');

% export figure to file in PNG format:
exportgraphics(f1, 'checkXcorr.png');

% unbiased with FT:
f2 = figure(2); clf;
plot(t2, corrTh, 'k', t2, rUnbiased, 'r--', ...
    mLagsFT*dt, corrGBig/T*M./(M-abs(mLagsFT)), 'g-.', 'LineWidth', 1.5);
grid('on');
xlabel('Time Lag [s]');
ylabel('Auto-Correlation');
legend('Theory', 'xcorr unbiased', 'FT Corr Pad');

% export figure to file in PNG format:
exportgraphics(f2, 'xCorrFT.png');