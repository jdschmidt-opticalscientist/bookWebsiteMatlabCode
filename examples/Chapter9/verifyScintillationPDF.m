% verifyScintillationPDF.m
% Simulation of a point source through turbulence with PDF verification.
% Compares simulated results to Log-Normal and Gamma-Gamma models.
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

example_pt_source_atmos_setup

analysis_pt_source_atmos_samp

example_pt_source_vac_prop

l0 = 0;     % inner scale [m]
L0 = 1e3;   % outer scale [m]

zt = [0 z];  % propagation plane locations
Delta_z = zt(2:n) - zt(1:n-1);    % propagation distances
% grid spacings
alpha = zt / zt(n);
delta = (1-alpha) * delta1 + alpha * deltan;

% initialize array for phase screens
phz = zeros(N, N, n);
nreals = 500;    % number of random realizations
% initialize arrays for propagated fields,
% aperture mask, and MCF
Uout = zeros(N);
mask = circ(xn/D2, yn/D2, 1);

% PSD Handle for ftShGaussianProc2
fm = 5.92 / l0 / (2*pi); f0 = 1 / L0;
psdAtm = @(fx, fy, r0_val) 0.023 * r0_val^(-5/3) ...
    * exp(-(sqrt(fx.^2+fy.^2)/fm).^2) ...
    ./ (fx.^2 + fy.^2 + f0^2).^(11/6);

% Initialize stats
map = logical(mask);
NMask = sum(map(:));
logAmp = zeros(NMask, nreals);
irr = zeros(NMask, nreals);

fprintf('Running %i realizations...\n', nreals);
for idxreal = 1 : nreals
    % 1. Generate Phase Screens using the new subharmonic function
    phz = zeros(N, N, n);
    for idxscr = 1:n
        % Using your new general function
        [phz_lo, phz_hi] = ftShGaussianProc2(N, delta(idxscr), ...
                                @(fx, fy) psdAtm(fx, fy, r0scrn(idxscr)));
        phz(:,:,idxscr) = phz_lo + phz_hi;
        % mean phase does not affect propagation, so subtract it:
        phz(:,:,idxscr) = phz(:,:,idxscr) - mean(mean(phz(:,:,idxscr)));
    end

    % 2. Turbulent Propagation
    % Pass the combined screens (sg = super-Gaussian absorbing boundary)
    [xn, yn, Uout] ...
        = ang_spec_multi_prop(pt, wvl, delta1, deltan, z, sg.*exp(1i*phz));
    
    % 3. Collimate and Mask
    Uout = Uout .* exp(-1i*pi/(wvl*R)*(xn.^2+yn.^2));
    
    % 4. Store perturbed fields
    % Normalized by vacuum mean log-amplitude
    irr(:,idxreal) = abs(Uout(map)).^2;
    logAmp(:,idxreal) = log(abs(Uout(map))); 
end

mLAV = mean(log(abs(Uvac(map)))); % mean vacuum log-amplitude
logAmp = logAmp - mLAV; % log-amplitude perturbation

% --- PDF Analysis ---
irrMean = mean(irr(:));
I_norm = irr(:) / irrMean; % Normalized Irradiance
I_bins = linspace(0.01, 6, 200);

% 1. Log-Normal Model (Weak Turbulence)
% Note: In weak turbulence, sigma_I^2 is approximately 4*sigma_chi^2
pdf_LogNormal = 1./(I_bins * sqrt(8*pi*rytov)) ...
    .* exp(-(log(I_bins) + 2*rytov).^2 / (8*rytov));

% 2. Gamma-Gamma Model (Strong/Saturated Turbulence)
% Parameters alpha and beta represent large and small scale variances
beta02 = 4 * rytov;
% Andrew & Phillips Ch. 9 Eq. (74)
sigI2 = exp(0.49*beta02/(1+0.56*beta02^(6/5))^(7/6) ...
    + 0.51*beta02/(1+0.69*beta02^(6/5))^(5/6)) - 1;
% from Andrews, Phillips, & Hopen, "Scintillation model for a satellite
% communication link at large zenith angles", Opt. Eng., Vol 39, No. 12,
% pp. 3272-3280 (1999):
% moments:
L = Dz;
z = linspace(0, L, 100);
xi = z/L;
Cn2Vec = repmat(Cn2, [1 100]);
mu0 = trapz(xi, Cn2Vec) * L; % used for large scale, Eq. (11)
mu3 = trapz(xi, Cn2Vec .* (xi.*(1-xi)).^(5/6)) * L; % Eq. (41)
mu4 = trapz(xi, Cn2Vec .* (xi.*(1-xi)).^(-1/3)) * L; % Eq. (50)

% variances:
sig22 = 4 * rytov; % also = beta0^2

% small-scale log-irr variance,  Eq. (46):
S2lnY = 0.51*sig22 ./ (1 + 0.69*sig22.^(6/5)).^(5/6);
% large-scale log-irr variance, Eq. (52):
S2lnX = 0.49*sig22 ...
    ./ (1 + 0.62*(mu3./mu4).^(6/7).*(mu0./mu3).^(6/5).*sig22.^(6/5)).^(7/6);

alph = 1 / (exp(S2lnX) - 1);
bet  = 1 / (exp(S2lnY) - 1);
pdf_GammaGamma = 2*(alph*bet)^((alph+bet)/2) / (gamma(alph)*gamma(bet)) ...
    * I_bins.^((alph+bet)/2-1) .* besselk(alph-bet, 2*sqrt(alph*bet*I_bins));

%% --- Plotting ---
f1 = figure(1);
t = tiledlayout(1,2,'TileSpacing','compact');

% Irradiance Histogram
nexttile;
histogram(I_norm, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8], ...
    'EdgeColor', 'none');
hold('on');
plot(I_bins, pdf_LogNormal, 'k-', 'LineWidth', 2);
plot(I_bins, pdf_GammaGamma, 'r--', 'LineWidth', 2);
xlim([0 6]);
grid('on');
xlabel('Normalized Irradiance I / <I>');
ylabel('Probability Density');
title('Irradiance Statistics');
legend('Simulation', 'Log-Normal', 'Gamma-Gamma');

% Log-Amplitude Histogram
nexttile;
histogram(logAmp(:), 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8], ...
    'EdgeColor', 'none');
hold('on');
chi = linspace(min(logAmp(:)), max(logAmp(:)), 100);
pdfChiR = 1/sqrt(2*pi*rytov) * exp(-(chi+rytov).^2/(2*rytov));
plot(chi, pdfChiR, 'k-', 'LineWidth', 2);
grid('on');
xlabel('Log-Amplitude \chi');
title('Log-Amplitude Statistics');

exportgraphics(f1, 'verifyScintillationPDF.png');