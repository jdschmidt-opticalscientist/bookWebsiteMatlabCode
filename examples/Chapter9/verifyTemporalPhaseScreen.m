% verifyTemporalPhaseScreen.m
% Computes temporal structure function averaged over space and realizations.
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

% 1. Simulation Parameters
N = 256;            % Number of grid points
delta = 0.002;      % Grid spacing [m]
r0 = 0.1;           % Fried parameter [m]
L0 = 20;            % Outer scale [m]
l0 = 0.005;         % Inner scale [m]
v_wind = [3, 1];    % Wind velocity [m/s]
v_mag = norm(v_wind); % Wind speed [m/s]
fG = 0.427 * v_mag / r0; % Greenwood frequency [Hz]
dt = min([1/(10*fG), delta/v_mag]); % time step [s]

NR = 5;             % Number of realizations
NS = 40;            % Time steps per realization

% 2. Modified von Karman PSD Function
f0 = 1/L0;
fm = 5.92 / (2 * pi * l0);
psd_fcn = @(fx, fy) 0.023 * r0^(-5/3) * ...
    (fx.^2 + fy.^2 + f0^2).^(-11/6) .* exp(-(fx.^2 + fy.^2) / fm^2);

tau = (0:NS-1) * dt; % time lag grid [s]
D_phi_total = zeros(1, NS); % inintialize sum of all structure functions

% 3. Main Realization Loop
fprintf('Starting realizations with spatial averaging...\n');
for idxR = 1:NR
    % Store a 3D cube of phase [N x N x Time]
    phz_cube = zeros(N, N, NS);
    coeffs = ftShGaussianProc2Coeffs(N, delta, psd_fcn);

    for k = 1:NS
        phz_cube(:,:,k) ...
            = ftShGaussianProc2Evolving(coeffs, N, delta, (k-1)*dt, v_wind);
    end

    % Compute Structure Function for this cube
    % Average over all pixels (x,y) and all valid time pairs for each lag
    D_phi_realization = zeros(1, NS); % initialize structure function
    for idx = 1 : N*N % loop over all grid points
        [idxRow, idxCol] = ind2sub([N N], idx); % get row & column
         % Get one grid point's history
        phz_history = squeeze(phz_cube(idxRow, idxCol, :));
        % compute str fcn and accumulate over grid points
        [D_pixel, lags] = str_fcn1_ft(phz_history);
        D_phi_realization = D_phi_realization + D_pixel.';
    end
    % accumulate over realizations
    D_phi_total = D_phi_total + D_phi_realization;
    fprintf('  Realization %d/%d complete.\n', idxR, NR);
end
D_phi_sim = D_phi_total/(NR*N^2); % from sum to average

% 4. Numerical Theory: Modified von Karman
fprintf('Computing numerical theory line...\n');
% D_phi(rho) = 4*pi * integral[ f * PSD(f) * (1 - J0(2*pi*f*rho)) df ]
psd_fcn = @(F) 0.023 * r0^(-5/3) * ...
    (F.^2 + f0^2).^(-11/6) .* exp(-F.^2 / fm^2);
rho = v_mag * tau;
integrand = @(F) 2*pi* F .* psd_fcn(F) .* (1-besselj(0, 2*pi*F.*rho.'));
D_theory_mvk = 2*integral(integrand, 0, Inf, 'ArrayValued', true);

% 5. Comparison Plot
f1 = figure('Color', 'w', 'Position', [100 100 700 500]);
loglog(tau, D_phi_sim, 'bo', 'MarkerSize', 5, 'DisplayName', ...
    'Simulation (Spatial Avg)');
hold on;
loglog(tau(2:end), D_theory_mvk(2:end), 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Theory (Mod. von Karman)');

% Add Kolmogorov for reference
D_kol = 6.88 * (rho / r0).^(5/3);
loglog(tau(2:end), D_kol(2:end), 'k--', ...
    'DisplayName', 'Theory (Kolmogorov)');

grid on;
xlabel('Time Lag \tau [s]'); ylabel('D_\phi(\tau) [rad^2]');
title('Temporal Structure Function: Space-Time Ensemble Average');
legend('Location', 'northwest');

% 6. Coherence Time (tau_0) Calculation
% Standard Kolmogorov analytical value
tau0_kol_theory = 0.314 * r0 / v_mag;

% Find tau_0 numerically for the Simulation
% We find where D_phi(tau) = 1 rad^2 using log-linear interpolation
tau0_sim = exp(interp1(log(D_phi_sim(2:end)), log(tau(2:end)), log(1)));

% Find tau_0 numerically for the MvK Theory
tau0_mvk_theory ...
    = exp(interp1(log(D_theory_mvk(2:end)), log(tau(2:end)), log(1)));

% Display Results
fprintf('\n--- Coherence Time (tau_0) Results ---\n');
fprintf('Kolmogorov Analytical: %.4f s\n', tau0_kol_theory);
fprintf('MvK Numerical Theory:  %.4f s\n', tau0_mvk_theory);
fprintf('Simulation Result:     %.4f s\n', tau0_sim);

% Add to the plot
loglog([tau(1) tau(end)], [1 1], 'k:', ...
    'HandleVisibility', 'off'); % 1 rad^2 line
loglog(tau0_sim, 1, 'ks', 'MarkerFaceColor', 'g', 'DisplayName', ...
    ['Sim \tau_0 = ' num2str(tau0_sim,3) 's']);

exportgraphics(f1, 'temporal_structure_function.png')