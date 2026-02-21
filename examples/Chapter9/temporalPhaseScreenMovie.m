% temporalPhaseScreenMovie.m
% Test script for evolving phase screens using modified von Karman PSD
% Creates a high-quality MP4 of evolving phase screens for the web.
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

% 1. Simulation Setup
N = 256;       % number of grid points per side
delta = 0.01;  % grid spacing [m]
D = N * delta; % grid size [m]
r0 = 0.1;      % Fried parameter of screens [m]
L0 = 10;       % outer scale [m]
l0 = 0.005;    % inner scale [m]
v_wind = [1.5, 0.3]; % wind velocity [m/s]
v_mag = norm(v_wind); % Wind speed [m/s]
fG = 0.427 * v_mag / r0; % Greenwood frequency [Hz]
dt = min([1/(10*fG), delta/v_mag]); % time step [s]
num_frames = 5/dt;    % 5 seconds of turbulence

% PSD Function Handle (Modified von Karman)
f0 = 1/L0; fm = 5.92 / (2 * pi * l0);
psd_fcn = @(fx, fy) 0.023 * r0^(-5/3) * ...
    (fx.^2 + fy.^2 + f0^2).^(-11/6) .* exp(-(fx.^2 + fy.^2) / fm^2);

% 2. Initialize Video Writer
video_filename = 'turbulence_evolution.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 30;
v.Quality = 95; % High quality for text/web
open(v);

% 3. Initialize Atmosphere
coeffs = ftShGaussianProc2Coeffs(N, delta, psd_fcn);

% 4. Generation Loop
fig = figure('Color', 'w', 'Position', [100, 100, 600, 550]);
x = ((-N/2) : (N/2-1)) * delta;

for k = 1:num_frames
    t = (k-1) * dt;
    phz = ftShGaussianProc2Evolving(coeffs, N, delta, t, v_wind);
    
    % --- High Quality Plotting ---
    imagesc(x, x, phz);
    axis('image', 'xy');
    colormap(parula); % Perceptually uniform
    cb = colorbar;
    ylabel(cb, 'Phase [rad]', 'FontSize', 12);
    
    % Limits should stay constant for video consistency
    clim([-45, 45]); 
    
    title(sprintf('Taylor Frozen Flow: t = %.2f s', t), 'FontSize', 14);
    xlabel('x [m]');
    ylabel('y [m]');
    set(gca, 'FontSize', 10, 'TickDir', 'out');

    % Capture frame and write
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);
fprintf('Movie saved as: %s\n', video_filename);