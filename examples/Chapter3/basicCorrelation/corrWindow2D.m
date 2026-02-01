% corrWindow2D.m Supplemental code for "Numerical Simulation of Optical 
% Wave Propagation with Examples in MATLAB"
% 
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ft.m) copyright SPIE.

% define window and spatial grid:
N = 512; % number of grid points per side
w = ones(N); % window that covers the entire grid
w(201:400,301:500) = 0;
dx = 1; % grid spacing [m]
x = (-N/2 : N/2-1) * dx; % spatial grid

% count the number of nonzero entries, output to command line:
fprintf('There are %i nonzero entries in w\n', sum(w(:)));

% show original window:
f1 = figure(1); clf;
imagesc(x, x, w);
axis('image', 'xy');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
exportgraphics(f1, 'mask2D.png');

% compute autocorrelation from original window using Matlab's xcorr2:
c1 = xcorr2(w, w);
xCor2 = (-N+1 : N-1)*dx; % spatial grid for xcorr2 result

% find the max value and its row & column:
[M, I] = max(c1(:));
[row, col] = ind2sub([2*N-1, 2*N-1], I);

% show window correlation from xcorr2:
f2 = figure(2); clf;
imagesc(xCor2, xCor2, c1);
axis('image', 'xy');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Max value is %i, in row %i and col %i', M, row, col));
exportgraphics(f2, 'maskCorr_xcorr2.png');

%% FT-based calculation of the correlation
% define zero-padded window and grid so that non-circular correlation is computed:
NPad = 2*N; 
wPad = zeros(NPad, NPad); 
idxOrig = (-N/2 : N/2-1) + N+1; 
wPad(idxOrig, idxOrig) = w; 
xPad = (-NPad/2 : NPad/2-1)*dx; 

% compute autocorrelation from padded window using Fourier transforms:
df = 1/(NPad*dx); % Corrected df for zero-padded grid
c2 = ift2(abs(ft2(wPad, dx)).^2, df);

% find the max value and its row & column:
[M2, I2] = max(c2(:));
[row2, col2] = ind2sub([2*N, 2*N], I2);

% show window correlation from Fourier transform:
f4 = figure(4); clf;
imagesc(xPad, xPad, abs(c2));
axis('image', 'xy');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('FT Max: %0.2f at (%i, %i)', M2, row2, col2));
exportgraphics(f4, 'maskCorr_FT.png');

%% Combined Figure for corrUnbiased2D.html
f5 = figure(5); clf;
set(f5, 'OuterPosition', [100 100 900 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(x, x, w);
axis image xy; colorbar;
xlabel('x [m]'); ylabel('y [m]');
title('Original Mask (w)');

nexttile;
% Using the FT result for the unbiased map
imagesc(xPad, xPad, abs(c2));
axis image xy; colorbar;
xlabel('x [m]'); ylabel('y [m]');
title('Autocorrelation of Mask (A)');

exportgraphics(f5, 'checkXcorr2_mask.png');