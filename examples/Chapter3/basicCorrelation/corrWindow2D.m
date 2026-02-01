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
sprintf('There are %i nonzero entries in w', sum(w(:)))

% show original window:
figure(1); clf;
imagesc(x, x, w);
axis('image', 'xy');
colorbar;
xlabel('x [m]');
ylabel('y [m]');

% compute autocorrelation from original window using Matlab's xcorr2:
c1 = xcorr2(w, w);
% size of c1 is 2N-1 by 2N-1
xCor2 = (-N+1 : N-1)*dx; % spatial grid for xcorr2 result
% find the max value and its row & column:
[M, I] = max(c1(:));
[row, col] = ind2sub([2*N-1, 2*N-1], I);

% show window correlation from xcorr2:
figure(2); clf;
imagesc(xCor2, xCor2, c1);
axis('image', 'xy');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Max value is %i, in row %i and col %i', M, row, col));

%% FT-based calculation of the correlation

% define zero-padded window and grid so that non-circular correlation is
% computed:
NPad = 2*N; % number of grid points in zero-padded array
wPad = zeros(NPad, NPad); % making zero-padded array
idxOrig = (-N/2 : N/2-1) + N+1; % indices for original window
wPad(idxOrig, idxOrig) = w; % fill center of zero-padded array
xPad = (-NPad/2 : NPad/2-1)*dx; % coordinates for zero-padded array

% show padded window:
figure(3); clf;
imagesc(xPad, xPad, wPad);
axis('image', 'xy');
colorbar;
xlabel('x [m]');
ylabel('y [m]');

% compute autocorrelation from padded window using Fourier transforms:
df = 1/NPad; % spatial frequency grid spacing
c2 = ift2(abs(ft2(wPad, dx)).^2, df);
% wFT = ft2(wPad, dx);
% wFtFlip = ft2(flip(flip(wFT,1),2), dx);
% c2 = ift2(wFT.*wFtFlip, df);
% find the max value and its row & column:
[M, I] = max(c2(:));
[row, col] = ind2sub([2*N, 2*N], I);

% show window correlation from Fourier transform:
figure(4); clf;
imagesc(xPad, xPad, abs(c2));
axis('image', 'xy');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Max value is %i, in row %i and col %i', M, row, col));