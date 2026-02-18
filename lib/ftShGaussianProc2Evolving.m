function phz = ftShGaussianProc2Evolving(coeffs, N, delta, t, v)
% FTSHGAUSSIANPROC2EVOLVING Synthesize a 2D Gaussian random process using 
% the FT and subharmonic method and shifts samples laterally based on
% velocity v and time t
%
% Usage: phz = ftShGaussianProc2Evolving(coeffs, N, delta, t, v)
%
% Inputs:
%   coeffs - Struct with coefficients and frequencies
%   N      - Number of samples per side (square grid)
%   delta  - Spatial sampling interval [m]
%   PSDFcn - Handle to a function that accepts two frequency vectors [1/m]
%            and returns the 2D Power Spectral Density
%   t      - time [s]
%   v      - velocity vector [m/s]
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.
% Code from the book (e.g., ift2.m) copyright SPIE.

% Apply the shift
dx = v(1)*t;
dy = v(2)*t;

shifted_cn = coeffs.cn ...
    .* exp(-2i * pi * (coeffs.fx * dx + coeffs.fy * dy));

% Synthesize the FFT portion:
% The factor of N^2 or similar is handled inside ift2 if using book definitions
phz_hi = real(ift2(shifted_cn, 1));

% Low-frequency component using subharmonics
x = (-N/2 : N/2-1) * delta;
[xx, yy] = meshgrid(x);
% --- Inside ftShGaussianProc2Evolving ---
phz_lo = zeros(N);
NP = 3;

for p = 1:NP
    % Access the specific level {p} from the cell arrays
    % Apply the shift to this specific level's frequencies
    shifted_sh_cn = coeffs.sh_cn{p} ...
        .* exp(-2i * pi * (coeffs.sh_fx{p} * dx + coeffs.sh_fy{p} * dy));
    
    SH = zeros(N);
    for ii = 1:3
        for jj = 1:3
            % Skip the center bin except for the very last level
            if (ii == 2 && jj == 2) && (p < NP)
                continue;
            end
            
            % Use the correct level's frequencies for the spatial reconstruction
            SH = SH + shifted_sh_cn(ii, jj) ...
                * exp(2i * pi * (coeffs.sh_fx{p}(ii, jj) * xx ...
                + coeffs.sh_fy{p}(ii, jj) * yy));
        end
    end
    phz_lo = phz_lo + SH;
end

% Final output: ensure real and combine ONLY ONCE after the loop
phz = real(phz_lo) + phz_hi;
end