function coeffs = ftShGaussianProc2Coeffs(N, delta, PSDFcn)
% FTSHGAUSSIANPROC2Coeffs Synthesize a 2D Gaussian random process 
% with subharmonics
%
% Usage: coeffs = ftShGaussianProc2Coeffs(N, delta, PSDFcn)
%
% This function uses additional frequencies beyond those available in the
% FFT to boost low-frequency content
%
% Inputs:
%   N      - Number of samples per side
%   delta  - Spatial sampling interval [m]
%   PSDFcn - Handle to a function for the 2D PSD
%
% Outputs:
%   coeffs - FFT and subharmonic Fourier-domain coefficients
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

% 1. High-frequency component from FFT method
df = 1/(N*delta);   % frequency grid spacing [1/m]
f = (-N/2 : N/2-1) * df;
[coeffs.fx, coeffs.fy] = meshgrid(f);

% Evaluate PSD on the 2D grid
PSD = PSDFcn(coeffs.fx, coeffs.fy);

% Random draws of Fourier coefficients:
% We multiply by df (sqrt(PSD * df_x * df_y))
coeffs.cn = complex(randn(N), randn(N)) .* sqrt(PSD) * df;
coeffs.cn(N/2+1, N/2+1) = 0;
% 2. Low-frequency component using subharmonics
D = N * delta;

% Number of subharmonic levels
NP = 3; 

% Pre-allocate cell arrays to prevent overwriting
coeffs.sh_cn = cell(1, NP);
coeffs.sh_fx = cell(1, NP);
coeffs.sh_fy = cell(1, NP);

for p = 1:NP
    df_p = 1 / (3^p * D); 
    f_p = [-1, 0, 1] * df_p;
    [sh_fx_p, sh_fy_p] = meshgrid(f_p); % Use local temp variables
    
    % Evaluate PSD
    PSD_p = PSDFcn(sh_fx_p, sh_fy_p);
    
    % Store in cell arrays using the index {p}
    coeffs.sh_fx{p} = sh_fx_p;
    coeffs.sh_fy{p} = sh_fy_p;
    coeffs.sh_cn{p} = complex(randn(3), randn(3)) .* sqrt(PSD_p) * df_p;
end

end