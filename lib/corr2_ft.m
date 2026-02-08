function [c, U12corr] = corr2_ft(u1_pad, u2_pad, mask_pad, delta)
% CORR2_FT Compute unbiased 2D cross-correlation using FTs
%
% This function assumes inputs are already zero-padded (e.g., to 2N x 2N)
% to ensure linear correlation rather than circular correlation.
%
% Inputs:
%   u1_pad, u2_pad - Masked and zero-padded signals
%   mask_pad       - Zero-padded binary mask
%   delta          - Spatial sampling interval [m]
%
% Outputs:
%   c       - Unbiased 2D correlation (Autocovariance)
%   U12corr - Raw (biased) correlation [Units: Signal^2 * Area]

    Npad = size(u1_pad, 1);
    delta_f = 1/(Npad * delta);  % frequency grid spacing [1/m]
    
    % 1. Compute the raw correlation (Numerator / Biased result)
    U1 = ft2(u1_pad, delta);    
    U2 = ft2(u2_pad, delta);
    U12corr = ift2(conj(U1) .* U2, delta_f);
    
    % 2. Compute the Area of Overlap (Denominator)
    M = ft2(mask_pad, delta);
    maskcorr = ift2(abs(M).^2, delta_f);
    
    % 3. Define the valid indices
    % Use the area of one pixel as the threshold to avoid noise at edges
    idx = real(maskcorr) >= (delta^2);
    
    % 4. Unbias
    c = zeros(Npad);
    % Taking the real part removes numerical imaginary residuals
    c(idx) = real(U12corr(idx)) ./ real(maskcorr(idx));
    
    % Ensure biased output is also real-valued for the ensemble average
    U12corr = real(U12corr);
end