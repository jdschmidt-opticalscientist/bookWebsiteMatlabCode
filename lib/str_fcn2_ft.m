function [D, maskCorr] = str_fcn2_ft(u_mask_pad, mask_pad, delta)
% STR_FCN2_FT Compute unbiased 2D structure function using FTs
%
% Usage: [D, maskCorr] = str_fcn2_ft(u_mask_pad, mask_pad, delta)
%
% This function assumes inputs are already zero-padded (e.g., to 2N x 2N)
% and the mask is applied to the signal.
%
% Inputs:
%   u_mask_pad - Masked and zero-padded signal
%   mask_pad   - Zero-padded binary mask
%   delta      - Spatial sampling interval [m]
%
% Outputs:
%   D          - Unbiased 2D structure function
%   maskCorr   - The mask autocorrelation (Area of Overlap)

    Npad = size(u_mask_pad, 1);
    df = 1 / (Npad * delta);
    
    % Fourier Transforms
    U = ft2(u_mask_pad, delta);       % F{u * W}
    U2 = ft2(u_mask_pad.^2, delta);   % F{u^2 * W}
    M = ft2(mask_pad, delta);         % F{W}
    
    % 1. Compute the raw (biased) structure function numerator
    % This handles the <u^2(x)> + <u^2(x+tau)> - 2<u(x)u(x+tau)> logic
    % through the correlation theorem.
    term12 = 2 * real(ift2(U2 .* conj(M), df));
    term3  = 2 * real(ift2(U .* conj(U), df));
    rawStructure = term12 - term3;
    
    % 2. Compute the Area of Overlap (Denominator)
    maskCorr = real(ift2(M .* conj(M), df));
    
    % 3. Define valid indices (threshold at one pixel area)
    idx = maskCorr >= (delta^2);
    
    % 4. Unbias
    D = zeros(Npad, Npad);
    D(idx) = rawStructure(idx) ./ maskCorr(idx);
end