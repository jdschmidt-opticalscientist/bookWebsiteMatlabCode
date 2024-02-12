function [phz_lo, phz_hi] ...
    = ftShGaussianProc1(M, dt, psdThFcn)
% function [phz_lo, phz_hi] ...
%     = ftShGaussianProc1(M, dt, psdTh)

    % high-frequency screen from FFT method:
    phz_hi = ftGaussianProc1(M, dt, psdThFcn);
    
    % subharmonics:
    x = (-M/2 : M/2-1).' * dt;
    D = M*dt;
    % initialize low-freq screen:
    phz_lo = zeros(M, 1);
    % loop over frequency grids with spacing 1/(3^p*D)
    NP = 5;
    for p = 1:NP
        % setup the PSD:
        del_f = 1 / (3^p*D); % frequency grid spacing [1/m]
        fx = linspace(-0.5, 0.5, 3) * del_f;
        PSD_phi = psdThFcn(fx);
        %PSD_phi(2) = 0;
        % random draws of Fourier coefficients:
        cn = complex(randn(1,3), randn(1,3)) .* sqrt(PSD_phi*del_f);
        SH = zeros(M, 1);
        % loop over frequencies on this grid:
        for ii = 1:3
            SH = SH + cn(ii) * exp(1i*2*pi*fx(ii)*x);
        end
        phz_lo = phz_lo + SH;   % accumulate subharmonics
    end
    phz_lo = real(phz_lo) - mean(real(phz_lo(:)));