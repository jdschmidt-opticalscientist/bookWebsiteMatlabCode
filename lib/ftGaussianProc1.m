function g = ftGaussianProc1(N, delta, PSDFcn)
% function g = ftGaussianProc1(N, delta, PSD)
% PSD must have N elements
    df = 1/(N*delta);   % frequency grid spacing [1/m]
    f = (-N/2 : N/2-1) * df;

    PSD = PSDFcn(f);
    PSD = PSD(:); % ensure that PSD is a column
    del_f = 1/(N*delta);   % frequency grid spacing [1/m]
    % random draws of Fourier coefficients:
    cn = complex(randn(N, 1), randn(N, 1)) .* sqrt(PSD*del_f);
    % synthesize the process:
    g = real(ift(cn, 1));