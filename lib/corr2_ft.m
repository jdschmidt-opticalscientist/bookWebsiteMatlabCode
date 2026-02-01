function [c, idx] = corr2_ft_new(u1, u2, mask, delta)
% function [c, idx] = corr2_ft_new(u1, u2, mask, delta)

    N = size(u1, 1);
    c = zeros(N);
    delta_f = 1/(N*delta);  % frequency grid spacing [m]
    
    U1 = ft2(u1 .* mask, delta);    % DFTs of signals
    U2 = ft2(u2 .* mask, delta);
    U12corr = ift2(conj(U1) .* U2, delta_f);

    areamask = sum(sum(mask)) * delta^2;
    maskcorr = ift2(abs(ft2(mask, delta)).^2, delta_f) / areamask;
    idx = maskcorr >= delta^2 / areamask;
    c(idx) =  U12corr(idx) ./ maskcorr(idx);