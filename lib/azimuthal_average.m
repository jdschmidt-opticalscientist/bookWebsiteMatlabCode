function [rvals, dat_avg] = azimuthal_average(x, y, data)
% azimuthal_average Averaging of 2-D data over the azimuthal angle
%   [rvals, dat_avg] = azimuthal_average(x, y, data) calculates the
%   mean value of 'data' for all points sharing the same radial distance
%   from the origin, effectively providing a radial profile.
%
%   Inputs:
%       x, y    - 2-D grids of Cartesian coordinates (from meshgrid)
%       data    - 2-D array of values to be averaged
%
%   Outputs:
%       rvals   - Unique radial distances from the origin
%       dat_avg - Average value of 'data' at each corresponding rval
%
% Supplemental code for "Numerical Simulation of Optical Wave 
% Propagation with Examples in MATLAB"
%
% Copyright (c) 2026, Jason D. Schmidt. Licensed under BSD 3-Clause.

    % Convert Cartesian coordinates to polar
    [~, r] = cart2pol(x, y);
    
    % Find unique radial distances and their mapping indices
    % ic is the index mapping each pixel in r(:) to a unique value in rvals
    [rvals, ~, ic] = unique(r(:));
    
    % Use accumarray to perform the mean operation on all data values 
    % that map to the same unique radial index
    dat_avg = accumarray(ic, data(:), size(rvals), @mean);
end