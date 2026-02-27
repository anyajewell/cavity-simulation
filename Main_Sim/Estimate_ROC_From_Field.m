function [Reff, a_fit, stats] = Estimate_ROC_From_Field(E, X, Y, k0, r_fit_max, Irel_thresh)
% Estimate effective radius of curvature at a plane by fitting phase ~ a*r^2 + b
% Reff = k0/a_fit (sign included)
%
% Inputs:
%   E           complex field at plane
%   X,Y         grids (same size as E) [m]
%   k0          2*pi/lambda (or n*k0 if in medium)
%   r_fit_max   max radius to fit over (e.g., 0.2*D) [m]
%   Irel_thresh relative intensity threshold for fit points (e.g., 0.2)
%
% Outputs:
%   Reff        effective ROC [m] (sign included)
%   a_fit       fitted quadratic phase coefficient [rad/m^2]
%   stats       struct with goodness-of-fit info

    r2 = X.^2 + Y.^2;
    I  = abs(E).^2;

    % Mask: near-axis region + sufficiently bright region
    Imax = max(I, [], 'all');
    mask = (r2 <= r_fit_max^2) & (I >= Irel_thresh*Imax);

    idx = find(mask);
    if numel(idx) < 50
        Reff = NaN; a_fit = NaN;
        stats = struct('rmse',NaN,'N',numel(idx));
        return
    end

    r2v = r2(idx);
    Ev  = E(idx);

    % Reference to reduce arbitrary piston phase
    [~, i0] = max(I(idx));               % pick brightest point in fit region
    Eref = Ev(i0);
    ph = angle(Ev .* conj(Eref));        % phase relative to reference
    % Unwrap in order of increasing r^2 (radial unwrap)
    [r2s, ord] = sort(r2v);
    phs = unwrap(ph(ord));

    % Weighted least squares fit: phs ≈ a*r2s + b
    % Use intensity weights to de-emphasize noisy low-I points
    ws = I(idx);
    ws = ws(ord);
    ws = ws / max(ws);

    A = [r2s, ones(size(r2s))];
    W = diag(ws);
    p = (A.'*W*A)\(A.'*W*phs);   % [a; b]
    a_fit = p(1);

    % Convert to effective ROC
    % If phase = a*r^2, then a = k0/Reff  (with sign)
    Reff = k0 / (2*a_fit);

    % Fit quality
    phhat = A*p;
    err = phs - phhat;
    stats = struct();
    stats.rmse = sqrt(mean(err.^2));
    stats.N = numel(r2s);
end