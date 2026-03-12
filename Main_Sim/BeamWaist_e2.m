function w = BeamWaist_e2(E, X, Y)

    I = abs(E).^2;

    % Centroid
    den = sum(I,'all');
    xc = sum(X.*I,'all') / den;
    yc = sum(Y.*I,'all') / den;

    % Nearest grid point to centroid
    [~,ix] = min(abs(X(1,:) - xc));
    [~,iy] = min(abs(Y(:,1) - yc));

    % Horizontal slice through centroid
    x = X(1,:);
    Ix = I(iy,:);

    % Normalize
    Ix = Ix / max(Ix);

    % Right-hand side from beam center outward
    xr = x(ix:end) - xc;
    Ir = Ix(ix:end);

    % Find end of central lobe: first place profile starts increasing
    dIr = diff(Ir);
    kRise = find(dIr > 0, 1, 'first');

    if isempty(kRise)
        kEnd = numel(Ir);
    else
        kEnd = kRise;
    end

    % Restrict search to central lobe only
    xr_core = xr(1:kEnd);
    Ir_core = Ir(1:kEnd);

    target = exp(-2);

    % Find first crossing within central lobe
    k = find(Ir_core <= target, 1, 'first');

    if isempty(k) || k == 1
        w = NaN;
        return
    end

    % Linear interpolation
    x1 = xr_core(k-1); x2 = xr_core(k);
    y1 = Ir_core(k-1); y2 = Ir_core(k);

    w = x1 + (target - y1)*(x2 - x1)/(y2 - y1);
    w = abs(w);

end