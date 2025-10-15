function amask = Create_Absorbing_Mask(W, X, Y)
    % Set up an absorbing mask
    % --- Absorbing mask B: radial profile smoothed with separable Gaussian convolution ---
    % Tunable params
    taper_width_m = 0.30 * W;   % make this large (0.25-0.40 * W)
    inner_radius = W - taper_width_m;
    outer_radius = W;           % full absorb at W
    poly_power = 2;             % polynomial fall (1 linear, 2 quadratic, etc.)
    
    % Make base mask (polynomial roll)
    ra = sqrt(X.^2 + Y.^2);
    amask_base = ones(size(ra));
    zone = (ra > inner_radius) & (ra < outer_radius);
    s = (ra(zone) - inner_radius) ./ (outer_radius - inner_radius); % 0..1
    amask_base(zone) = (1 - s).^poly_power;
    amask_base(ra >= outer_radius) = 0;
    %amask = 1;
    
    % Gaussian smoothing by 1D separable kernel (cheap & robust)
    % choose sigma in grid points (not meters). More sigma => smoother.
    sigma_m = 0.08 * W;                      % smoothing width in meters (try 0.05-0.12*W)
    sigma_px = max(1, round(sigma_m / dx));  % convert to pixels
    % Create 1D Gaussian kernel
    k_half = ceil(4 * sigma_px);
    xk = -k_half:k_half;
    g1d = exp( - (xk.^2) / (2 * sigma_px^2) );
    g1d = g1d / sum(g1d); % normalize
    
    % separable convolution: first along x, then y
    amask_sm = conv2(g1d, g1d, amask_base, 'same');
    
    % enforce bounds and numerical floor
    amask_sm(amask_sm<1e-12) = 0;
    amask = amask_sm;
    
    % visualize
    figure(101); imagesc(x, y, amask); axis equal tight; colorbar;
    title(sprintf('Smoothed mask (taper %.3fm, sigma_px=%d)', taper_width_m, sigma_px));
    drawnow;
end