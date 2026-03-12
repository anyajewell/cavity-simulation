function w = BeamWaist_e2(E, X, Y)

% Intensity
I = abs(E).^2;

% Centroid
den = sum(I,'all');
xc = sum(X.*I,'all')/den;
yc = sum(Y.*I,'all')/den;

% Find grid indices nearest centroid
[~,ix] = min(abs(X(1,:) - xc));
[~,iy] = min(abs(Y(:,1) - yc));

% Horizontal intensity slice
Ix = I(iy,:);

% Target intensity level
Imax = max(Ix);
target = Imax*exp(-2);

% Search right side of beam
[~,j] = min(abs(Ix(ix:end) - target));

% Radius
xvals = X(1,:);
w = abs(xvals(ix+j-1) - xc);

end