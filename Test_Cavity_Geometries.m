lambda = 1.064e-6; % m
L = 500; % m
R = 2*L; % m

% ABCD matrices
FS = @(L) [1, L; 0, 1];
Mmir = @(R) [1, 0; -2/R, 1]; % mirror reflection matrix (paraxial)

M = Mmir(R) * FS(L) * Mmir(R) * FS(L); % one round-trip ABCD
A = M(1,1); B = M(1,2); C = M(2,1); D = M(2,2);

% Solve for q: C*q^2 + (D-A)*q - B = 0
coeffs = [C, (D-A), -B];
q_roots = roots(coeffs);

% Pick root with positive imaginary part
q = q_roots(find(imag(q_roots) > 0, 1));

% Derive waist and spot at mirror
invq = 1/q;
w0 = sqrt(-lambda/(pi * imag(invq))); % waist at cavity center
zR = pi*w0^2 / lambda;
w_mirror = w0 * sqrt(1 + ( (L/2)/zR )^2);

fprintf('For L=%.3g m, R=%.3g m, lambda=%.3g m:\n', L, R, lambda);
fprintf('  w0 = %.4g m = %.3f mm\n', w0, w0*1e3);
fprintf('  zR = %.4g m\n', zR);
fprintf('  w_mirror (radius) = %.4g m = %.3f mm\n', w_mirror, w_mirror*1e3);
