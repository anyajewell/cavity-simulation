function err = Compare_Field_Error(E_numeric, E_analytic, intensity_threshold)
% Compare analytic and numeric complex fields
% Removes best-fit complex scalar alpha before comparing.
%
% intensity_threshold:
%   optional relative threshold for phase comparison mask
%   default = 1e-3

    if nargin < 3
        intensity_threshold = 1e-3;
    end

    % Flatten
    E_num = E_numeric(:);
    E_an  = E_analytic(:);

    % Best-fit complex scalar: E_num ~ alpha * E_an
    alpha = (E_an' * E_num) / (E_an' * E_an);

    % Aligned analytic field
    E_an_fit = alpha * E_an;

    % Complex field error
    diffE = E_num - E_an_fit;
    err.field_L2  = norm(diffE) / norm(E_num);
    err.field_max = max(abs(diffE));

    % Intensity error
    I_num = abs(E_num).^2;
    I_fit = abs(E_an_fit).^2;

    err.intensity_L2  = norm(I_num - I_fit) / norm(I_num);
    err.intensity_max = max(abs(I_num - I_fit));

    % Phase error only where beam actually has intensity
    I_mask = I_num / max(I_num);
    mask = I_mask > intensity_threshold;

    phase_num = angle(E_num(mask));
    phase_fit = angle(E_an_fit(mask));

    phase_diff = angle(exp(1i*(phase_num - phase_fit)));

    err.phase_rms = sqrt(mean(phase_diff.^2));
    err.phase_max = max(abs(phase_diff));

    % Report fit scalar
    err.alpha = alpha;
    err.alpha_mag = abs(alpha);
    err.alpha_phase = angle(alpha);

    % Optional correlation metric
    err.overlap = abs(E_an' * E_num) / (norm(E_an) * norm(E_num));

    fprintf('------------------------------------\n');
    fprintf('Field comparison error metrics\n');
    fprintf('------------------------------------\n');
    fprintf('Best-fit |alpha|        : %.6e\n', err.alpha_mag);
    fprintf('Best-fit phase(alpha)   : %.6e rad\n', err.alpha_phase);
    fprintf('Field L2 error          : %.6e\n', err.field_L2);
    fprintf('Field max error         : %.6e\n', err.field_max);
    fprintf('Intensity L2 error      : %.6e\n', err.intensity_L2);
    fprintf('Intensity max error     : %.6e\n', err.intensity_max);
    fprintf('Phase RMS error (rad)   : %.6e\n', err.phase_rms);
    fprintf('Phase max error (rad)   : %.6e\n', err.phase_max);
    fprintf('Normalized overlap      : %.6e\n', err.overlap);
    fprintf('------------------------------------\n');
end