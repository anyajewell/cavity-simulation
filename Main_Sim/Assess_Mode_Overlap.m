function [] = Assess_Mode_Overlap(result2D, result3D)

    E2D = result2D;
    E3D = result3D;
    
    assert(isequal(size(E2D), size(E3D)), ...
        'Fields must use the same transverse grid.');
    
    % Normalized complex-field overlap.
    inner_product = sum(conj(E2D).*E3D, 'all');
    norm_product = sqrt(sum(abs(E2D).^2, 'all') ...
        * sum(abs(E3D).^2, 'all'));
    
    rho = abs(inner_product)/norm_product;
    power_overlap = rho^2;
    
    fprintf('Field overlap = %.8f\n', rho);
    fprintf('Power overlap = %.6f %%\n', 100*power_overlap);

end