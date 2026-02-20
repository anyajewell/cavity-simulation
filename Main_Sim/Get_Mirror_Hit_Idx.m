function [hit_idx, m1_idx, m2_idx] = Get_Mirror_Hit_Idx(zs, L, tol)
    if nargin < 3 || isempty(tol)
        tol = 1e-9 * max(1, abs(L));
    end
    m1_idx = find(abs(zs - L/2) < tol);
    m2_idx = find(abs(zs + L/2) < tol);
    hit_idx = sort([m1_idx(:); m2_idx(:)]);
end