function L_phi_y = construct_lphi_y(params)
% CONSTRUCT_LPHI_Y  Build the Laplacian for phi_y (vector potential, y-component).
%
%   L_phi_y = construct_lphi_y(params)
%
%   Constructs a sparse matrix for the diffusion terms in the phi_y equation.
%   phi_y has diffusion in the x and z directions (cross-directions):
%
%       (L_phi_y * phi_y)_m = kappa^2/hx^2 * [phi_y(m+1) - 2*phi_y(m) + phi_y(m-1)]
%                            + kappa^2/hz^2 * [phi_y(m+stride_k) - 2*phi_y(m) + phi_y(m-stride_k)]
%
%   Inputs:
%     params : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     L_phi_y : sparse matrix (n_nodes_full × n_nodes_full).
%
%   See also CONSTRUCT_LPHI_X, CONSTRUCT_LPHI_Z

    n_full   = params.n_nodes_full;
    stride_k = params.stride_k;
    kappa    = params.kappa;
    hx       = params.hx;
    hz       = params.hz;
    m        = params.interior_to_full;

    coeff_x = kappa^2 / hx^2;
    coeff_z = kappa^2 / hz^2;

    % Diagonal
    L_phi_y = sparse(m, m, -2 * (coeff_x + coeff_z), n_full, n_full);

    % x-direction neighbors
    if params.Nx > 1
        L_phi_y = L_phi_y + sparse(m, m + 1, coeff_x, n_full, n_full);
        L_phi_y = L_phi_y + sparse(m, m - 1, coeff_x, n_full, n_full);
    end

    % z-direction neighbors
    if params.Nz > 1
        L_phi_y = L_phi_y + sparse(m, m + stride_k, coeff_z, n_full, n_full);
        L_phi_y = L_phi_y + sparse(m, m - stride_k, coeff_z, n_full, n_full);
    end
end
