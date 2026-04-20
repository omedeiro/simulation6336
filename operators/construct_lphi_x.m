function L_phi_x = construct_lphi_x(params)
% CONSTRUCT_LPHI_X  Build the Laplacian for phi_x (vector potential, x-component).
%
%   L_phi_x = construct_lphi_x(params)
%
%   Constructs a sparse matrix for the diffusion terms in the phi_x equation.
%   phi_x has diffusion in the y and z directions (cross-directions), so:
%
%       (L_phi_x * phi_x)_m = kappa^2/hy^2 * [phi_x(m+stride_j) - 2*phi_x(m) + phi_x(m-stride_j)]
%                            + kappa^2/hz^2 * [phi_x(m+stride_k) - 2*phi_x(m) + phi_x(m-stride_k)]
%
%   Inputs:
%     params : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     L_phi_x : sparse matrix (n_nodes_full × n_nodes_full).
%
%   See also CONSTRUCT_LPHI_Y, CONSTRUCT_LPHI_Z

    n_full   = params.n_nodes_full;
    stride_j = params.stride_j;
    stride_k = params.stride_k;
    kappa    = params.kappa;
    hy       = params.hy;
    hz       = params.hz;
    m        = params.interior_to_full;

    coeff_y = kappa^2 / hy^2;
    coeff_z = kappa^2 / hz^2;

    % Diagonal
    L_phi_x = sparse(m, m, -2 * (coeff_y + coeff_z), n_full, n_full);

    % y-direction neighbors
    if params.Ny > 1
        L_phi_x = L_phi_x + sparse(m, m + stride_j, coeff_y, n_full, n_full);
        L_phi_x = L_phi_x + sparse(m, m - stride_j, coeff_y, n_full, n_full);
    end

    % z-direction neighbors
    if params.Nz > 1
        L_phi_x = L_phi_x + sparse(m, m + stride_k, coeff_z, n_full, n_full);
        L_phi_x = L_phi_x + sparse(m, m - stride_k, coeff_z, n_full, n_full);
    end
end
