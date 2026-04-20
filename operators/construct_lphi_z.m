function L_phi_z = construct_lphi_z(params)
% CONSTRUCT_LPHI_Z  Build the Laplacian for phi_z (vector potential, z-component).
%
%   L_phi_z = construct_lphi_z(params)
%
%   Constructs a sparse matrix for the diffusion terms in the phi_z equation.
%   phi_z has diffusion in the x and y directions (cross-directions):
%
%       (L_phi_z * phi_z)_m = kappa^2/hx^2 * [phi_z(m+1) - 2*phi_z(m) + phi_z(m-1)]
%                            + kappa^2/hy^2 * [phi_z(m+stride_j) - 2*phi_z(m) + phi_z(m-stride_j)]
%
%   Inputs:
%     params : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     L_phi_z : sparse matrix (n_nodes_full × n_nodes_full).
%
%   See also CONSTRUCT_LPHI_X, CONSTRUCT_LPHI_Y

    n_full   = params.n_nodes_full;
    stride_j = params.stride_j;
    kappa    = params.kappa;
    hx       = params.hx;
    hy       = params.hy;
    m        = params.interior_to_full;

    coeff_x = kappa^2 / hx^2;
    coeff_y = kappa^2 / hy^2;

    % Diagonal
    L_phi_z = sparse(m, m, -2 * (coeff_x + coeff_y), n_full, n_full);

    % x-direction neighbors
    if params.Nx > 1
        L_phi_z = L_phi_z + sparse(m, m + 1, coeff_x, n_full, n_full);
        L_phi_z = L_phi_z + sparse(m, m - 1, coeff_x, n_full, n_full);
    end

    % y-direction neighbors
    if params.Ny > 1
        L_phi_z = L_phi_z + sparse(m, m + stride_j, coeff_y, n_full, n_full);
        L_phi_z = L_phi_z + sparse(m, m - stride_j, coeff_y, n_full, n_full);
    end
end
