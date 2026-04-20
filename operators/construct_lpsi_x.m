function L_psi_x = construct_lpsi_x(phi_x_full, params)
% CONSTRUCT_LPSI_X  Build the x-component of the covariant Laplacian for psi.
%
%   L_psi_x = construct_lpsi_x(phi_x_full, params)
%
%   Constructs a sparse matrix representing the x-direction part of the
%   gauge-covariant Laplacian (nabla - i*A)^2 acting on the order parameter
%   psi, discretized with Peierls phases:
%
%       (L_psi_x * psi)_m = exp(-i * phi_x(m-1)) * psi(m-1)
%                          - 2 * psi(m)
%                          + exp(+i * phi_x(m))   * psi(m+1)
%
%   divided by hx^2 (the caller handles the 1/h^2 factor).
%
%   Inputs:
%     phi_x_full : double vector, length n_nodes_full.
%                  x-component of vector potential on the full grid.
%     params     : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     L_psi_x : sparse matrix (n_nodes_full × n_nodes_full).
%               Only rows corresponding to interior nodes are populated.
%
%   See also CONSTRUCT_LPSI_Y, CONSTRUCT_LPSI_Z

    n_full = params.n_nodes_full;
    m      = params.interior_to_full;  % interior node indices in full grid

    % Off-diagonal: left neighbor (i-1) with Peierls phase exp(-i * phi_x(m-1))
    L_psi_x = sparse(m, m - 1, exp(-1i * phi_x_full(m - 1)), n_full, n_full);

    % Diagonal: -2
    L_psi_x = L_psi_x + sparse(m, m, -2, n_full, n_full);

    % Off-diagonal: right neighbor (i+1) with Peierls phase exp(+i * phi_x(m))
    L_psi_x = L_psi_x + sparse(m, m + 1, exp(1i * phi_x_full(m)), n_full, n_full);
end
