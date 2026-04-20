function L_psi_y = construct_lpsi_y(phi_y_full, params)
% CONSTRUCT_LPSI_Y  Build the y-component of the covariant Laplacian for psi.
%
%   L_psi_y = construct_lpsi_y(phi_y_full, params)
%
%   Constructs a sparse matrix representing the y-direction part of the
%   gauge-covariant Laplacian for psi using Peierls phases:
%
%       (L_psi_y * psi)_m = exp(-i * phi_y(m-stride_j)) * psi(m-stride_j)
%                          - 2 * psi(m)
%                          + exp(+i * phi_y(m))          * psi(m+stride_j)
%
%   where stride_j = Nx + 1.
%
%   Inputs:
%     phi_y_full : double vector, length n_nodes_full.
%                  y-component of vector potential on the full grid.
%     params     : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     L_psi_y : sparse matrix (n_nodes_full × n_nodes_full).
%
%   See also CONSTRUCT_LPSI_X, CONSTRUCT_LPSI_Z

    n_full   = params.n_nodes_full;
    stride_j = params.stride_j;
    m        = params.interior_to_full;

    % Diagonal: -2
    L_psi_y = sparse(m, m, -2, n_full, n_full);

    % Off-diagonal terms only if Ny > 1
    if params.Ny > 1
        % Forward neighbor (j+1)
        L_psi_y = L_psi_y + sparse(m, m + stride_j, ...
                    exp(1i * phi_y_full(m)), n_full, n_full);
        % Backward neighbor (j-1)
        L_psi_y = L_psi_y + sparse(m, m - stride_j, ...
                    exp(-1i * phi_y_full(m - stride_j)), n_full, n_full);
    end
end
