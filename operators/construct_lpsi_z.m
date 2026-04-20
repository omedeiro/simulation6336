function L_psi_z = construct_lpsi_z(phi_z_full, params)
% CONSTRUCT_LPSI_Z  Build the z-component of the covariant Laplacian for psi.
%
%   L_psi_z = construct_lpsi_z(phi_z_full, params)
%
%   Constructs a sparse matrix representing the z-direction part of the
%   gauge-covariant Laplacian for psi using Peierls phases:
%
%       (L_psi_z * psi)_m = exp(-i * phi_z(m-stride_k)) * psi(m-stride_k)
%                          - 2 * psi(m)
%                          + exp(+i * phi_z(m))          * psi(m+stride_k)
%
%   where stride_k = (Nx+1)*(Ny+1).
%
%   Returns a zero-contribution matrix when Nz == 1 (quasi-2D case).
%
%   Inputs:
%     phi_z_full : double vector, length n_nodes_full.
%                  z-component of vector potential on the full grid.
%     params     : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     L_psi_z : sparse matrix (n_nodes_full × n_nodes_full).
%
%   See also CONSTRUCT_LPSI_X, CONSTRUCT_LPSI_Y

    n_full   = params.n_nodes_full;
    stride_k = params.stride_k;
    m        = params.interior_to_full;

    % Diagonal: -2
    L_psi_z = sparse(m, m, -2, n_full, n_full);

    % Off-diagonal terms only if Nz > 1
    if params.Nz > 1
        % Forward neighbor (k+1)
        L_psi_z = L_psi_z + sparse(m, m + stride_k, ...
                    exp(1i * phi_z_full(m)), n_full, n_full);
        % Backward neighbor (k-1)
        L_psi_z = L_psi_z + sparse(m, m - stride_k, ...
                    exp(-1i * phi_z_full(m - stride_k)), n_full, n_full);
    end
end
