function F_psi = construct_fpsi(psi_full, params)
% CONSTRUCT_FPSI  Build the nonlinear forcing term for the psi equation.
%
%   F_psi = construct_fpsi(psi_full, params)
%
%   Computes the Ginzburg-Landau nonlinear term:
%
%       F_psi(m) = (1 - |psi(m)|^2) * psi(m)
%
%   at each interior node, returned as a sparse column vector indexed
%   by the interior numbering (1 : n_interior).
%
%   Inputs:
%     psi_full : complex vector, length n_nodes_full.
%                Order parameter on the full grid.
%     params   : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     F_psi : sparse vector (n_interior × 1).
%
%   See also CONSTRUCT_FPHI_X, CONSTRUCT_FPHI_Y, CONSTRUCT_FPHI_Z

    m_full    = params.interior_to_full;
    m_compact = params.interior_numbering;

    psi_interior = psi_full(m_full);
    nonlinear    = (1 - conj(psi_interior) .* psi_interior) .* psi_interior;

    F_psi = sparse(m_compact, 1, nonlinear);
end
