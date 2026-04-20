function F_phi_x = construct_fphi_x(psi_full, phi_x_full, phi_y_full, phi_z_full, params)
% CONSTRUCT_FPHI_X  Build the nonlinear forcing term for the phi_x equation.
%
%   F_phi_x = construct_fphi_x(psi_full, phi_x_full, phi_y_full, phi_z_full, params)
%
%   Computes the curl-curl and supercurrent contributions to dphi_x/dt:
%
%       F_phi_x(m) = kappa^2/hy^2 * [curlz contribution from phi_y]
%                  + kappa^2/hz^2 * [curly contribution from phi_z]  (3D only)
%                  + Im[exp(-i*phi_x(m)) * conj(psi(m)) * psi(m+1)]
%
%   The curl-curl terms enforce the discrete curl(curl(A)) identity.
%   The last term is the gauge-invariant supercurrent.
%
%   Inputs:
%     psi_full, phi_x_full, phi_y_full, phi_z_full :
%         Field arrays on the full grid (length n_nodes_full).
%     params : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     F_phi_x : sparse vector (n_interior × 1).
%
%   See also CONSTRUCT_FPHI_Y, CONSTRUCT_FPHI_Z, CONSTRUCT_FPSI

    stride_j = params.stride_j;
    stride_k = params.stride_k;
    kappa    = params.kappa;
    hy       = params.hy;
    hz       = params.hz;
    m_full   = params.interior_to_full;
    m_compact = params.interior_numbering;

    % Supercurrent: Im[exp(-i*phi_x(m)) * psi*(m) * psi(m+1)]
    supercurrent = imag(exp(-1i * phi_x_full(m_full)) .* ...
                        conj(psi_full(m_full)) .* psi_full(m_full + 1));

    % Curl-curl from phi_y (always present)
    curl_y = (kappa^2 / hy^2) * ( ...
        -phi_y_full(m_full + 1) + phi_y_full(m_full) + ...
         phi_y_full(m_full + 1 - stride_j) - phi_y_full(m_full - stride_j));

    if params.Nz > 1
        % Curl-curl from phi_z
        curl_z = (kappa^2 / hz^2) * ( ...
            -phi_z_full(m_full + 1) + phi_z_full(m_full) + ...
             phi_z_full(m_full + 1 - stride_k) - phi_z_full(m_full - stride_k));

        values = curl_y + curl_z + supercurrent;
    else
        values = curl_y + supercurrent;
    end

    F_phi_x = sparse(m_compact, 1, values);
end
