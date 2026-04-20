function params = construct_grid_indices(params)
% CONSTRUCT_GRID_INDICES  Build all index arrays for the 3D TDGL grid.
%
%   params = construct_grid_indices(params)
%
%   Populates the parameter structure with index arrays that map between
%   the compact interior numbering and the full (Nx+1)×(Ny+1)×(Nz+1) grid.
%
%   Index Naming Conventions
%   ------------------------
%   For each axis (x, y, z) the following arrays are created:
%
%     Face indices ("all" — every node on that face):
%       x_face_lo, x_face_lo_first, x_face_hi_last, x_face_hi
%       (similarly for y and z)
%
%     Face indices ("inner" — excluding edges/corners):
%       x_face_lo_inner, x_first_inner, x_last_inner, x_face_hi_inner
%       (similarly for y and z)
%
%     Edge indices (for checkerboard boundary-error monitoring):
%       x_edge_face, x_edge_first_inner, x_edge_last_inner
%       (similarly for y and z)
%
%   Grid-wide arrays:
%     interior_to_full : int array (1×n_interior)
%         Maps interior node number → full-grid linear index.
%     bfield_interior  : int array
%         Interior nodes that are one layer further inward (for B = curl A).
%
%   The full-grid linear index for node (i, j, k) is:
%       m = i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*(Ny+1)
%   where i, j, k are 1-based.
%
%   See also SETUP_PARAMETERS, LINEAR_INDEX

    Nx = params.Nx;
    Ny = params.Ny;
    Nz = params.Nz;
    stride_j = Nx + 1;
    stride_k = (Nx + 1) * (Ny + 1);

    % Determine z-loop limit: for quasi-2D (Nz==1), only one z-plane
    if Nz == 1
        k_max = 1;
    else
        k_max = Nz + 1;
    end

    % --- Pre-allocate with dynamic lists (converted to arrays at end) ---
    % Interior nodes
    interior_list       = [];

    % X-face arrays
    x_face_lo_all       = [];  % i=1 face, all (j,k)
    x_first_all         = [];  % i=2 (first interior), all (j,k)
    x_last_all          = [];  % i=Nx (last interior), all (j,k)
    x_face_hi_all       = [];  % i=Nx+1 face, all (j,k)

    x_face_lo_inner     = [];  % i=1, inner (j,k) only
    x_first_inner       = [];  % i=2, inner (j,k)
    x_last_inner        = [];  % i=Nx, inner (j,k)
    x_face_hi_inner     = [];  % i=Nx+1, inner (j,k)

    x_edge_face         = [];  % i=1,Nx+1 edge nodes
    x_edge_first_1      = [];  % i=2 edge-row 1
    x_edge_last_1       = [];  % i=Nx edge-row 1
    x_edge_first_2      = [];  % i=2 edge-row 2
    x_edge_last_2       = [];  % i=Nx edge-row 2

    % Checkerboard arrays for x
    x_first_cb1         = [];
    x_last_cb1          = [];
    x_first_cb2         = [];
    x_last_cb2          = [];

    % Y-face arrays
    y_face_lo_all       = [];
    y_first_all         = [];
    y_last_all          = [];
    y_face_hi_all       = [];

    y_face_lo_inner     = [];
    y_first_inner       = [];
    y_last_inner        = [];
    y_face_hi_inner     = [];

    y_edge_face         = [];
    y_edge_first_1      = [];
    y_edge_last_1       = [];
    y_edge_first_2      = [];
    y_edge_last_2       = [];

    y_first_cb1         = [];
    y_last_cb1          = [];
    y_first_cb2         = [];
    y_last_cb2          = [];

    % Z-face arrays
    z_face_lo_all       = [];
    z_first_all         = [];
    z_last_all          = [];
    z_face_hi_all       = [];

    z_face_lo_inner     = [];
    z_first_inner       = [];
    z_last_inner        = [];
    z_face_hi_inner     = [];

    z_edge_face         = [];
    z_edge_first_1      = [];
    z_edge_last_1       = [];
    z_edge_first_2      = [];
    z_edge_last_2       = [];

    z_first_cb1         = [];
    z_last_cb1          = [];
    z_first_cb2         = [];
    z_last_cb2          = [];

    % === Main triple loop ===
    for k = 1:k_max
        for j = 1:(Ny + 1)
            for i = 1:(Nx + 1)
                m = i + (j - 1) * stride_j + (k - 1) * stride_k;

                % --- Interior nodes ---
                is_interior = (i >= 2) && (i <= Nx) && (j >= 2) && (j <= Ny);
                if Nz > 1
                    is_interior = is_interior && (k >= 2) && (k <= Nz);
                end
                if is_interior
                    interior_list(end+1) = m; %#ok<AGROW>
                end

                % --- X-direction faces (i == 1) ---
                if i == 1
                    x_face_lo_all(end+1)  = m; %#ok<AGROW>
                    m_first               = m + 1;          % i=2
                    m_last                = m + Nx - 1;     % i=Nx
                    m_hi                  = m + Nx;         % i=Nx+1

                    x_first_all(end+1)    = m_first; %#ok<AGROW>
                    x_last_all(end+1)     = m_last; %#ok<AGROW>
                    x_face_hi_all(end+1)  = m_hi; %#ok<AGROW>

                    % Inner condition: not on y or z boundary
                    if Nz > 1
                        is_inner_x = (k >= 2) && (k <= Nz) && (j >= 2) && (j <= Nx);
                    else
                        is_inner_x = (j >= 2) && (j <= Nx);
                    end

                    if is_inner_x
                        x_face_lo_inner(end+1)  = m; %#ok<AGROW>
                        x_first_inner(end+1)    = m_first; %#ok<AGROW>
                        x_last_inner(end+1)     = m_last; %#ok<AGROW>
                        x_face_hi_inner(end+1)  = m_hi; %#ok<AGROW>

                        if mod(j, 2) == 0
                            x_first_cb1(end+1) = m_first; %#ok<AGROW>
                            x_last_cb1(end+1)  = m_last; %#ok<AGROW>
                        else
                            x_first_cb2(end+1) = m_first; %#ok<AGROW>
                            x_last_cb2(end+1)  = m_last; %#ok<AGROW>
                        end
                    end

                    % Edge nodes
                    if (j == 1) || (j == Ny + 1)
                        x_edge_face(end+1) = m; %#ok<AGROW>
                    end

                    if Nz > 1
                        cond_e1 = (j == 2) && (k >= 2) && (k <= Nz);
                    else
                        cond_e1 = (j == 2);
                    end
                    if cond_e1
                        x_edge_first_1(end+1) = m_first; %#ok<AGROW>
                        x_edge_last_1(end+1)  = m_last; %#ok<AGROW>
                    end

                    if Nz > 1
                        cond_e2 = (j == Ny) && (k >= 2) && (k <= Nz);
                    else
                        cond_e2 = (j == Ny);
                    end
                    if cond_e2
                        x_edge_first_2(end+1) = m_first; %#ok<AGROW>
                        x_edge_last_2(end+1)  = m_last; %#ok<AGROW>
                    end
                end

                % --- Y-direction faces (j == 1) ---
                if j == 1
                    y_face_lo_all(end+1) = m; %#ok<AGROW>
                    m_first              = m + stride_j;
                    m_last               = m + stride_j * (Ny - 1);
                    m_hi                 = m + stride_j * Ny;

                    y_first_all(end+1)   = m_first; %#ok<AGROW>
                    y_last_all(end+1)    = m_last; %#ok<AGROW>
                    y_face_hi_all(end+1) = m_hi; %#ok<AGROW>

                    if Nz > 1
                        is_inner_y = (k >= 2) && (k <= Nz) && (i >= 2) && (i <= Ny);
                    else
                        is_inner_y = (i >= 2) && (i <= Ny);
                    end

                    if is_inner_y
                        y_face_lo_inner(end+1)  = m; %#ok<AGROW>
                        y_first_inner(end+1)    = m_first; %#ok<AGROW>
                        y_last_inner(end+1)     = m_last; %#ok<AGROW>
                        y_face_hi_inner(end+1)  = m_hi; %#ok<AGROW>

                        if mod(k, 2) == 0
                            y_first_cb1(end+1) = m_first; %#ok<AGROW>
                            y_last_cb1(end+1)  = m_last; %#ok<AGROW>
                        else
                            y_first_cb2(end+1) = m_first; %#ok<AGROW>
                            y_last_cb2(end+1)  = m_last; %#ok<AGROW>
                        end
                    end

                    if (k == 1) || (k == Nz + 1)
                        y_edge_face(end+1) = m; %#ok<AGROW>
                    end

                    if (k == 2) && (i >= 2) && (i <= Ny)
                        y_edge_first_1(end+1) = m_first; %#ok<AGROW>
                        y_edge_last_1(end+1)  = m_last; %#ok<AGROW>
                    end
                    if (k == Nz) && (i >= 2) && (i <= Ny)
                        y_edge_first_2(end+1) = m_first; %#ok<AGROW>
                        y_edge_last_2(end+1)  = m_last; %#ok<AGROW>
                    end
                end

                % --- Z-direction faces (k == 1) ---
                if k == 1
                    z_face_lo_all(end+1) = m; %#ok<AGROW>
                    m_first              = m + stride_k;
                    m_last               = m + stride_k * (Nz - 1);
                    m_hi                 = m + stride_k * Nz;

                    z_first_all(end+1)   = m_first; %#ok<AGROW>
                    z_last_all(end+1)    = m_last; %#ok<AGROW>
                    z_face_hi_all(end+1) = m_hi; %#ok<AGROW>

                    is_inner_z = (j >= 2) && (j <= Ny) && (i >= 2) && (i <= Nx);
                    if is_inner_z
                        z_face_lo_inner(end+1)  = m; %#ok<AGROW>
                        z_first_inner(end+1)    = m_first; %#ok<AGROW>
                        z_last_inner(end+1)     = m_last; %#ok<AGROW>
                        z_face_hi_inner(end+1)  = m_hi; %#ok<AGROW>

                        if mod(i, 2) == 0
                            z_first_cb1(end+1) = m_first; %#ok<AGROW>
                            z_last_cb1(end+1)  = m_last; %#ok<AGROW>
                        else
                            z_first_cb2(end+1) = m_first; %#ok<AGROW>
                            z_last_cb2(end+1)  = m_last; %#ok<AGROW>
                        end
                    end

                    if (i == 1) || (i == Nx + 1)
                        z_edge_face(end+1) = m; %#ok<AGROW>
                    end
                    if (i == 2) && (j >= 2) && (j <= Ny)
                        z_edge_first_1(end+1) = m_first; %#ok<AGROW>
                        z_edge_last_1(end+1)  = m_last; %#ok<AGROW>
                    end
                    if (i == Nx) && (j >= 2) && (j <= Ny)
                        z_edge_first_2(end+1) = m_first; %#ok<AGROW>
                        z_edge_last_2(end+1)  = m_last; %#ok<AGROW>
                    end
                end

            end % i
        end % j
    end % k

    % === B-field interior (one more layer inward) ===
    bfield_interior = [];
    for k = 1:max(Nz - 1, 1)
        for j = 1:(Ny - 1)
            for i = 1:(Nx - 1)
                if Nz > 1
                    is_bfield = (k <= Nz - 2) && (j <= Ny - 2) && (i <= Nx - 2);
                else
                    is_bfield = (j <= Ny - 2) && (i <= Nx - 2);
                end
                if is_bfield
                    bfield_interior(end+1) = i + (Nx - 1) * (j - 1) + (Nx - 1) * (Ny - 1) * (k - 1); %#ok<AGROW>
                end
            end
        end
    end

    % === Store all index arrays ===
    params.interior_to_full = interior_list;

    % X-face indices
    params.x_face_lo        = x_face_lo_all;
    params.x_first          = x_first_all;
    params.x_last           = x_last_all;
    params.x_face_hi        = x_face_hi_all;
    params.x_face_lo_inner  = x_face_lo_inner;
    params.x_first_inner    = x_first_inner;
    params.x_last_inner     = x_last_inner;
    params.x_face_hi_inner  = x_face_hi_inner;
    params.x_edge_face      = x_edge_face;
    params.x_edge_first_1   = x_edge_first_1;
    params.x_edge_last_1    = x_edge_last_1;
    params.x_edge_first_2   = x_edge_first_2;
    params.x_edge_last_2    = x_edge_last_2;
    params.x_first_cb1      = x_first_cb1;
    params.x_last_cb1       = x_last_cb1;
    params.x_first_cb2      = x_first_cb2;
    params.x_last_cb2       = x_last_cb2;

    % Y-face indices
    params.y_face_lo        = y_face_lo_all;
    params.y_first          = y_first_all;
    params.y_last           = y_last_all;
    params.y_face_hi        = y_face_hi_all;
    params.y_face_lo_inner  = y_face_lo_inner;
    params.y_first_inner    = y_first_inner;
    params.y_last_inner     = y_last_inner;
    params.y_face_hi_inner  = y_face_hi_inner;
    params.y_edge_face      = y_edge_face;
    params.y_edge_first_1   = y_edge_first_1;
    params.y_edge_last_1    = y_edge_last_1;
    params.y_edge_first_2   = y_edge_first_2;
    params.y_edge_last_2    = y_edge_last_2;
    params.y_first_cb1      = y_first_cb1;
    params.y_last_cb1       = y_last_cb1;
    params.y_first_cb2      = y_first_cb2;
    params.y_last_cb2       = y_last_cb2;

    % Z-face indices
    params.z_face_lo        = z_face_lo_all;
    params.z_first          = z_first_all;
    params.z_last           = z_last_all;
    params.z_face_hi        = z_face_hi_all;
    params.z_face_lo_inner  = z_face_lo_inner;
    params.z_first_inner    = z_first_inner;
    params.z_last_inner     = z_last_inner;
    params.z_face_hi_inner  = z_face_hi_inner;
    params.z_edge_face      = z_edge_face;
    params.z_edge_first_1   = z_edge_first_1;
    params.z_edge_last_1    = z_edge_last_1;
    params.z_edge_first_2   = z_edge_first_2;
    params.z_edge_last_2    = z_edge_last_2;
    params.z_first_cb1      = z_first_cb1;
    params.z_last_cb1       = z_last_cb1;
    params.z_first_cb2      = z_first_cb2;
    params.z_last_cb2       = z_last_cb2;

    % B-field interior
    params.bfield_interior  = bfield_interior;

    % Interior numbering (compact 1:n_interior)
    params.interior_numbering = 1:params.n_interior;
end
