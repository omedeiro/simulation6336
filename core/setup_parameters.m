function params = setup_parameters(varargin)
% SETUP_PARAMETERS  Create and validate the simulation parameter structure.
%
%   params = setup_parameters('Nx', 10, 'Ny', 10, 'Nz', 3, ...)
%
%   All parameters are specified as Name-Value pairs.  Unspecified
%   parameters take sensible defaults.
%
%   Grid Parameters
%   ---------------
%     Nx, Ny, Nz : int
%         Number of grid cells in each direction.  The full grid has
%         (Nx+1)×(Ny+1)×(Nz+1) nodes (indices 1 … Nx+1 in MATLAB).
%         Set Nz = 1 for a quasi-2D simulation.  (default: 10, 10, 3)
%
%     hx, hy, hz : double
%         Grid spacing in each direction, in units of the coherence
%         length ξ.  (default: 1.0)
%
%   Physics Parameters
%   ------------------
%     kappa : double
%         Ginzburg–Landau parameter κ = λ/ξ.  (default: 5.0)
%
%   Applied Field Magnitudes
%   ------------------------
%     applied_Bx, applied_By, applied_Bz : double
%         Magnitude of the applied magnetic field along each axis,
%         in units of Φ₀/(2πξ²).  (default: 0.0, 0.0, 0.0)
%
%   Boundary Conditions
%   -------------------
%     periodic_x, periodic_y, periodic_z : logical
%         If true, periodic boundary conditions are applied along that
%         axis.  Otherwise zero-current (natural) BCs are used.
%         (default: false)
%
%   Time Integration
%   ----------------
%     t_start, t_stop : double
%         Simulation time window.  (default: 0.0, 100.0)
%
%     dt : double
%         Time step.  For Forward Euler the CFL condition requires
%         dt < h^2 / (4 κ^2).  (default: 0.1)
%
%   Visualization
%   -------------
%     slice_z : int
%         z-plane index (1-based) for 2D slice visualization.
%         (default: 1)
%
%     frames : int
%         Number of frames to save for animation.  (default: 100)
%
%     visualize_save : logical
%         If true, save animated GIF during visualization.  (default: false)
%
%     linearize : int
%         If 1, enable linearization analysis mode.  (default: 0)
%
%   Returns
%   -------
%     params : struct
%         Complete parameter structure ready for construct_grid_indices().
%
%   Example
%   -------
%     params = setup_parameters('Nx', 20, 'Ny', 20, 'Nz', 4, ...
%                               'kappa', 5, 'applied_Bz', 0.6);
%
%   See also CONSTRUCT_GRID_INDICES

    % ---- Parse inputs ----
    ip = inputParser;
    ip.addParameter('Nx',           10,     @(x) isnumeric(x) && x >= 2);
    ip.addParameter('Ny',           10,     @(x) isnumeric(x) && x >= 2);
    ip.addParameter('Nz',            3,     @(x) isnumeric(x) && x >= 1);
    ip.addParameter('hx',          1.0,     @isnumeric);
    ip.addParameter('hy',          1.0,     @isnumeric);
    ip.addParameter('hz',          1.0,     @isnumeric);
    ip.addParameter('kappa',       5.0,     @isnumeric);
    ip.addParameter('applied_Bx',  0.0,     @isnumeric);
    ip.addParameter('applied_By',  0.0,     @isnumeric);
    ip.addParameter('applied_Bz',  0.0,     @isnumeric);
    ip.addParameter('periodic_x', false,    @islogical);
    ip.addParameter('periodic_y', false,    @islogical);
    ip.addParameter('periodic_z', false,    @islogical);
    ip.addParameter('t_start',     0.0,     @isnumeric);
    ip.addParameter('t_stop',    100.0,     @isnumeric);
    ip.addParameter('dt',         0.1,      @isnumeric);
    ip.addParameter('slice_z',          1,      @isnumeric);
    ip.addParameter('frames',          100,     @isnumeric);
    ip.addParameter('visualize_save', false,    @islogical);
    ip.addParameter('linearize',        0,      @isnumeric);
    ip.parse(varargin{:});
    params = ip.Results;

    % ---- Derived quantities ----
    params.n_nodes_full     = (params.Nx + 1) * (params.Ny + 1) * (params.Nz + 1);
    params.stride_j         = params.Nx + 1;           % index stride along y
    params.stride_k         = (params.Nx + 1) * (params.Ny + 1); % stride along z

    if params.Nz > 1
        params.n_interior   = (params.Nx - 1) * (params.Ny - 1) * (params.Nz - 1);
    else
        params.n_interior   = (params.Nx - 1) * (params.Ny - 1);
    end

    params.n_state          = 4 * params.n_interior;    % [ψ; φ_x; φ_y; φ_z]

    % ---- Version stamp ----
    params.version = tdgl_version();

    % ---- Initialize tracking fields ----
    params.B_real_x = 0;
    params.B_real_y = 0;
    params.B_real_z = 0;
    params.B_real_x_history = [];
    params.B_real_y_history = [];
    params.B_real_z_history = [];
    params.t = [];

    % ---- CFL warning ----
    h_min = min([params.hx, params.hy, params.hz]);
    dt_cfl = h_min^2 / (4 * params.kappa^2);
    if params.dt > dt_cfl
        warning('TDGL:CFL', ...
            'dt = %.4g exceeds CFL limit %.4g for Forward Euler (h=%.2g, κ=%.2g).', ...
            params.dt, dt_cfl, h_min, params.kappa);
    end
end
