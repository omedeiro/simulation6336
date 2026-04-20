function state = create_initial_state(params, varargin)
% CREATE_INITIAL_STATE  Build initial state vector for TDGL simulation.
%
%   state = create_initial_state(params)
%   state = create_initial_state(params, 'perturbation', 0.01)
%
%   Creates a column vector [psi; phi_x; phi_y; phi_z] of length 4*n_interior.
%
%   Default initial condition:
%       psi   = 1  (uniform superconducting state)
%       phi_x = 0
%       phi_y = 0
%       phi_z = 0
%
%   Optional Name-Value Pairs:
%     'perturbation'  : double (default 0). Amplitude of random real
%                       perturbation added to psi.
%     'seed'          : non-negative integer (default []). RNG seed for
%                       reproducibility. Leave empty for non-deterministic.
%
%   See also SETUP_PARAMETERS, CONSTRUCT_GRID_INDICES

    ip = inputParser;
    addRequired(ip,  'params', @isstruct);
    addParameter(ip, 'perturbation', 0,  @(x) isnumeric(x) && x >= 0);
    addParameter(ip, 'seed',         [], @(x) isempty(x) || (isnumeric(x) && x >= 0));
    parse(ip, params, varargin{:});

    n_int = params.n_interior;
    pert  = ip.Results.perturbation;
    seed  = ip.Results.seed;

    % Seed RNG if requested
    if ~isempty(seed)
        rng(seed);
    end

    % Order parameter: start in superconducting state
    psi_init = ones(n_int, 1);
    if pert > 0
        psi_init = psi_init + pert * (2 * rand(n_int, 1) - 1);
    end

    % Vector potential components: start at zero
    phi_x_init = zeros(n_int, 1);
    phi_y_init = zeros(n_int, 1);
    phi_z_init = zeros(n_int, 1);

    % Assemble state vector
    state = [psi_init; phi_x_init; phi_y_init; phi_z_init];
end
