# simulation6336 — 3D TDGL Superconductor Vortex Simulator

[![Version](https://img.shields.io/badge/version-0.2.0-blue.svg)](VERSION)

MATLAB implementation of the 3D Time-Dependent Ginzburg-Landau (TDGL) equations
for simulating magnetic vortex dynamics in Type-II superconductors.

## Quick Start

```matlab
% 1. Add paths
run('setup_paths.m')

% 2. Set up parameters
params = setup_parameters('Nx', 60, 'Ny', 60, 'Nz', 3, 'kappa', 5, 'applied_Bz', 0.6);
params = construct_grid_indices(params);

% 3. Create initial state and run
X0 = create_initial_state(params, 'perturbation', 0.01);
[X, params] = trapezoidal_solver(@evaluate_rhs, X0, params, @evaluate_applied_field, 0, 20, 0.1, []);

% 4. Visualize
plot_fields_2d(X, params);
```

See `examples/run_trapezoidal.m` and `examples/run_forward_euler.m` for complete runnable demos.

## Equations

The TDGL system in the zero-electric-potential gauge:

```
∂ψ/∂t = (∇ − iA)²ψ + (1 − |ψ|²)ψ
∂A/∂t = κ² ∇×(∇×A) − Im[ψ*(∇ − iA)ψ]
```

where ψ is the superconducting order parameter, A = (φ_x, φ_y, φ_z) is the
magnetic vector potential, and κ is the Ginzburg-Landau parameter.

## Directory Structure

```
simulation6336/
├── setup_paths.m               — Add all subdirectories to MATLAB path
├── VERSION                     — Semantic version (0.2.0)
├── CHANGELOG.md                — Release history
├── README.md
├── LICENSE
│
├── core/                       — Parameter setup and grid construction
│   ├── setup_parameters.m      — inputParser-based parameter initialization
│   ├── construct_grid_indices.m — Build interior ↔ full-grid index maps
│   ├── create_initial_state.m  — Generate initial state vector
│   └── tdgl_version.m          — Return current version string
│
├── operators/                  — Sparse matrix operator construction
│   ├── construct_lpsi_{x,y,z}.m — Covariant Laplacian for ψ (Peierls phases)
│   ├── construct_lphi_{x,y,z}.m — Standard Laplacian for φ components (κ²/h²)
│   ├── construct_fpsi.m         — Nonlinear GL term: (1-|ψ|²)ψ
│   └── construct_fphi_{x,y,z}.m — Curl-curl + supercurrent forcing
│
├── physics/                    — RHS evaluation and boundary conditions
│   ├── evaluate_rhs.m          — Full TDGL right-hand side dX/dt = f(X)
│   ├── apply_boundary_conditions.m — Zero-current + magnetic field BCs
│   ├── evaluate_applied_field.m    — Time-dependent applied field ramp
│   └── evaluate_bfield.m          — B = curl(A) from state
│
├── solvers/                    — Time integration and Newton solvers
│   ├── forward_euler.m             — Explicit Euler integrator
│   ├── trapezoidal_solver.m        — Implicit trapezoidal with Newton-GCR
│   ├── newton_gcr_trapezoidal.m    — Newton iteration for trapezoidal system
│   ├── tgcr_matrix_free.m          — GCR linear solver (matrix-free)
│   └── tgcr_matrix_free_trapezoidal.m — GCR for trapezoidal Jacobian
│
├── visualization/              — Plotting and animation
│   ├── plot_fields_2d.m        — 2D imagesc plots (|ψ|², Bz) with GIF export
│   └── plot_fields_3d.m        — 3D scatter plots (|ψ|², Bx, By, Bz)
│
├── utils/                      — Small utility functions
│   ├── cube_to_column.m        — Reshape 3D array → column vector
│   ├── column_to_cube.m        — Reshape column vector → 3D array
│   └── linear_index.m          — (i,j,k) → flat linear index
│
├── examples/                   — Runnable demo scripts
│   ├── run_trapezoidal.m       — 60×60×3 thin film with trapezoidal
│   └── run_forward_euler.m     — 30×30×30 cube with forward Euler
│
├── archive/                    — Original (pre-restructuring) files
│   ├── contruct_indices.m      — Original index builder
│   ├── eval_f.m                — Original RHS
│   ├── Trapezoidal.m           — Original solver
│   └── ...                     — All other original .m files
│
├── gifs/                       — Saved GIF animations
├── export/                     — Exported data files
└── GUI/                        — GUI-specific files (legacy)
```

## State Vector Layout

The state vector `X` has length `4 * n_interior` where
`n_interior = (Nx-1) * (Ny-1) * (Nz-1)`:

| Segment | Indices | Content |
|---------|---------|---------|
| 1 | `1 : n_int` | ψ (complex order parameter) |
| 2 | `n_int+1 : 2*n_int` | φ_x (vector potential, x) |
| 3 | `2*n_int+1 : 3*n_int` | φ_y (vector potential, y) |
| 4 | `3*n_int+1 : 4*n_int` | φ_z (vector potential, z) |

## Grid Conventions

- Full grid: `(Nx+1) × (Ny+1) × (Nz+1)` nodes
- Interior grid: `(Nx-1) × (Ny-1) × (Nz-1)` nodes (excludes boundaries)
- Linear index: `m = i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*(Ny+1)` (1-based)
- Boundary conditions applied on face nodes

## Key Improvements Over Original Code

- **Modular directory structure** — functions grouped by purpose
- **Explicit variable names** — `psi`, `phi_x`, `interior_to_full` instead of `x`, `y1`, `M2`
- **Full documentation** — every function has a complete help block
- **Input validation** — `inputParser` with type/range checking
- **Consistent naming** — `snake_case` for all new files and functions
- **Precomputed operators** — constant Laplacians built once, not per RHS call
- **Original files preserved** — everything in `archive/` for reference

## Versioning

This project uses [Semantic Versioning](https://semver.org/). The current
version is stored in the `VERSION` file at the project root and can be
queried programmatically:

```matlab
>> tdgl_version()
ans = '0.2.0'
```

See [CHANGELOG.md](CHANGELOG.md) for release history.

