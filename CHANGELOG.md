# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.2.0] — 2026-04-12

### Added
- **Modular directory structure**: `core/`, `operators/`, `physics/`, `solvers/`,
  `visualization/`, `utils/`, `examples/`, `analysis/`, `tests/`, `archive/`.
- `setup_paths.m` — single script to configure MATLAB path.
- `core/setup_parameters.m` — `inputParser`-based parameter initialization with
  validation, CFL warnings, and explicit variable names.
- `core/construct_grid_indices.m` — documented rewrite of index construction with
  descriptive field names (`interior_to_full`, `x_face_lo_inner`, etc.).
- `core/create_initial_state.m` — reproducible initial condition generator.
- `core/tdgl_version.m` — programmatic version query.
- `operators/construct_lpsi_{x,y,z}.m` — covariant Laplacian operators for ψ.
- `operators/construct_lphi_{x,y,z}.m` — standard Laplacian operators for φ.
- `operators/construct_fpsi.m` — nonlinear Ginzburg-Landau forcing term.
- `operators/construct_fphi_{x,y,z}.m` — curl-curl + supercurrent forcing.
- `physics/evaluate_rhs.m` — full TDGL right-hand side with clear sectioned logic.
- `physics/apply_boundary_conditions.m` — zero-current and magnetic field BCs.
- `physics/evaluate_applied_field.m` — time-dependent field ramp.
- `physics/evaluate_bfield.m` — B = curl(A) computation.
- `solvers/forward_euler.m` — explicit Euler integrator.
- `solvers/trapezoidal_solver.m` — implicit trapezoidal with adaptive dt.
- `solvers/newton_gcr_trapezoidal.m` — Newton solver for implicit system.
- `solvers/tgcr_matrix_free.m` — matrix-free GCR linear solver.
- `solvers/tgcr_matrix_free_trapezoidal.m` — GCR for trapezoidal Jacobian.
- `visualization/plot_fields_2d.m` — 2D imagesc with GIF export.
- `visualization/plot_fields_3d.m` — 3D scatter plots.
- `utils/cube_to_column.m`, `utils/column_to_cube.m`, `utils/linear_index.m`.
- `examples/run_trapezoidal.m`, `examples/run_forward_euler.m`.
- `VERSION` file and `CHANGELOG.md`.

### Changed
- All variables renamed to be descriptive: `psi`, `phi_x`, `interior_to_full`
  instead of `x`, `y1`, `M2`.
- All functions documented with full help blocks.
- Constant φ Laplacians precomputed once instead of rebuilt every RHS call.

### Deprecated
- Original flat-directory `.m` files moved to `archive/`. These will be removed
  in a future major release.

## [0.1.0] — 2021-12-06

### Added
- Initial MATLAB implementation of 3D TDGL equations.
- Forward Euler and trapezoidal time integrators.
- Newton-GCR nonlinear solver with matrix-free directional derivatives.
- Zero-current and magnetic field boundary conditions.
- 2D and 3D visualization with GIF export.
- GUI prototype.
