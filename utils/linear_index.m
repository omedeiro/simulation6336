function m = linear_index(i, j, k, Nx, Ny)
% LINEAR_INDEX  Convert (i,j,k) grid coordinates to a flat linear index.
%
%   m = linear_index(i, j, k, Nx, Ny)
%
%   Computes the 1-based linear index for node (i, j, k) on a grid of
%   size (Nx+1) × (Ny+1) × (Nz+1):
%
%       m = i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*(Ny+1)
%
%   If i, j, or k is out of bounds it is clamped to [1, Nx+1], [1, Ny+1],
%   or [1, ∞] respectively (for safe neighbor lookups at boundaries).
%
%   Inputs:
%     i  : x-index (1-based), scalar or array
%     j  : y-index (1-based), scalar or array
%     k  : z-index (1-based), scalar or array
%     Nx : number of interior x-intervals (grid has Nx+1 x-nodes)
%     Ny : number of interior y-intervals (grid has Ny+1 y-nodes)
%
%   Output:
%     m  : linear index (same size as input arrays)
%
%   See also CONSTRUCT_GRID_INDICES

    % Clamp to valid range
    i = max(1, min(i, Nx + 1));
    j = max(1, min(j, Ny + 1));
    k = max(1, k);  % no upper bound enforced for k

    stride_j = Nx + 1;
    stride_k = (Nx + 1) * (Ny + 1);

    m = i + (j - 1) .* stride_j + (k - 1) .* stride_k;
end
