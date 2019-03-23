function G = matrix2grid(X, G)

% rewrap a matrix file to a grid

% dimensions of the input matrix
Xdims = size(X);

% spatial dimensions of the grid
spatial_dims_rh = size(G.grid_data{1});
spatial_dims_rh = spatial_dims_rh(1:2);
spatial_dims_lh = size(G.grid_data{2});
spatial_dims_lh = spatial_dims_lh(1:2);

% first column of the matrix should be the unwrapped and concatenated spatial
% dimensions of the right and left hemisphere
assert(Xdims(1) == (prod(spatial_dims_rh) + prod(spatial_dims_lh)));

% reshape right hemisphere
G.grid_data{1} = ...
    reshape(X(1:prod(spatial_dims_rh),:), [spatial_dims_rh, Xdims(2:end)]);

% reshape left hemisphere
G.grid_data{2} = ...
    reshape(X(prod(spatial_dims_rh)+1:end,:), [spatial_dims_lh, Xdims(2:end)]);