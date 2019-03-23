function Y = fillin_NaN(X, elements_without_NaN, dim_to_adjust)

% -- Examples --
% 
% fillin_NaN(1:3, [1 0 0 1 1], 2)
% fillin_NaN(rand(2,3), [0 0 1 1 0 0], 1)
% 
% Helper function for glm_event_regression.m and sigav_glm.m (amongst others)
% 
% 2016-09-21: Renamed, Sam NH

assert(size(X,dim_to_adjust) == sum(elements_without_NaN));
dims = size(X);
dims(dim_to_adjust) = length(elements_without_NaN);

xi = cell(1,length(dims));
for i = 1:length(dims)
    if i == dim_to_adjust
        xi{i} = find(elements_without_NaN);
    else
        xi{i} = 1:dims(i);
    end
end

Y = nan(dims);
Y(xi{:}) = X;