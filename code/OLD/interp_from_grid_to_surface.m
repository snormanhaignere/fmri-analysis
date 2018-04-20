function Y_surf = interp_from_grid_to_surface(Y, voxel_selection, grid_x, grid_y, vi, vras, fwhm_surf, varargin)
% function Y_surf = interp_from_grid_to_surface(Y, voxel_selection, grid_x, grid_y, vi, vras, fwhm_surf, varargin)
% 
% Interpolates a matrix of voxel responses (Y: M conditions by N voxels * S subjects) to the cortical surface in freesurfer.
% Inverse of interp_from_surface_to_grid function.
% 
% voxel_selection is a vector that was used to select voxels from the original voxel grid returned by interp_from_surface_to_grid
% 
% grid_x and grid_y are matrices specifying the RAS positions of the voxels (1 matrix per hemisphere)
% vras is a vector of RAS positions of the target surface vertices (1 per hemisphere)
% vi is a vector of indices for each RAS position, which should match those of the fsaverage template brain (1 per hemisphere)
% the responses can be optionally smoothed, by specifying a smoothing kernel, fwhm_surf in mm
% 
% If smoothing is applied voxels with NaN values are set to zero prior to smoothing.
% It is also possible to just smooth the voxel with NaN values, thereby interpolating undefined points, 
% by adding the flag 'just_smooth_NaNs' as an optional argument.
% 
% Modified on 1/5/2014 to use 3D information in the interpolation
% 
% Last edited by Sam Norman-Haigner on 1/5/2014

% dimensions of the matrix
[M,N] = size(Y);

% parameters
nsurfpts = 163842;

% number of voxels per hemisphere and subject
n_voxels_rh = size(grid_x{1},1)*size(grid_x{1},2);
n_voxels_lh = size(grid_x{2},1)*size(grid_x{2},2);
n_voxels = n_voxels_lh + n_voxels_rh;

% number of subjects
n_subjects = length(voxel_selection) / n_voxels;

% grid voxel weights
Y_gridded = cell(1,2);
Y_gridded{1} = nan([size(grid_x{1}), n_subjects, M]);
Y_gridded{2} = nan([size(grid_x{2}), n_subjects, M]);
for k = 1:M
    % n_voxels x subject voxel weight matrix
    x = nan(n_voxels*n_subjects,1);
    x(voxel_selection) = Y(k,:)';
    voxel_weights_reshaped = reshape(x, [n_voxels, n_subjects]);
    hemis = {'rh','lh'};
    for i = 1:length(hemis)
        % convert to grid: xgrid x ygrid x n_subjects
        switch hemis{i}
            case 'rh'
                Y_gridded{i}(:,:,:,k) = reshape(voxel_weights_reshaped(1:n_voxels_rh, :), [size(grid_x{1}), n_subjects]);
                
            case 'lh'
                Y_gridded{i}(:,:,:,k) = reshape(voxel_weights_reshaped((1:n_voxels_lh) + n_voxels_rh, :), [size(grid_x{2}), n_subjects]);
        end
    end
end

% optionally smooth on grid
if fwhm_surf > 0
    
    kernel_size = round(3*fwhm_surf/2); % kernel is 3x FWHM
    kernel_size = kernel_size-mod(kernel_size,2)+1; % make odd
    gaussian_kernel = fspecial('gaussian', [kernel_size,kernel_size], fwhm2std(fwhm_surf/2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
    
    Y_gridded_smoothed = cell(1,2);
    Y_gridded_smoothed{1} = nan(size(Y_gridded{1}));
    Y_gridded_smoothed{2} = nan(size(Y_gridded{2}));
    for k = 1:M
        for i = 1:length(hemis{i})
            % interpolate to surface
            for j = 1:n_subjects
                if optInputs(varargin, 'just_smooth_NaNs')
                    x1 = Y_gridded{i}(:,:,j,k);
                    x2 = conv2_setNaNs_tozero(x1, gaussian_kernel);
                    x1(isnan(x1)) = x2(isnan(x1));
                    Y_gridded_smoothed{i}(:,:,j,k) = x1;
                else
                    Y_gridded_smoothed{i}(:,:,j,k) = conv2_setNaNs_tozero(Y_gridded{i}(:,:,j,k), gaussian_kernel);
                end
            end
        end
    end
    Y_gridded = Y_gridded_smoothed;
end

% interpolate voxel weights to surface
Y_surf = nan(M, nsurfpts*2, n_subjects);
for k = 1:M
    for i = 1:length(hemis{i})
        for j = 1:n_subjects
            Y_surf(k,vi{i} + 1 + (i-1)*nsurfpts,j) = interp2(grid_x{i},grid_y{i},Y_gridded{i}(:,:,j,k),vras{i}(:,1),vras{i}(:,2));
        end
    end
end