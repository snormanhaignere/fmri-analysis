function G = interp_from_surface_to_grid(surface_rh, surface_lh, patch_rh, patch_lh, roi_rh_label, roi_lh_label, grid_spacing_mm, plot_figures)
% function interp_from_surface_to_grid(surface_rh, surface_lh, patch_rh, patch_lh, roi_rh_label, roi_lh_label, grid_spacing_mm, plot_figures)
% 
% Interpolates a set of surface points to a 2-D grid, using a flattened surface patch
% 
% surface_rh, surface_lh are paths to the right and left hemispheres of a surface file in the freesurfers fsaverage template brain
% patch_rh, and patch_lh are paths to flattened surface patches
% roi_rh_label and roi_lh_label are paths to roi labels
% only surface points within constraint region are interpolated
% grid_spacing_mm defines the spacing between grid-points in mm (default is 2mm)
% plot_figures if set to true will plot some figures that are useful for making sure the interpolation has worked as expected
% 
% Example
% surface_rh = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/fla/naturalsound_us1/main_v3_combined_r1/rh.sigav_demean_300mm_hpfilt-Infsec-order2_tp3-6.mgz';
% surface_lh = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/fla/naturalsound_us1/main_v3_combined_r1/lh.sigav_demean_300mm_hpfilt-Infsec-order2_tp3-6.mgz';
% patch_rh = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/surf/rh.cortex.patch.flat';
% patch_lh = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/surf/lh.cortex2.patch.flat';
% roi_rh_label = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/rh.stp.label';
% roi_lh_label = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/lh.stp.label';
% grid_spacing_mm = 2;
% plot_figures = 1;
% interp_from_surface_to_grid(surface_rh, surface_lh, patch_rh, patch_lh, roi_rh_label, roi_lh_label, grid_spacing_mm, plot_figures);
% 
% Last modified by Sam Norman-Haignere on 12/31/2014

% directory this file is contained in
source_directory = strrep(which('interp_from_surface_to_grid.m'),'interp_from_surface_to_grid.m',''); 

% directory of matlab scripts for manipulating freesurfer files
addpath([source_directory '/fs']);

% simple check, to avoid accidentally swapping left and right hemispheres
if isempty(strfind(roi_rh_label, 'rh.')) && ~isempty(strfind(roi_rh_label, 'lh.'))
    error('Right and Left ROIs are likely swapped.');
end
if isempty(strfind(roi_lh_label, 'lh.')) && ~isempty(strfind(roi_lh_label, 'rh.'))
    error('Right and Left ROIs are likely swapped.');
end
if isempty(strfind(patch_rh, 'rh.')) && ~isempty(strfind(patch_rh, 'lh.'))
    error('Right and Left Patches are likely swapped.');
end
if isempty(strfind(patch_lh, 'lh.')) && ~isempty(strfind(patch_lh, 'rh.'))
    error('Right and Left Patches are likely swapped.');
end
if isempty(strfind(surface_rh, 'rh.')) && ~isempty(strfind(surface_rh, 'lh.'))
    error('Right and Left Surface files are likely swapped.');
end
if isempty(strfind(surface_lh, 'lh.')) && ~isempty(strfind(surface_lh, 'rh.'))
    error('Right and Left Surface files are likely swapped.');
end

% initialize
hemis = {'rh','lh'};
G.patch = cell(1,2); G.vras = cell(1,2); G.vi = cell(1,2); G.roi_mask = cell(1,2);
G.grid_x = cell(1,2); G.grid_y = cell(1,2); bounds = cell(1,2);
surf_data = cell(1,2); G.grid_data = cell(1,2);
for i = 1:2
    
    % P.patch coordinates and surface voxel indices
    if strcmp(hemis{i},'lh')
        G.patch{i} = read_patch(patch_lh);
        G.roi_mask{i} = read_label(roi_lh_label);
        surf_data{i} = MRIread(surface_lh);
    elseif strcmp(hemis{i}, 'rh')
        G.patch{i} = read_patch(patch_rh);
        G.roi_mask{i} = read_label(roi_rh_label);
        surf_data{i} = MRIread(surface_rh);
    end
    
    % vertices
    % sign of P.patch vertices indicates whether face is on border of
    % P.patch, which is not relevant for selecting the vertices
    [~,ai,ib] = intersect(abs(G.patch{i}.vnums),G.roi_mask{i}.vnums);
    G.vras{i} = G.patch{i}.vras(ai,1:2); % 2-D coordinates
    G.vi{i} = G.roi_mask{i}.vnums(ib); % indices within roi
    
    % bounds for regridding, 2 mm resolution
    bounds{i} = [floor(min(G.vras{i})); ceil(max(G.vras{i}))]; % min and max of 2-d coordinates
    [G.grid_x{i},G.grid_y{i}] = meshgrid(bounds{i}(1,1) : grid_spacing_mm : bounds{i}(2,1),  bounds{i}(1,2) : grid_spacing_mm : bounds{i}(2,2)); % regrid in units of mm
    
    % initialize voxel grid
    % vgrid{i} = nan(size(P.grid_x{i},1), size(P.grid_y{i},2), length(usubs), 165, maxruns);
        
    % matrix of voxel responses for each point in roi
    vr = squeeze(surf_data{i}.vol(:,G.vi{i}+1,:,:))';
    M = size(vr,1);
    
    % add boundary points as NaN values, useful for interpolation
    vras_boundary_pts = [G.vras{i}; bounds{i}(1,:); bounds{i}(2,:); bounds{i}(1,1), bounds{i}(2,2); bounds{i}(2,1), bounds{i}(1,2)];
    vr_boundary_pts = [vr,nan(M,4)];
    
    % grided nan mask
    mask = all(~isnan(vr_boundary_pts));
    grid_mask = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),double(mask),G.grid_x{i},G.grid_y{i},'natural');
        
    % interpolate to grid
    G.grid_data{i} = nan(size(G.grid_x{i},1), size(G.grid_y{i},2), M);
    for q = 1:M
        vr_gridded = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),vr_boundary_pts(q,:),G.grid_x{i},G.grid_y{i},'natural');
        vr_gridded(grid_mask<0.99) = NaN;
        G.grid_data{i}(:,:,q) = vr_gridded;
        
        % plot figures for the first measured response to check that things worked
        if plot_figures && q == 1
            figure;
            set(gcf,'Position',[0 0 1440 400]);
            
            subplot(1,3,1);
            colormap('jet');
            scatter(vras_boundary_pts(:,1),vras_boundary_pts(:,2),3,vr_boundary_pts(q,:));
            color_bounds = nanmedian(vr_boundary_pts(q,:)) + nanstd(vr_boundary_pts(q,:))*[-1,1]*3;
            caxis(color_bounds);
            title(sprintf('surface data, %s', hemis{i}));
            
            subplot(1,3,2);
            colormap('jet');
            scatter(G.grid_x{i}(:),G.grid_y{i}(:),3,vr_gridded(:),'LineWidth',ceil(5000/numel(vr_gridded)));
            caxis(color_bounds);
            title(sprintf('gridded data, %s', hemis{i}));
            
            subplot(1,3,3);
            colormap('jet');
            interp_back_to_grid_mask = interp2(G.grid_x{i},G.grid_y{i},double(~isnan(vr_gridded)),G.vras{i}(:,1),G.vras{i}(:,2),'linear');
            vr_gridded_withzeros = vr_gridded;
            vr_gridded_withzeros(isnan(vr_gridded_withzeros)) = 0;
            interp_back_to_grid = interp2(G.grid_x{i},G.grid_y{i},vr_gridded_withzeros,G.vras{i}(:,1),G.vras{i}(:,2));
            interp_back_to_grid(interp_back_to_grid_mask<0.001) = NaN;
            scatter(G.vras{i}(:,1),G.vras{i}(:,2),3,interp_back_to_grid);
            caxis(color_bounds);
            title(sprintf('gridded data interpolated\nback to the surface, %s', hemis{i}));
        end
    end    
end
