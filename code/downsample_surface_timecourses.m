function downsample_surface_timecourses(exp, us, runtype, r, input_fname, grid_spacing_mm, plot_figures, varargin)
% function downsample_surface_timecourses(exp, us, runtype, r, input_fname, grid_spacing_mm, plot_figures, varargin)
% 
% Downsamples preprocessed timecourse files to a grid on the flattened
% cortical surface. Useful because Freesurfer uses an excessively fine
% tessellation for most purposes. 
% 
% The function relies on the subfunction "interp_from_surface_to_grid.m"
% to perform the interpolation. 
% 
% Example human subject:
% exp = 'pitch_localizer_human';
% us = 1;
% runtype = 'pitchloc2_combined_split';
% r = 1;
% input_fname = 'motcorr_smooth300mm';
% grid_spacing_mm = 1.5;
% varargin = {'overwrite'};
% downsample_surface_timecourses(exp, us, runtype, r, input_fname, grid_spacing_mm, plot_figures, varargin{:})
% 
% Example monkey subject:
% exp = 'pitch_localizer_monkey';
% us = 158;
% runtype = 'pitchloc2_combined_split';
% r = 10;
% input_fname = 'motcorr_smooth286mm';
% grid_spacing_mm = 2.86*0.5;
% varargin = {'monkey','overwrite'};
% downsample_surface_timecourses(exp, us, runtype, r, input_fname, grid_spacing_mm, plot_figures, varargin{:})

% preprocessing directory
subjid = [exp '_us' num2str(us)];
if optInputs(varargin, 'monkey')
    preprocdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
else
    preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
end

% roi name
if optInputs(varargin, 'monkey')
    roi = 'hand-audctx';
else
    roi = 'hand-stp-stg';
end
if optInputs(varargin, 'roi')
    roi = varargin{optInputs(varargin, 'roi')+1};
end

% flattened patches and constraint ROIs for monkeys and humans
if optInputs(varargin, 'monkey')
    patch_rh = ['/mindhive/nklab/u/svnh/freesurfer/' subjid '/surf/rh.cortex.patch.flat'];
    patch_lh = ['/mindhive/nklab/u/svnh/freesurfer/' subjid '/surf/lh.cortex.patch.flat'];
    roi_rh_label = ['/mindhive/nklab/u/svnh/freesurfer/' subjid '/label/rh.' roi '.label'];
    roi_lh_label = ['/mindhive/nklab/u/svnh/freesurfer/' subjid '/label/lh.' roi '.label'];
else
    patch_rh = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/surf/rh.cortex.patch.flat';
    patch_lh = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/surf/lh.cortex2.patch.flat';
    roi_rh_label = ['/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/rh.' roi '.label'];
    roi_lh_label = ['/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/lh.' roi '.label'];
end

% mat file with gridded values
grid_mat_file = [preprocdir input_fname '_grid_' roi '_' num2str(grid_spacing_mm) 'mm.mat'];
if ~exist(grid_mat_file, 'file') || optInputs(varargin, 'overwrite')
    
    tic;

    % surface files using Freesurfer's fine tessellation
    surface_rh = [preprocdir 'rh.' input_fname '.mgz'];
    surface_lh = [preprocdir 'lh.' input_fname '.mgz'];
    
    
    % interpolate to the grid
    G = interp_from_surface_to_grid(surface_rh, surface_lh, patch_rh, patch_lh, roi_rh_label, roi_lh_label, grid_spacing_mm, plot_figures);
    
    % save to the mat file
    save(grid_mat_file, 'G');
    
    toc;
end