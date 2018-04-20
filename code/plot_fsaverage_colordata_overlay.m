function [patch_handle, light_handle] = plot_fsaverage_colordata_overlay(colordata, hemi, figh, varargin)

% Plots a surface overlay on the fsaverage template brain using the matlab function patch
% 
% -- Inputs --
% 
% surface_values: vector of values, one per vertex, to plot, NaN values are not plotted
% 
% hemi: whether the left or right hemisphere is being plotted
% 
% color_map_to_plot: the name of the colormap to use (default is 'parula' if not specified), can
% also be a N x 3 matrix of values to interpolate within the color range.
% 
% color_range: range of values to plot, [lower_bound, upperbound] (if not specified the central 95 of the distribution of values is plotted)
% 
% -- Outputs --
% 
% patch_handle: handle to the patch object created
% 
% light_handle: handle to the light object created
% 
% -- Example: Plots significance map for music component --  
% 
% hemi = 'rh';
% surf = MRIread([params('rootdir') 'fmri-analysis/test_data/' hemi '.ICA_pmap_music.mgz']);
% surface_values = surf.vol;
% colormapname = 'parula';
% plot_fsaverage_surface(surface_values, hemi, colormapname);
% 
% Created by Sam NH on 8/12/2015

% read vertices and faces
% nvertices = 163842; 
[vertices, faces] = freesurfer_read_surf(ref_image);

% create figure of specified size
if  nargin < 3 || isempty(figh)
    figh = figure;
    pos = get(figh,'Position');
    set(figh, 'Position', [pos(1:2), 800 800]);
else
    clf(figh);
end

% create the patch object
patch_handle = patch('vertices', vertices, 'Faces', faces, 'FaceVertexCData', colordata,'FaceLighting','gouraud','SpecularStrength',0,'DiffuseStrength',0.7);
shading interp;

% adjust viewing angle 
switch hemi
    case {'rh'}
%         view([115 10])
        camup([0 0.5 1]);
        camva(4.532)
        campos([1.3604e+03 1.0309e+03 363.9286]);
        camtarget([0 15 -10]);
        xlim(2.2*40*[-1 1]); ylim(2.2*65*[-1 1]); zlim(2.2*55*[-1 1]);
        light_handle = camlight('right','infinite');
        set(light_handle, 'Position', [0.33 0.33 0.33]);

    case {'lh'}
%         view([-105 0]);
        camup([0 0.5 1]);
        % camzoom(2.2);
        camtarget([-10 10 -5]);
        %         campos([-1449.9 631.312 363.929]);
        campos([-1449.9*1 631.312*1.3 363.929*1]);
        camva(4);
        % campos([1.3604e3 1.0309e3 0.3639e3])
        xlim(2.2*40*[-1 1]); ylim(2.2*65*[-1 1]); zlim(2.2*55*[-1 1]);
        light_handle = camlight('left','infinite');
        set(light_handle, 'Position', [-1, 1, 0.33]);
        
    otherwise
        error('hemi should be "rh" or "lh", not %s',hemi);
end
