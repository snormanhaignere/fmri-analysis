function [patch_handle, light_handle] = plot_inflated_colordata_overlay(subjid, hemi, colordata, figh, varargin)

% Plots a surface overlay on a freesurfer-inflated brain. Similar to the function
% plot_fsaverage_1D_overlay.m, but doesn't assume the overlay is being plotted on the fsaverage
% template brain. The subjid is used to customize the camera for that particular subject.
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

addpath(genpath('/mindhive/nklab/u/svnh/fmri-analysis/code'));

% unique id
x = regexp(subjid, '_us\d*','match');
if ~isempty(x)
    us = str2double(x{1}(4:end));
end

% read vertices and faces
if strcmp(subjid, 'fsaverage')
    inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated'];
else
    switch us
        case 158
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated_dist3'];
        case 157
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated_dist4'];
        case 170
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated_dist2'];
        otherwise
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated'];
    end
end
[vertices, faces] = freesurfer_read_surf(inflated_file);
nvertices = size(vertices,1);

% set default color data based on gyral/sulcal divisions
curv = read_curv([params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.curv']);
try
    colordata(any(isnan(colordata),2) & curv>0,:) = 0.3;
    colordata(any(isnan(colordata),2) & curv<0,:) = 0.5;
catch
    colordata(any(isnan(colordata),2),:) = 0.3;
    warning('Curvature file not read properly.');
end

% create figure of specified size
if  nargin < 4 || isempty(figh)
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
        
        % set the viewing angle depending on the subject
        if strcmp(subjid, 'fsaverage')
            camup([0 0.5 1]);
            camva(4.532);
            camtarget([0 15 -10]);
            magfac = 1;
            campos([1.3604e+03 1.0309e+03 363.9286]);
        else
            switch us
                case 158
                    camup([-0.2283    0.4328    0.8721]);
                    camva(5.8442);
                    camtarget([25 15 -10]);
                    campos(1.0e+03 * [1.4598    0.7155    0.0180]);
                    magfac = 1.6;
                case 157
                    camup([-0.2413    0.4600    0.8545]);
                    camva(5.8442);
                    camtarget([35 15 -10]);
                    campos(1.0e+03 * [1.4598    0.7155    0.0180]);
                    magfac = 1.6;
                    
                case 170
                    camup([-0.1410    0.4865    0.8622]);
                    camva(4.3);
                    camtarget([10 15 -10]);
                    campos(1.0e+03 * [1.3897    0.9535   -0.3139]);
                    magfac = 1.4;
                    
                case 373
                    camup([-0.1410    0.4865    0.8622]);
                    camva(4.3);
                    camtarget([10 15 -10]);
                    campos(1.0e+03 * [1.3897    0.9535   -0.3139]);
                    magfac = 1.4;
                otherwise
                    error('No matching freesurfer subject');
            end
        end
        light_handle = camlight('headlamp','infinite');
        xlim(2.2*40*[-1 1]*magfac); ylim(2.2*65*[-1 1]*magfac); zlim(2.2*55*[-1 1]*magfac);
        drawnow;
        
    case {'lh'}
        
        % set the viewing angle depending on the subject
        if strcmp(subjid, 'fsaverage')
            camup([0 0.5 1]);
            campos([-1449.9*1 631.312*1.3 363.929*1]);
            camva(4);
            camtarget([-10 10 -5]);
            magfac = 1;
            fprintf('lh-v1\n');
        else
            switch us
                case 158
                    camup([0.3174    0.3964    0.8615]);
                    campos(1.0e+03*[-1.4818    0.7493    0.1934]);
                    camva(5.5535);
                    camtarget([-20 10 -5]);
                    magfac = 1.5;
                    
                case 157
                    camup([0.1071    0.4897    0.8653]);
                    campos(1.0e+03*[-1.5119    0.5919   -0.1478]);
                    camva(5.5535);
                    camtarget([-45 10 0]);
                    magfac = 1.5;
                    fprintf('v5');
                case 170
                    camup([0.3850    0.3973    0.8330]);
                    campos(1.0e+03*[-1.4300    0.9813    0.1788]);
                    camva(4.3);
                    camtarget([-30 10 -5]);
                    magfac = 1.5;
                case 373
                    camup([0.3850    0.3973    0.8330]);
                    campos(1.0e+03*[-1.4300    0.9813    0.1788]);
                    camva(4.3);
                    camtarget([-30 10 -5]);
                    magfac = 1.5;
                otherwise
                    error('No matching freesurfer subject');
            end
        end
        
        xlim(2.2*40*[-1 1]*magfac); ylim(2.2*65*[-1 1]*magfac); zlim(2.2*55*[-1 1]*magfac);
        light_handle = camlight('headlamp','infinite');
        
        % set(light_handle, 'Position', [-1, 1, 0.33]);
        
    otherwise
        error('hemi should be "rh" or "lh", not %s',hemi);
end
