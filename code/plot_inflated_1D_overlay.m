function [color_data, patch_handle, light_handle] = ...
    plot_inflated_1D_overlay(subjid, hemi, surface_values, color_map_to_plot, color_range, figh, varargin)

% function [color_data, patch_handle, light_handle] = plot_inflated_1D_overlay(subjid, surface_values, hemi, color_map_to_plot, color_range, figh, varargin)
%
% Plots a surface overlay on a freesurfer-inflated brain. Similar to the function
% plot_fsaverage_1D_overlay.m, but doesn't assume the overlay is being plotted on the fsaverage
% template brain. The subjid is used to customize the camera for that particular subject.
%
% -- Inputs --
%
% subjid: freesurfer subject id
%
% surface_values: vector of values, one per vertex, to plot, NaN values are not plotted
%
% hemi: whether the left or right hemisphere is being plotted
%
% color_map_to_plot (optional): the name of the colormap to use (default is 'parula' if not specified), can
% also be a N x 3 matrix of values to interpolate within the color range.
%
% color_range (optional): range of values to plot, [lower_bound, upperbound] (if not specified the central 95 of the distribution of values is plotted)
%
% figh (optional): matlab handle of the figure to plot the surface in, if unspecified a new figure
% handle is created (i.e. figh = figure)
%
% -- Outputs --
%
% color_data: N x 3 matrix specifying the RGB color of each vertex
%
% patch_handle: handle to the patch object created
%
% light_handle: handle to the light object created
%
% -- Example: Plots significance map for music component discovered by ICA --
%
% hemi = 'rh';
% surf = MRIread(['/mindhive/nklab/u/svnh/fmri-analysis/test_data/' hemi '.ICA_pmap_music.mgz']);
% surface_values = surf.vol;
% colormapname = 'parula';
% plot_inflated_1D_overlay('fsaverage', surface_values, hemi, colormapname);
%
% Modified by Sam NH on 9/3/2015

addpath(genpath('/mindhive/nklab/u/svnh/fmri-analysis/code'));

% unique id
x = regexp(subjid, '_us\d*','match');
if ~isempty(x)
    us = str2double(x{1}(4:end));
end

% default colormap is parula
if nargin < 4
    color_map_to_plot = 'parula';
end

% by default plots the central 95% of the distribution
if nargin < 5
    if any(~isnan(surface_values))
        [Nx,x] = hist(surface_values(~isnan(surface_values) & surface_values ~= 0),1000);
        Cx = cumsum(Nx/sum(Nx));
        [~,xi] = unique(Cx);
        x = x(xi);
        Cx = Cx(xi);
        color_range = interp1(Cx,x,[0.025 0.975]);
    else
        color_range = [];
    end
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
color_data = 0.5*ones(nvertices,3);
if length(curv) == size(color_data,1)
    color_data(curv>0,:) = 0.3;
end

% read in a colormap
if ischar(color_map_to_plot)
    h = figure;
    cmap = colormap(color_map_to_plot);
    close(h);
else
    cmap = color_map_to_plot;
end

if ~isempty(color_range)
    % interpolate surface values to colormap
    surface_values_bounded = surface_values;
    surface_values_bounded(surface_values_bounded < color_range(1)) = color_range(1);
    surface_values_bounded(surface_values_bounded > color_range(2)) = color_range(2);
    x = linspace(color_range(1),color_range(2),size(cmap,1))';
    for i = 1:3
        color_data(~isnan(surface_values_bounded),i) = interp1(x, cmap(:,i), surface_values_bounded(~isnan(surface_values_bounded))','pchip');
    end
end

% return after calculating color data without plotting
if optInputs(varargin,'noplot')
    return;
end

% create figure of specified size
if nargin < 6 || isempty(figh)
    figh = figure;
    pos = get(figh,'Position');
    set(figh, 'Position', [pos(1:2), 800 800]);
    clf(figh);
else
    clf(figh);
end

% create the patch object
patch_handle = patch('vertices', vertices, 'Faces', faces, 'FaceVertexCData', color_data,'FaceLighting','gouraud','SpecularStrength',0,'DiffuseStrength',0.7);
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
        
    otherwise
        error('hemi should be "rh" or "lh", not %s',hemi);
end

if ~isempty(color_range)
    
    % colorbar
    colormap(color_map_to_plot);
    colorbar_handle = colorbar('Location','South');
    colorbar_labels_str = cell(1,5);
    colorbar_labels_num = linspace(color_range(1),color_range(2),5);
    for i = 1:5
        num_sig_digits = round(log10(colorbar_labels_num(i)));
        if num_sig_digits > 3 || num_sig_digits < -3
            colorbar_labels_str{i} = sprintf('%.2e',colorbar_labels_num(i));
        else
            colorbar_labels_str{i} = sprintf('%.2f',colorbar_labels_num(i));
        end
    end
    set(colorbar_handle, 'XTick', linspace(0,1,5), 'XTickLabel', colorbar_labels_str,'FontSize',20,'Position',[0.1469 0.05 0.7438 0.0312]);
end