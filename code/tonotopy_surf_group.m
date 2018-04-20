function [bestfreq_group, outputfile_colorwarp] = tonotopy_surf_group(expsub, hemi, fwhm, anova_thresh, min_subjects_per_voxel, varargin)

projfrac = [0 1];

switch expsub{1,1}
    case {'pitch_resthr_v2','pitch_resthr_v4','naturalsound'}
        runtype = 'localizer';
        model = 'block';
        freqs = [200 400 800 1600 3200 6400];
        basefreq = 100;
        freeview_control_points = [log2(200/basefreq), log2(1131/basefreq), log2(6400/basefreq)];
    otherwise
        error('No valid experiment');
end

fwhm_vol = read_smooth(expsub{1,1}, varargin{:});

% color warping to avoid green zone
inputvals = [1 3 4 6];
% outputvals = [1 2.66, 4.33, 6];
outputvals = [1 3 4 6];

emptysurf = MRIread(['~/freesurfer/fsaverage/surf/' hemi '.area.mgh']);
emptysurf.fspec = '';
emptysurf.vol = zeros(size(emptysurf.vol));
nverts = length(emptysurf.vol);

bestfreq_sub = emptysurf;
nsubs = emptysurf;

labelstring = '';
if optInputs(varargin,'mask_label')
    label = varargin{optInputs(varargin,'mask_label')+1};
    labelstring = ['_' label];
end

for i = 1:length(expsub)
    % write to file
    fprintf('%s, sub %d\n',expsub{i,1},expsub{i,2});
    subjid = [expsub{i,1} '_us' num2str(expsub{i,2})];
    runnum = read_runs(expsub{i,1},expsub{i,2},runtype,varargin{:});
    idstring = ['fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm_anovathresh' num2str(anova_thresh) labelstring];
    bestfreq_file = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '.bestfreq_' idstring '.mgh'];
    if ~exist(bestfreq_file,'file') || optInputs(varargin,'overwrite');
        tonotopy_surf(expsub{i,1}, expsub{i,2}, hemi, fwhm, anova_thresh, varargin{:});
    end
    bfs = MRIread(bestfreq_file);
    bfs.vol = bfs.vol(:);
    try
        bestfreq_sub.vol(bfs.vol > 0) = bestfreq_sub.vol(bfs.vol > 0) + bfs.vol(bfs.vol > 0)';
        nsubs.vol(bfs.vol > 0) = nsubs.vol(bfs.vol > 0) + 1;
    catch
        keyboard;
    end
end

bestfreq_group = emptysurf;
vertices_above_thresh = nsubs.vol > min_subjects_per_voxel - 1e-3;
bestfreq_group.vol(vertices_above_thresh) = bestfreq_sub.vol(vertices_above_thresh) ./ nsubs.vol(vertices_above_thresh);

% write to file
exp_unique = unique(expsub(:,1));
outputdir = ['~/freesurfer/fsaverage/group_tonotopy/' sprintf('%s_', exp_unique{:}) 'us' sprintf('%d',expsub{:,2}) '/'];
if ~exist(outputdir,'dir');
    mkdir(outputdir);
end

idstring = ['fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_anovathresh' num2str(anova_thresh) labelstring '_' runtype '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm'];
outputfile = [outputdir hemi '.bestfreq_' idstring '.mgh'];
bestfreq_group.fspec = outputfile;
MRIwrite(bestfreq_group, outputfile);

outputfile_colorwarp = [outputdir hemi '.bestfreq_' idstring '_colorwarp.mgh'];
bestfreq_group.fspec = outputfile_colorwarp;
bestfreq_group.vol(bestfreq_group.vol > 0) = interp1(inputvals,outputvals,bestfreq_group.vol(bestfreq_group.vol > 0), 'cubic');
bestfreq_group.vol(bestfreq_group.vol > 0 & bestfreq_group.vol < 2) = 2;
% bestfreq_group.vol(bestfreq_group.vol > 0 & bestfreq_group.vol > 5) = 5;
bestfreq_group.vol(120077+1) = 6; % sets voxel in the corpus collosum to this value, useful for plotting
MRIwrite(bestfreq_group, outputfile_colorwarp);

% keyboard;
% % if optInputs(varargin,'freeview');
% figure_directory = [params('rootdir') expsub{1,1} '/figures/tonotopy/' sprintf('%s_', exp_unique{:}) 'us' sprintf('%d',expsub{:,2}) '/'];
% if ~exist(figure_directory,'dir');
%     mkdir(figure_directory);
% end
% screenshot_file = [figure_directory  idstring '_' hemi '.png'];

if optInputs(varargin, 'noplot')
    return;
end
freeview3('fsaverage',hemi, 'overlay',outputfile_colorwarp, 'overlay_threshold', [1.99 3.5 6.01]);

% end

outputfile_colorwarp

% colormapping = [1 1.5 2.5 4 5 6];
% x1 = linspace(selectivity_range(1),selectivity_range(2),4);
% x2 = linspace(x1(2),x1(3),6);
% inputvals = [x1(1), x2([2,5]), x1(4), 100];
% outputvals = [x1, 100];


%     annotfile = ['/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/' hemi '.stp-pm2al-5lab-recolor.annot'];
%     inflated_file = ['~/freesurfer/fsaverage/surf/' hemi '.inflated'];
%     if strcmp(hemi,'rh')
%         cam = [200 20 20];
%     elseif strcmp(hemi,'lh')
%         cam = [-20 -20 15];
%     end
%     unix_freesurfer5p2(['freeview --surface ' inflated_file ':overlay=' outputfile ':overlay_threshold=' num2str(freeview_control_points(1)) ',' num2str(freeview_control_points(2)) ',' num2str(freeview_control_points(3)) ':annot=' annotfile ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' Zoom 1.2 &']);
