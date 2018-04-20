function tonotopy_selectivity_group(expsub, hemi, fwhm, con, pthresh, min_subjects_per_voxel, selectivity_range, varargin)

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

% colormapping = [1 1.5 2.5 4 5 6];
fwhm_vol = read_smooth(expsub{1,1}, varargin{:});
emptysurf = MRIread(['~/freesurfer/fsaverage/surf/' hemi '.area.mgh']);
emptysurf.fspec = '';
emptysurf.vol = zeros(size(emptysurf.vol));
nverts = length(emptysurf.vol);

freqsel_sub = emptysurf;
nsubs = emptysurf;

labelstring = '';
if optInputs(varargin,'mask_label')
    label = varargin{optInputs(varargin,'mask_label')+1};
    lb = read_label_SNH(['~/freesurfer/fsaverage/label/' hemi '.' label '.label']);
    freqsel.vol(lb.vnums+1) = 0;
    labelstring = ['_' label];
end

for i = 1:length(expsub)
    % write to file
    fprintf('%s, sub %d\n',expsub{i,1},expsub{i,2});
    subjid = [expsub{i,1} '_us' num2str(expsub{i,2})];
    runnum = read_runs(expsub{i,1},expsub{i,2},runtype,varargin{:});
    idstring = ['fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm_' con '_pthresh' num2str(pthresh) labelstring];
    freqsel_file = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '.freqsel_' idstring '.mgh'];
    if ~exist(freqsel_file,'file') || optInputs(varargin,'overwrite');
        tonotopy_selectivity(expsub{i,1}, expsub{i,2}, hemi, fwhm, con, pthresh, selectivity_range, varargin{:});
    end
    freqsel = MRIread(freqsel_file);
    freqsel.vol = freqsel.vol(:);
    freqsel_sub.vol(freqsel.vol > 1e-3) = freqsel_sub.vol(freqsel.vol > 1e-3) + freqsel.vol(freqsel.vol > 1e-3)';
    nsubs.vol(freqsel.vol > 1e-3) = nsubs.vol(freqsel.vol > 1e-3) + 1;
end

freqsel_group = emptysurf;
vertices_above_thresh = nsubs.vol > min_subjects_per_voxel - 1e-3;
freqsel_group.vol(vertices_above_thresh) = freqsel_sub.vol(vertices_above_thresh) ./ nsubs.vol(vertices_above_thresh);
freqsel_group.vol(freqsel_group.vol > 1e-3 & freqsel_group.vol < selectivity_range(1)) = selectivity_range(1);

% write to file
exp_unique = unique(expsub(:,1));
outputdir = ['~/freesurfer/fsaverage/group_tonotopy/' sprintf('%s_', exp_unique{:}) 'us' sprintf('%d',expsub{:,2}) '/'];
if ~exist(outputdir,'dir');
    mkdir(outputdir);
end

idstring = ['fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm_' con '_pthresh' num2str(pthresh) '_selectivity_range' num2str(selectivity_range(1)) '-' num2str(selectivity_range(2))];
outputfile = [outputdir hemi '.freqsel_' idstring '.mgh'];
freqsel_group.fspec = outputfile;
MRIwrite(freqsel_group, outputfile);

% keyboard;
outputfile_colorwarp = [outputdir hemi '.freqsel_' idstring '_colorwarp.mgh'];
freqsel_group.fspec = outputfile_colorwarp;
x1 = linspace(selectivity_range(1),selectivity_range(2),4);
x2 = linspace(x1(2),x1(3),6);
inputvals = [x1(1), x2([2,5]), x1(4), max(freqsel_group.vol)];
outputvals = [x1, max(freqsel_group.vol)];
% y = [x(1), linspace(x(2),x(4),4), x(end)]
freqsel_group.vol(freqsel_group.vol > 0) = interp1(inputvals,outputvals,freqsel_group.vol(freqsel_group.vol > 0), 'cubic');
MRIwrite(freqsel_group, outputfile_colorwarp);

if optInputs(varargin,'freeview');
    %freeview(outputfile_colorwarp, hemi, [selectivity_range(1)-0.01, mean(selectivity_range), selectivity_range(2)]);
    freeview3('fsaverage',hemi, 'overlay',outputfile_colorwarp, 'overlay_threshold', [selectivity_range(1)-0.01, mean(selectivity_range), selectivity_range(2)]);
end

%     annotfile = ['/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/' hemi '.stp-pm2al-5lab-recolor.annot'];
%     inflated_file = ['~/freesurfer/fsaverage/surf/' hemi '.inflated'];
%     if strcmp(hemi,'rh')
%         cam = [200 20 20];
%     elseif strcmp(hemi,'lh')
%         cam = [-20 -20 15];
%     end
%     unix_freesurfer5p2(['freeview --surface ' inflated_file ':overlay=' outputfile ':overlay_threshold=' num2str(freeview_control_points(1)) ',' num2str(freeview_control_points(2)) ',' num2str(freeview_control_points(3)) ':annot=' annotfile ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' Zoom 1.2 &']);
