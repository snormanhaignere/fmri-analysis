function [bfmap, bestfreq_smooth_thresh_colorwarp_file] = tonotopy_surf(exp, us, hemi, fwhm, anova_thresh, varargin)

% Parameters
% anova_thresh = 3;
% pval_thresh_direct_contrasts = 1.3;
projfrac = [0 1];
path_to_this_file = which('tonotopy_surf');
addpath(strrep(path_to_this_file, 'tonotopy_surf.m','fs'));

freesurfer_version = read_freesurfer_version(exp);

switch exp
    case {'pitch_resthr_v2','pitch_resthr_v4','naturalsound'}
        runtype = 'localizer';
        model = 'block';
        freqs = [200 400 800 1600 3200 6400];
        condlevels = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
        basefreq = 100;
        freeview_control_points = [log2(200/basefreq), log2(1131/basefreq), log2(6400/basefreq)];
    case {'tono-pitch-localizer'}
        runtype = 'tono_pitch_localizer';
        model = 'block';
        freqs = [200 400 800 1600 3200 6400];
        condlevels = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
        basefreq = 100;
        freeview_control_points = [log2(200/basefreq), log2(1131/basefreq), log2(6400/basefreq)];
    case {'tono-localizer'}
        runtype = 'tono_localizer';
        model = 'block';
        freqs = [200 400 800 1600 3200 6400];
        condlevels = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
        basefreq = 100;
        freeview_control_points = [log2(200/basefreq), log2(1131/basefreq), log2(6400/basefreq)];
    otherwise
        error('No valid experiment');
end

fwhm_vol = read_smooth(exp, varargin{:});
inputvals = [1 3 4 6];
% outputvals = [1 2.66, 4.33, 6];
outputvals = [1 3, 4, 6];

subjid = [exp '_us' num2str(us)];
runnum = read_runs(exp,us,runtype,varargin{:});

% directory used to store first level surface files
sladir = ['~/freesurfer/' subjid '/sla/'];
if ~exist(sladir,'dir');
    mkdir(sladir);
end

emptysurf = MRIread(['~/freesurfer/fsaverage/surf/' hemi '.area.mgh']);
emptysurf.fspec = '';
emptysurf.vol = zeros(size(emptysurf.vol));
nverts = length(emptysurf.vol);

% Betas for each condition
betas = nan(nverts, length(freqs));
for i = 1:length(freqs)
    idstring = ['freq' num2str(freqs(i)) '_fwhm0_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm'];
    glmdir = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '_' idstring '_r' sprintf('%d',runnum)];
    allcopes_file = [glmdir '/allcopes.mgz'];
    if ~exist(allcopes_file,'file') || optInputs(varargin,'overwrite')
        sla_surf(exp,us,['freq' num2str(freqs(i))],0,[],[],[],[],hemi,runtype,'block','nocorrection','fsaverage',varargin{:});
    end
    cope = MRIread([glmdir '/allcopes.mgz']);
    betas(:,i) = squeeze(mean(cope.vol,4));
end

% Maximum beta for each vertex
[~,xi] = max(betas,[],2);
maxfreq = xi;
maxfreq(all(betas==0,2)) = 0;

% write to file
idstring = ['fwhm0_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm'];
bestfreq_file = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '.bestfreq_' idstring '.mgh'];
bestfreq = emptysurf;
bestfreq.vol = maxfreq;
bestfreq.fspec = bestfreq_file;
MRIwrite(bestfreq, bestfreq_file);

if fwhm > 0;
    bestfreq_smooth_file = strrep(bestfreq_file,'fwhm0',['fwhm' num2str(fwhm)]);
    unix_freesurfer_version(freesurfer_version, ['mri_surf2surf --s fsaverage --sval ' bestfreq_file ' --tval ' bestfreq_smooth_file ' --hemi ' hemi ' --fwhm-src ' num2str(fwhm)]);
else
    bestfreq_smooth_file = bestfreq_file;
end
% anova_surf_subjrecon = [sladir hemi '.' 'anova_' sprintf('%s',condlevels{:}) '_fwhm0_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_r' sprintf('%d',runnum) '_' model  '_' params('smooth') 'mm.mgz'];

if optInputs(varargin, 'cluscorr')
    anova_surf = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '.' 'anova_' sprintf('%s',condlevels{:}) '_cluscorr_fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm.mgz'];
else
    anova_surf = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '.' 'anova_' sprintf('%s',condlevels{:}) '_fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm.mgz'];
end
if ~exist(anova_surf,'file') || optInputs(varargin,'overwrite')
    if optInputs(varargin, 'cluscorr')
        anova_vol = [params('rootdir') exp '/analysis/sla/usub' num2str(us) '/' runtype '_anova_' sprintf('%s',condlevels{:}) '_r' sprintf('%d',runnum) '_' model '_' num2str(fwhm_vol*100,'%.0f') 'mm.gfeat/cope1.feat/thresh_zfstat1.nii.gz'];
    else
        anova_vol = [params('rootdir') exp '/analysis/sla/usub' num2str(us) '/' runtype '_anova_' sprintf('%s',condlevels{:}) '_r' sprintf('%d',runnum) '_' model '_' num2str(fwhm_vol*100,'%.0f') 'mm.gfeat/cope1.feat/stats/zfstat1.nii.gz'];
    end
    regfile = regfile_stdLAS2anatRAS(exp, us, model, varargin{:});
    refvol = ['~/freesurfer/' subjid '/mri/orig.mgz'];
    unix_freesurfer_version(freesurfer_version, ['mri_vol2surf --ref ' refvol ' --mov ' anova_vol ' --reg ' regfile ' --o ' anova_surf ' --hemi ' hemi ' --trgsubject fsaverage --surf white --interp trilinear --projfrac-avg ' num2str(projfrac(1)) ' ' num2str(projfrac(2)) ' 0.05 '  ifelse(fwhm > 0, [' --surf-fwhm ' num2str(fwhm) ' '], '')]);
end
anova_map = MRIread(anova_surf);
anovaP = -log10(normcdf(-abs(anova_map.vol))) .* sign(anova_map.vol);
anovaP(anovaP == inf) = max(anovaP(:));
anovaP(anovaP == -inf) = min(anovaP(:));

bfmap = MRIread(bestfreq_smooth_file);
bfmap.vol(anovaP < anova_thresh) = 0;

labelstring = '';
if optInputs(varargin,'mask_label')
    label = varargin{optInputs(varargin,'mask_label')+1};
    lb = read_label(['~/freesurfer/fsaverage/label/' hemi '.' label '.label']);
    xi = setdiff(1:length(bfmap.vol), lb.vnums+1);
    bfmap.vol(xi) = 0;
    labelstring = ['_' label];
end

bestfreq_smooth_thresh_file = strrep(bestfreq_smooth_file,'.mgh',['_anovathresh' num2str(anova_thresh) labelstring '.mgh']);
MRIwrite(bfmap, bestfreq_smooth_thresh_file);

bestfreq_smooth_thresh_colorwarp_file = strrep(bestfreq_smooth_file,'.mgh','_colorwarp.mgh');
bfmap.vol(bfmap.vol > 0) = interp1(inputvals,outputvals,bfmap.vol(bfmap.vol > 0),'cubic');
MRIwrite(bfmap, bestfreq_smooth_thresh_colorwarp_file);

% plot
if optInputs(varargin,'freeview');
%     labelfile = ['~/freesurfer/fsaverage/sla/pitch_resthr_v4_us' num2str(us) '/' hemi '_harm-f0333_vs_noise-f0333_fwhm0_projfrac0-1_trilinear_localizer_block_3mm_r' sprintf('%d',runnum) '/osgm/sig_stp_v3.bin.smooth3.bin.label'];
%     freeview(bestfreq_smooth_thresh_colorwarp_file, hemi, [0.9 3.5 6.1], 'label', labelfile);
    freeview3('fsaverage',hemi,'overlay',bestfreq_smooth_thresh_colorwarp_file, 'overlay_threshold', [0.9 3.5 6.1]);
end



% % Best-frequency of vertices that exceed p-value threshold
% bestfreq = emptysurf;
% bestfreq.vol = zeros(size(bestfreq.vol));
% % f = 1:0.1:6;
% for i = 1:nverts
%     if anovaP(i) > anova_thresh;
%         %         bf = spline(1:6,betas(i,:),f);
%         %         [~,xi] = max(bf);
%         %         bestfreq.vol(i) = f(xi);
%         [~,xi] = max(betas(i,:));
%         bestfreq.vol(i) = xi;
%     end
% end

% % write to file
% idstring = ['fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_anovathresh' num2str(anova_thresh) '_' runtype '_r' sprintf('%d',runnum) '_' model  '_' params('smooth') 'mm'];
% bestfreq_file = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '.bestfreq_' idstring '.mgh'];
% bestfreq.fspec = bestfreq_file;
% MRIwrite(bestfreq, bestfreq_file);


% 
% % minimum beta for each vertex
% [~,xi] = min(betas,[],2);
% minfreq = xi;
% minfreq(all(betas==0,2)) = 0;

% % P-value associated with maxfreq > mean(all_other_freqs)
% peak_pval = nan(size(maxfreq));
% for i = 1:length(freqs)
%     
%     x = sprintf('-%d',setdiff(freqs,freqs(i)));
%     con = ['freq' num2str(freqs(i)) '_vs_freq' x(2:end)];
%     sla_surf(exp,us,con,fwhm,[],[],[],[],hemi,'localizer','block','nocorrection','fsaverage',varargin{:});
%     
%     idstring = [con '_fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_' model  '_' params('smooth') 'mm'];
%     glmdir = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '_' idstring '_r' sprintf('%d',runnum)];
%     sig = MRIread([glmdir '/osgm/sig.mgh']);
%     peak_pval(maxfreq == freqs(i)) = sig.vol(maxfreq == freqs(i));
%     
% end

% % P-value associated with all direct contrasts between two frequencies
% pval_all_freq_pairs = nan(nverts, length(freqs), length(freqs));
% for i = 1:length(freqs)
%     for j = i+1:length(freqs)
%         
%         con = ['freq' num2str(freqs(i)) '_vs_freq' num2str(freqs(j))];
%         sla_surf(exp,us,con,fwhm,[],[],[],[],hemi,'localizer','block','nocorrection','fsaverage',varargin{:});
%         idstring = [con '_fwhm' num2str(fwhm) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_' model  '_' params('smooth') 'mm'];
%         glmdir = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '_' idstring '_r' sprintf('%d',runnum)];
%         sig = MRIread([glmdir '/osgm/sig.mgh']);
%         pval_all_freq_pairs(:,i,j) = sig.vol;
%         
%     end
% end

    %     annotfile = ['/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/' hemi '.stp-pm2al-5lab-recolor.annot'];
    %     inflated_file = ['~/freesurfer/fsaverage/surf/' hemi '.inflated'];
    %     if strcmp(hemi,'rh')
    %         cam = [200 20 20];
    %     elseif strcmp(hemi,'lh')
    %         cam = [-20 -20 15];
    %     end
    %     unix_freesurfer5p2(['freeview --surface ' inflated_file ':overlay=' bestfreq_file ':overlay_threshold=' num2str(freeview_control_points(1)) ',' num2str(freeview_control_points(2)) ',' num2str(freeview_control_points(3)) ':annot=' annotfile ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' Zoom 1.2 &']);

%
%     area_file = ['~/freesurfer/fsaverage/surf/' hemi '.area2.mgh'];
% anatarg = [' -annotation ~/freesurfer/fsaverage/label/' hemi '.stp-pm2al-5lab-recolor.annot -labels-under '];
% titlestring = ['us' num2str(us) '-tonotopy'];
% if ~optInputs(varargin,'tksurfer');
%     unix_freesurfer(['tksurfer fsaverage ' hemi ' inflated -overlay ' bestfreq_file ' -fminmax 0 2.5 -title ' titlestring ' ' anatarg ' &']);
% end


% midpt = log2(1131/100);
% bounds = ;

%
% area_file = ['~/freesurfer/fsaverage/surf/' hemi '.area.mgh'];
% area_copy_file = ['~/freesurfer/fsaverage/surf/' hemi '.area2.mgh'];
% area_br = MRIread(area_file);
% area_br.fspec = area_copy_file;
% delete(area_copy_file);
% MRIwrite(area_br, area_copy_file);
% area_file = ['~/freesurfer/fsaverage/surf/' hemi '.area2.mgh'];
% inflated_file = ['~/freesurfer/fsaverage/surf/' hemi '.inflated'];
% unix(['freeview --surface ' inflated_file ':overlay=' area_copy_file ]);
%%
% 
% addpath('/software/Freesurfer/emf-5.1.0/matlab');
% addpath('/fs');