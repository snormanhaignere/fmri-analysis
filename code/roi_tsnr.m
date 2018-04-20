function tsnr_mean = roi_tsnr(exp, us, roi, runtype, type, varargin)

% [psc testruns_used analysis_string] = roi_psc(s, testruntype, testmodel, masktype, roi, varargin)
% For anatomical + contrast: roi_psc(s, testruntype, testmodel, masktype, roi, locruntype, locmodel, contrast)
%
% Main function for computing the psc value for an roi. If localizer and test runs are overlapping, uses leave-one-out style analysis.
% Specificy a subject, a test runtype (e.g. 'main', 'localizer'), a test model (e.g. 'block' or 'event'), and then a type of mask: anat, anat_contrast, and cluster
%
% For anat you only need to specify an roi (needs to be updated)
% For cluster and anat_contrast, you need to specify an anatomical roi, a localizer contrast, a localizer runtype, & a localizer model (next two arguments respectively)
%
% Sample
% [psc testruns_used analysis_string] = roi_psc(1, 'main', 'block', 'anat_contrast', 'destrieux_pp', 'harm-highfreq_vs_noise-highfreq', 'localizer', 'block')
% roi_psc(1, 'main', 'block', 'anat', 'hvprob_pp', 'harm_vs_noise', 'main', 'block')

%% Directories
analysisdir = [params('rootdir') exp '/analysis/'];
tsnrdir = [analysisdir 'tsnr/usub' num2str(us) '/'];
roidir = [analysisdir 'roi/'];
if ~exist(roidir,'dir')
    mkdir(roidir,'dir');
end
model = 'block';

%% Registration Algorithm Used

% types of registration to move between
% functional, highres, and standard spaces.
% only used for registering anatomical masks.

% defaults are bbregister to go from functional
% to highres (flirt is alternative)
% and flirt to go from highres to standard (fnirt is alternative).

% if freesurfer roi is being used
% then there is no highres to standard transformation
% and so this is just set to 'freesurfer'

gm_mask = false;
if optInputs(varargin, 'mask')
    gm_mask = true;
end

highres2standard = 'flirt';
if optInputs(varargin,'highres2standard');
    highres2standard = varargin{optInputs(varargin,'highres2standard')+1};
end

if ~isempty(findstr('destrieux_',roi)) || ~isempty(findstr('mylabel_',roi)) % destrieux is from freesurfer and so doesnt really have a standard to highres transformation
    highres2standard = 'freesurfer';
end

% string used to remember this specific analysis
analysis_string = ['tsnr_us' num2str(us) '_' roi '_' runtype '_' type '_' model '_highres2standard_' highres2standard '_' params('smooth') 'mm' '_' DataHash(varargin)];

if optInputs(varargin, 'allruns')
    matfile = [roidir analysis_string ifelse(~gm_mask,'_nomask','') '_allruns.mat'];
else
    matfile = [roidir analysis_string ifelse(~gm_mask,'_nomask','') '.mat'];
end

if exist(matfile,'file') && ~optInputs(varargin,'overwrite')
    load(matfile);
    return;
end


%% Calculate Mask

if strcmp(type, 'raw')
    modelstr = '';
elseif strcmp(type, 'preproc')
    modelstr = ['_' model '_' params('smooth') 'mm'];
end

% roi
anatfile = [analysisdir 'sla/usub' num2str(us) '/masks/' roi  '_highres_highres2standard_' highres2standard  '.nii.gz'];
anatbr = readmr(anatfile,'NOPROGRESSBAR');

% grey matter mask
if gm_mask
    gmfile = [tsnrdir 'greymatter.nii.gz'];
    gmbr = readmr(gmfile,'NOPROGRESSBAR');
    maskbr = anatbr;
    maskbr.data = anatbr.data .* gmbr.data;
else
    maskbr = anatbr;
end

if any(maskbr.data(:) < 0) 
    fprintf('Mask should not have negative values\n');
    drawnow;
    keyboard;
end

sum(maskbr.data(maskbr.data>0.1))/sum(maskbr.data(:))

if optInputs(varargin, 'allruns')
    allruns = read_runs(exp,us,runtype,varargin{:});
    tsnr_mean = nan(length(allruns),1);
    for i = 1:length(allruns)
        switch type
            case 'raw'
                tsnrbr = readmr([tsnrdir 'tsnr_' type '_' runtype '_r' num2str(allruns(i)) '_highres.nii.gz'],'NOPROGRESSBAR');
            case 'preproc'
                tsnrbr = readmr([tsnrdir 'tsnr_' type '_' runtype '_r' num2str(allruns(i)) '_' model '_' params('smooth') 'mm' '_highres.nii.gz'],'NOPROGRESSBAR');
        end
        x = maskbr.data .* tsnrbr.data;
        y = maskbr.data;
        tsnr_mean(i) = sum(x(:))/sum(y(:));
        fprintf('Run %d: %.2f\n',allruns(i),tsnr_mean(i)); drawnow;
    end
else
    
    tsnrfile = [tsnrdir 'tsnr_' type '_mean_' runtype modelstr '_highres.nii.gz'];
    
    tsnrbr = readmr(tsnrfile,'NOPROGRESSBAR');
    x = maskbr.data .* tsnrbr.data;
    y = maskbr.data;
    tsnr_mean = sum(x(:))/sum(y(:));
end

save(matfile,'tsnr_mean');

