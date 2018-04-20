function ipfsl(exp,us,runtype,r,varargin)

% function ipfsl(exp,us,runtype,r,varargin)
% smoothing in volume using freesurfer function mri_fwhm

freesurfer_version = read_freesurfer_version(exp);

preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

% smoothing kernel used for the experiment
% for old experiments "actual_fwhm" and "fwhm"
% differ because of an issue with formatting
% should be the same for more recent experiments
[fwhm,actual_fwhm] = read_smooth(exp, varargin{:});
if optInputs(varargin, 'fwhm');
    fwhm = varargin{optInputs(varargin, 'fwhm')+1};
    actual_fwhm = fwhm;
end

% input file is brain-extracted file
input_file = [preprocdir 'brain_thresh.nii.gz'];
if optInputs(varargin, 'input-file')
    input_file = [preprocdir varargin{optInputs(varargin, 'input-file')+1} '.nii.gz'];
end

smooth_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm.nii.gz'];
% example_func_bet = [preprocdir 'example_func_bet.nii.gz'];
mask_file = [preprocdir 'brain_mask.nii.gz'];

if ~exist(smooth_file,'file') || optInputs(varargin, 'overwrite')
    if actual_fwhm == 0
        copyfile(input_file,smooth_file);
    else
        % call susan
        unix_freesurfer_version(freesurfer_version, ['mri_fwhm --i ' input_file ' --o ' smooth_file ' --fwhm ' num2str(actual_fwhm) ' --smooth-only --mask ' mask_file]);
    end
end