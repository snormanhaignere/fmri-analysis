function register_label2highres(exp,us,label,hemi,varargin)

% function register_labels_highres(s,label,labelnum,varargin)
%
% transforms freesurfer labels from the surface to the highres volume
% uses mri_annotation2label and mri_label2vol

%% paramaters, directories, files, etc.

analysisdir = [params('rootdir') exp '/analysis/'];
datadir = [params('rootdir') exp '/data/brain/nifti/usub' num2str(us) '/' ];
voldir = [analysisdir 'sla/usub' num2str(us) '/masks/'];
subjid = [exp '_us' num2str(us)];
% fwhm = read_smooth(exp, varargin{:});

fillthresh = 0;
start = 0;
stop = 1;
delta = 0.01;

krn = 2/2.355; % kernel size in mm
cutoff = 0.5; % cutoff after smoothing

orig = ['~/freesurfer/' subjid '/mri/orig.mgz'];
% orig = ['~/pitch_dp/analysis/freesurfer/sub1/mri/orig.mgz'];

struct = [analysisdir 'preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz'];
if ~exist(struct,'file')
    struct = [datadir 'struct_r1_LAS.nii.gz'];
end

if ~exist(voldir,'dir')
    mkdir(voldir);
end

fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

%%

if optInputs(varargin,'subspec')
    inputlabel = ['~/freesurfer/fsaverage/label/us' num2str(us) '/' hemi '.' label '.label'];
else
    inputlabel = ['~/freesurfer/fsaverage/label/' hemi '.' label '.label'];
end

outputlabel = ['~/freesurfer/' subjid '/label/' hemi '.' label '.label' ];
if ~exist(outputlabel,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ['mri_label2label --srclabel ' inputlabel ' --srcsubject fsaverage --trglabel ' outputlabel ' --trgsubject ' subjid ' --regmethod surface --hemi ' hemi]);
end

%%

labelfile = ['~/freesurfer/' subjid '/label/' hemi '.' label '.label' ];
volfile = [voldir 'mylabel_' hemi '-' label '_highres_highres2standard_freesurfer.nii.gz'];

oldvol = [voldir 'mylabel_' hemi '-' label '.nii.gz'];
if (~exist(volfile,'file') || optInputs(varargin,'overwrite')) && exist(oldvol,'file');
    copyfile(oldvol,volfile);
end    

if ~exist(volfile,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ['mri_label2vol --subject ' subjid ' --label ' labelfile ' --hemi ' hemi ' --fillthresh ' num2str(fillthresh) ' --proj frac ' num2str(start) ' ' num2str(stop) ' ' num2str(delta) ' --temp ' struct ' --regheader ' orig ' --o ' volfile]);
    unix_fsl(fsl_version, ['fslmaths ' volfile ' -kernel gauss ' num2str(krn) ' -fmean ' volfile]);%
end