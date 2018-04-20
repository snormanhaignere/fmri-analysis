function freesurfer(exp,us,varargin)

% function freesurfer(s)
% 
% runs recon-all on a subject, should only be called once per subject
% the command is run with nohup, so it is fine to exit the terminal
% after starting the reconstruction.

% directory setup
datadir = [params('rootdir') exp '/data/brain/nifti/usub' num2str(us) '/'];
freesurfdir = '~/freesurfer';
logdir = [freesurfdir '/log/'];
if ~exist(freesurfdir,'dir'); mkdir(freesurfdir); end
if ~exist(logdir,'dir'); mkdir(logdir); end
freesurfer_version = read_freesurfer_version(exp, us, varargin);

subjid = [exp '_us' num2str(us)];

structfile = [datadir 'struct_r1_LAS.nii.gz'];
logfile = [logdir 'recon-all_' subjid '.txt'];

if ~exist([freesurfdir '/' subjid '/'],'dir');
    unix_freesurfer_version(freesurfer_version, ['recon-all -i ' structfile ' -subjid ' subjid]);
end

if optInputs(varargin,'nohup')
    unix_freesurfer_version(freesurfer_version, ['nohup recon-all -all -subjid ' subjid ' > ' logfile ' &']);
elseif optInputs(varargin,'&')
    unix_freesurfer_version(freesurfer_version, ['recon-all -all -subjid ' subjid ' > ' logfile ' &']);
else
    unix_freesurfer_version(freesurfer_version, ['recon-all -all -subjid ' subjid ' > ' logfile]);
end