function resample_func2surface(exp,us,runtype,r,hemi,volname,regtype,varargin)

% if us == 170
%     projfrac = [-0.5 1.5];
% end
projfrac = [0 1];

interpmethod = 'trilinear';
subjid = [exp '_us' num2str(us)];
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% volume file to resample
voldir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(voldir,'dir');
    mkdir(voldir);
end
volfile = [voldir volname '.nii.gz'];

% surface_file
surfdir = ['~/freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(surfdir,'dir');
    mkdir(surfdir);
end
surf_file = [surfdir hemi '.' volname '.mgz'];

% resample to surface
regfile = [voldir 'reg_' regtype '/example_func2orig.dat'];
refvol = ['~/freesurfer/' subjid '/mri/orig.mgz'];
if ~exist(surf_file,'file') || optInputs(varargin,'overwrite')
  unix_freesurfer_version(freesurfer_version, ['mri_vol2surf --ref ' refvol ' --mov ' volfile ' --reg ' regfile ' --o ' surf_file ' --hemi ' hemi ' --surf white --interp ' interpmethod ' --projfrac-avg ' num2str(projfrac(1)) ' ' num2str(projfrac(2)) ' 0.05']);
end

% don't use fsaverage for monkey brain
if optInputs(varargin, 'monkey')
    return;
end

% map to fsaverage template
fsaverage_surfdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
if ~exist(fsaverage_surfdir,'dir')
    mkdir(fsaverage_surfdir)
end

fsaverage_file = [fsaverage_surfdir hemi '.' volname '.mgz'];
if ~exist(fsaverage_file,'file') || optInputs(varargin,'overwrite')
  unix_freesurfer_version(freesurfer_version, ['mri_surf2surf --srcsubject ' subjid ' --sval ' surf_file ' --trgsubject fsaverage --tval ' fsaverage_file  ' --hemi ' hemi]);
end