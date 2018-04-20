function register_func2highres_apply(exp, us, runtype, r, reg_type, fname, varargin)

% function register_func2highres_handtune(exp, us, runtype, r, init_reg_type, varargin)
% applyies registration to functional volume, resampling to downsampled structural

% FSL and freesurfer versions
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% preprocessind directory
preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];

% directory with initial registration files
regdir = [preprocdir 'reg_' reg_type '/'];
if ~exist(regdir, 'dir')
  mkdir(regdir);
end
 
% downsampled structural
highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz'];
highres_2mm = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain_2mm.nii.gz'];
if ~exist(highres_2mm,'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ['fslmaths ' highres ' -subsamp2 ' highres_2mm]);
end

% input functional & anatomical image
func = [preprocdir fname '.nii.gz'];
func_highres_2mm = [preprocdir fname '_highres_2mm.nii.gz'];
exfunc2highres = [regdir 'example_func2highres.mat'];

% exfunc2highres image
if ~exist(func_highres_2mm,'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ['flirt -interp trilinear -in ' func ' -ref ' highres_2mm ' -applyxfm -init ' exfunc2highres ' -out ' func_highres_2mm]);
end

% optionally view image
if optInputs(varargin, 'freeview')
    unix_freesurfer_version(freesurfer_version,['freeview ' highres ' ' func_highres_2mm '&']);
end

