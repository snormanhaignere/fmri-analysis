function register_whitematter_highres_2mm(exp,us,varargin)

analysisdir = [params('rootdir') exp '/analysis/'];
subjid = [exp '_us' num2str(us)];
freesurferdir = [params('rootdir') 'freesurfer/' subjid '/'];
fsl_version = read_fsl_version(exp, varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

aseg_fsurf_mgz = [freesurferdir 'mri/aseg.mgz'];
wm_fsurf_mgz = [freesurferdir 'mri/wm_bin.mgz'];
if ~exist(wm_fsurf_mgz,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ['mri_binarize --i ' aseg_fsurf_mgz ' --o ' wm_fsurf_mgz ' --match 41 2']);
end

highres_fsl_niigz = [analysisdir 'preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz'];
highres_fsl_2mm_niigz = [analysisdir 'preprocess/usub' num2str(us) '/struct_r1/brain_2mm.nii.gz'];
if ~exist(highres_fsl_2mm_niigz,'file')
    unix_fsl(fsl_version, ['fslmaths ' highres_fsl_niigz ' -subsamp2 ' highres_fsl_2mm_niigz]);
end

highres_fsurf_mgz = [freesurferdir 'mri/orig.mgz'];
highres_fsurf_to_highres_fsl_dat = [freesurferdir 'mri/transforms/highres_fsurf_to_highres_2mm_fsl.dat'];
% highres_fsurf_to_highres_fsl_dat = [freesurferdir 'mri/transforms/highres_fsurf_to_highres_fsl.dat'];
% unix_freesurfer(['tkregister2 --s ' subjid ' --mov ' highres_fsurf_mgz ' --targ ' highres_fsl_niigz ' --reg ' highres_fsurf_to_highres_fsl_dat ' --regheader --noedit ']);
if ~exist(highres_fsurf_to_highres_fsl_dat,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' highres_fsurf_mgz ' --targ ' highres_fsl_2mm_niigz ' --reg ' highres_fsurf_to_highres_fsl_dat ' --regheader --noedit ']);
end

wm_fsl_niigz = [analysisdir 'preprocess/usub' num2str(us) '/struct_r1/white_matter_2mm_highres.nii.gz'];
if ~exist(wm_fsl_niigz,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ['mri_vol2vol --mov ' wm_fsurf_mgz ' --targ ' highres_fsl_2mm_niigz ' --reg ' highres_fsurf_to_highres_fsl_dat ' --o ' wm_fsl_niigz ' --interp trilin ']);
end

% unix_freesurfer(['freeview ' highres_fsl_2mm_niigz ' ' wm_fsl_niigz ' ']);