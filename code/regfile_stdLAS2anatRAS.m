function standardLAS2highresRAS = regfile_stdLAS2anatRAS(exp, us, model, varargin)

freesurfer_version = read_freesurfer_version(exp,varargin{:});
subjid = [exp '_us' num2str(us)];
runtypes = read_runtypes(exp, us, varargin{:});
runs = read_runs(exp,us,runtypes{1},varargin{:});
fwhm = read_smooth(exp, varargin{:});

sample_feat_dir = [params('rootdir') exp '/analysis/fla/usub' num2str(us) '/' runtypes{1} '_r' num2str(runs(1)) '_' model '_' num2str(fwhm*100,'%.0f') 'mm.feat/'];

standardLAS = [sample_feat_dir 'reg/standard.nii.gz'];
highresLAS = [sample_feat_dir 'reg/highres.nii.gz'];
highresRAS = ['~/freesurfer/' subjid '/mri/orig.mgz'];

standardLAS2highresLAS_fsl = [sample_feat_dir 'reg/standard2highres.mat'];
standardLAS2highresLAS_freesurf = [sample_feat_dir 'reg/standard2highres.dat'];
highresLAS2highresRAS = ['~/freesurfer/' subjid '/mri/transforms/highresLAS2highresRAS.dat'];
standardLAS2highresRAS = ['~/freesurfer/' subjid '/mri/transforms/standardLAS2highresRAS.dat'];

if ~exist(standardLAS2highresRAS, 'file') || optInputs(varargin, 'overwrite')
    unix_freesurfer_version(freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' standardLAS ' --targ ' highresLAS ' --fsl ' standardLAS2highresLAS_fsl ' --reg ' standardLAS2highresLAS_freesurf ' --noedit ']);
    unix_freesurfer_version(freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' highresLAS ' --targ ' highresRAS ' --reg ' highresLAS2highresRAS ' --noedit --regheader ']);
    unix_freesurfer_version(freesurfer_version, ['mri_matrix_multiply -im ' standardLAS2highresLAS_freesurf ' -im ' highresLAS2highresRAS ' -om ' standardLAS2highresRAS]);
end