function funcLAS2highresRAS = regfile_funcLAS2anatRAS(exp, us, runtype, r, model, varargin)

freesurfer_version = read_freesurfer_version(exp,varargin{:});
fwhm = read_smooth(exp, varargin{:});
analysisdir = [params('rootdir') exp '/analysis/'];
subjid = [exp '_us' num2str(us)];

funcLAS = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/example_func.nii.gz'];
highresLAS = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/highres.nii.gz'];
highresRAS = ['~/freesurfer/' subjid '/mri/orig.mgz'];

funcLAS2highresLAS_fsl = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/example_func2highres.mat'];
funcLAS2highresLAS_freesurf = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/example_func2highresLAS.dat'];
highresLAS2highresRAS = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/highresLAS2highresRAS.dat'];
funcLAS2highresRAS = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/example_func2highresRAS.dat'];

if ~exist(funcLAS2highresRAS, 'file') || optInputs(varargin, 'overwrite')
    unix_freesurfer_version(freesurfer_version,['tkregister2 --s ' subjid ' --mov ' funcLAS ' --targ ' highresLAS ' --fsl ' funcLAS2highresLAS_fsl ' --reg ' funcLAS2highresLAS_freesurf ' --noedit ']);
    unix_freesurfer_version(freesurfer_version,['tkregister2 --s ' subjid ' --mov ' highresLAS ' --targ ' highresRAS ' --reg ' highresLAS2highresRAS ' --noedit --regheader ']);
    unix_freesurfer_version(freesurfer_version,['mri_matrix_multiply -im ' funcLAS2highresLAS_freesurf ' -im ' highresLAS2highresRAS ' -om ' funcLAS2highresRAS]);
end