function custom_feat_reg_highres_2mm(exp, us, runtype, r, varargin)

freesurfer_version = read_freesurfer_version(exp, varargin{:});
fwhm = read_smooth(exp, us, runtype, varargin{:});
model = 'block';
feat_directory = [params('rootdir') exp '/analysis/fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model  '_' num2str(fwhm*100,'%.0f') 'mm.feat/'];
reg_directory = [feat_directory 'reg/'];

if exist(reg_directory,'dir') && ~optInputs(varargin, 'overwrite')
    return;
end

if ~exist(reg_directory,'dir')
    mkdir(reg_directory);
end

subjid = [exp '_us' num2str(us)];

% highres 2mm to highres
highres_2mm = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain_2mm.nii.gz'];
highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz'];
highres_2mm_to_highres_mat = [reg_directory 'example_func2highres.mat'];
highres_2mm_to_highres_dat = [reg_directory 'example_func2highres.dat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' highres_2mm ' --targ ' highres ' --reg ' highres_2mm_to_highres_dat ' --fslregout ' highres_2mm_to_highres_mat ' --regheader --noedit ']);
delete(highres_2mm_to_highres_dat);

% identity matrix for highres to standard
highres_to_highres_mat = [reg_directory 'highres2standard.mat'];
highres_to_highres_dat = [reg_directory 'highres2standard.dat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' highres ' --targ ' highres ' --reg ' highres_to_highres_dat ' --fslregout ' highres_to_highres_mat ' --regheader --noedit ']);
delete(highres_to_highres_dat);

copyfile(highres,[reg_directory 'highres.nii.gz'],'f')
copyfile(highres,[reg_directory 'standard.nii.gz'],'f')
