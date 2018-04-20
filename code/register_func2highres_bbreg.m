function register_func2highres_bbreg(exp,us,runtype,r,init_reg_type,varargin)

% function register_func2highres_handtune(exp, us, runtype, r, init_reg_type, varargin)
% fine-tunes an initial registration using bbregister

% FSL and freesurfer versions
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% preprocessind directory
preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];

% directory with initial registration files
init_reg_dir = [preprocdir 'reg_' init_reg_type '/'];
if ~exist(init_reg_dir, 'dir')
  mkdir(init_reg_dir);
end
 
% input functional & anatomical image
exfunc = [preprocdir 'example_func_bet.nii.gz'];
highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz'];

% freesurfer subject id
subjid = [exp '_us' num2str(us)];

% anatomical volume in freesurfer format
orig = [params('rootdir') 'freesurfer/' subjid '/mri/orig.mgz']; % orig is the anatomical volume in the freesurfer space
if ~exist(orig, 'file')
  orig = [params('rootdir') 'freesurfer/' subjid '/mri/brain.mgz']; % orig is the anatomical volume in the freesurfer space
  fprintf('Using brain.mgz instead of orig.mgz\n'); drawnow;
end

% directory with hand-tuned files
bbreg_dir = [preprocdir 'reg_bbreg/'];
if ~exist(bbreg_dir,'dir')
  mkdir(bbreg_dir);
end

exfunc2orig_init_dat = [init_reg_dir 'example_func2orig.dat'];
exfunc2orig_bbreg_dat = [bbreg_dir 'example_func2orig.dat'];

%%

% registration via bbregister
unix_freesurfer_version( freesurfer_version, ['bbregister --s ' subjid ' --mov ' exfunc ' --t2 --init-reg ' exfunc2orig_init_dat ' --reg ' exfunc2orig_bbreg_dat]);

% compare registration to original with tkregister2
tkregister_title = ['us' num2str(us) '-r' num2str(r) '-' init_reg_type];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --reg ' exfunc2orig_init_dat ' --title ' tkregister_title ' --surf &']);
tkregister_title = ['us' num2str(us) '-r' num2str(r) '-bbreg'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --reg ' exfunc2orig_bbreg_dat ' --title ' tkregister_title ' --surf &']);

% exfunc -> orig in fsl format
exfunc2orig_bbreg_mat = [bbreg_dir 'example_func2orig.mat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_bbreg_dat ' --fslregout ' exfunc2orig_bbreg_mat ' --noedit ']);

% orig -> highres, based on headers
orig2highres_mat = [bbreg_dir 'orig2highres.mat'];
orig2highres_dat = [bbreg_dir 'orig2highres.dat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' orig ' --targ ' highres ' --reg ' orig2highres_dat ' --fslregout ' orig2highres_mat ' --regheader --noedit ']);

% exfunc -> orig + orig -> highres
exfunc2highres_bbreg_mat = [bbreg_dir 'example_func2highres.mat'];
unix_fsl(fsl_version, ['convert_xfm -omat ' exfunc2highres_bbreg_mat ' -concat ' orig2highres_mat ' ' exfunc2orig_bbreg_mat]);

% exfunc2highres image
exfunc2highres_bbreg_nii = [bbreg_dir 'example_func2highres.nii.gz'];
unix_fsl(fsl_version, ['flirt -interp trilinear -in ' exfunc ' -ref ' highres ' -applyxfm -init ' exfunc2highres_bbreg_mat ' -out ' exfunc2highres_bbreg_nii]);

% view image
if optInputs(varargin, 'monkey')
    unix_freesurfer_version(freesurfer_version,['freeview ' highres ':grayscale=10,150 ' exfunc2highres_bbreg_nii ':grayscale=0,2000 &'])
else
    unix_freesurfer_version(freesurfer_version,['freeview ' highres ':grayscale=10,500 ' exfunc2highres_bbreg_nii ':grayscale=10,2000 &'])
end


