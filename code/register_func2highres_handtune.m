function register_func2highres_handtune(exp, us, runtype, r, init_reg_type, varargin)

% function register_func2highres_handtune(exp, us, runtype, r, init_reg_type, varargin)
% hand-tunes an initial registration
% with tkregister2

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
handtune_reg_dir = [preprocdir 'reg_handtune/'];
if ~exist(handtune_reg_dir,'dir')
  mkdir(handtune_reg_dir);
end

% copy over initial registration file if not already present
if ~exist([handtune_reg_dir 'example_func2orig.dat'],'file') || optInputs(varargin, 'start_from_init')
  % copy registration file to new filename already in hand-tuned directory
  % if present, prevents overwriting
  if optInputs(varargin, 'start_from_init') && exist([handtune_reg_dir 'example_func2orig.dat'],'file');
    copyfile([handtune_reg_dir 'example_func2orig.dat'], [handtune_reg_dir 'example_func2orig-' strrep(strrep(datestr(clock,0),' ','-'),':','-') '.dat']);
  end
  
  % copy over initial registration file, erasing prior file if present
  copyfile([init_reg_dir 'example_func2orig.dat'],[handtune_reg_dir 'example_func2orig.dat'],'f');
end

% perform manual registration with tkregister2
tmp = ['us' num2str(us) '-r' num2str(r)];
exfunc2orig_freesurf = [handtune_reg_dir 'example_func2orig.dat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_freesurf ' --title ' tmp ' --surf']);

%% Create auxillary registration files based on hand-tuned registration

% exfunc -> orig in fsl format
exfunc2orig_fsl = [handtune_reg_dir 'example_func2orig.mat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_freesurf ' --fslregout ' exfunc2orig_fsl ' --noedit ']);

% orig -> highres, based on headers
orig2highres_fsl = [handtune_reg_dir 'orig2highres.mat'];
orig2highres_freesurf = [handtune_reg_dir 'orig2highres.dat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' orig ' --targ ' highres ' --reg ' orig2highres_freesurf ' --fslregout ' orig2highres_fsl ' --regheader --noedit ']);

% exfunc -> orig + orig -> highres
exfunc2highres_fsl = [handtune_reg_dir 'example_func2highres.mat'];
unix_fsl(fsl_version, ['convert_xfm -omat ' exfunc2highres_fsl ' -concat ' orig2highres_fsl ' ' exfunc2orig_fsl]);

% exfunc2highres image
exfunc2highres_nii = [handtune_reg_dir 'example_func2highres.nii.gz'];
unix_fsl(fsl_version, ['flirt -interp trilinear -in ' exfunc ' -ref ' highres ' -applyxfm -init ' exfunc2highres_fsl ' -out ' exfunc2highres_nii]);

% view image
unix_freesurfer_version(freesurfer_version,['freeview ' highres ':grayscale=10,150 ' exfunc2highres_nii ':grayscale=0,2000 &'])

