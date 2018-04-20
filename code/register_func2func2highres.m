function register_func2func2highres(exp, us, runtype, r1, r2, func2anat_regtype, varargin)

% Registers a functional to an anatomical by registering the functional to
% another functional and registering that functional to an anatomical
% func r1 -> func r2 -> highres anatomical


%% Registration using flirt

% combine two runs into 2-element array
assert(isscalar(r1) && isscalar(r2));
runs = [r1, r2];
assert(r1 ~= r2)

% addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% preprocessing directory for the target and reference
preprocdir = cell(1,2);
for i = 1:2
    preprocdir{i} = [params('rootdir') exp '/analysis/preprocess' ...
        '/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/'];
end

% create the output registration directory if not present
if ~exist([preprocdir{1} 'reg_' func2anat_regtype], 'dir')
    mkdir([preprocdir{1} 'reg_' func2anat_regtype]);
end

exfunc1 = [preprocdir{1} 'example_func.nii.gz'];
reg_func1_to_func2_mat = [preprocdir{1} 'reg_flirt/example_func2example_func_r' num2str(r2) '.mat']; 
reg_func2_to_highres_mat = [preprocdir{2} 'reg_' func2anat_regtype '/example_func2highres.mat']; 
reg_func1_to_highres_mat = [preprocdir{1} 'reg_' func2anat_regtype '/example_func2highres.mat']; 
reg_func2_to_orig_mat = [preprocdir{2} 'reg_' func2anat_regtype '/example_func2orig.mat']; 
reg_func1_to_orig_mat = [preprocdir{1} 'reg_' func2anat_regtype '/example_func2orig.mat']; 

% func1 -> func2 -> anat (highres and orig)
unix_fsl(fsl_version, ['convert_xfm -omat ' reg_func1_to_highres_mat ' -concat ' ...
    reg_func2_to_highres_mat ' ' reg_func1_to_func2_mat]);
unix_fsl(fsl_version, ['convert_xfm -omat ' reg_func1_to_orig_mat ' -concat ' ...
    reg_func2_to_orig_mat ' ' reg_func1_to_func2_mat]);

% create freesurfer-style matrix
subjid = [exp '_us' num2str(us)];
reg_func1_to_highres_dat = strrep(reg_func1_to_highres_mat, '.mat', '.dat');
reg_func1_to_orig_dat = strrep(reg_func1_to_orig_mat, '.mat', '.dat');
highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz'];
orig = [params('rootdir') 'freesurfer/' subjid '/mri/orig.mgz']; % orig is the anatomical volume in the freesurfer format
if ~exist(orig, 'file')
  orig = [params('rootdir') 'freesurfer/' subjid '/mri/brain.mgz']; % orig is the anatomical volume in the freesurfer format
  fprintf('Using brain.mgz instead of orig.mgz\n'); drawnow;
end
unix_freesurfer_version( freesurfer_version, ...
    ['tkregister2 --s ' subjid ' --mov ' exfunc1 ' --targ ' highres ...
    ' --reg ' reg_func1_to_highres_dat ' --fsl ' reg_func1_to_highres_mat ' --noedit ']);
unix_freesurfer_version( freesurfer_version, ...
    ['tkregister2 --s ' subjid ' --mov ' exfunc1 ' --targ ' orig ...
    ' --reg ' reg_func1_to_orig_dat ' --fsl ' reg_func1_to_orig_mat ' --noedit ']);

% view
if optInputs(varargin, 'tkregister2')
    tmp = ['us' num2str(us) '-r' num2str(r1) '-to-anat'];
    unix_freesurfer_version( freesurfer_version, ...
        ['tkregister2 --s ' subjid ' --mov ' exfunc1 ' --targ ' orig ...
        ' --reg ' reg_func1_to_orig_dat ' --title ' tmp ' --surf']);
end


