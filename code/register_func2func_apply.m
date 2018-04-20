function register_func2func_apply(exp, us, runtype, r1, r2, reg_type, fname, varargin)

% Applies a functional to functional registration to an input image

% combine two runs into 2-element array
assert(isscalar(r1) && isscalar(r2));
runs = [r1, r2];
clear r1 r2;

% FSL and freesurfer versions
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% preprocessing directory for the target and reference
preprocdir = cell(1,2);
for i = 1:2
    % preprocessing directory
    preprocdir{i} = [params('rootdir') exp '/analysis/preprocess' ...
        '/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/'];    
end

% directory with initial registration files
regmat = [preprocdir{1} 'reg_' reg_type '/' ...
    'example_func2example_func_r' num2str(runs(2)) '.mat'];

% reference functional
func_ref = [preprocdir{2} 'example_func_bet.nii.gz']; 

% input functional & registered image
func = [preprocdir{1} fname '.nii.gz'];
func_registered = [preprocdir{1} fname '_reg_to_r' num2str(runs(2)) '.nii.gz'];

% exfunc2highres image
if ~exist(func_registered,'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ...
        ['flirt -interp trilinear -in ' func ' -ref ' func_ref ...
        ' -applyxfm -init ' regmat ' -out ' func_registered]);
end

% optionally view image
if optInputs(varargin, 'freeview')
    unix_freesurfer_version(freesurfer_version,...
        ['freeview ' func_ref ' ' func_registered '&']);
end

