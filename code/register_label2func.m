function register_label2func(exp,roi,us,runtype,r,model,varargin)

%% directories, paramaters, etc
analysisdir = [params('rootdir') exp '/analysis/'];

func2highres = 'bbreg'; % default registration from func to highres is boundary-basedfl
if optInputs(varargin,'func2highres');
    func2highres = varargin{optInputs(varargin,'func2highres')+1}; % default registration from func to highres is boundary-based
end

fsl_version = read_fsl_version(exp,varargin{:});
fwhm = read_smooth(exp, varargin{:});

%% main

slamaskdir = [analysisdir 'sla/usub' num2str(us) '/masks/'];
featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/'];
flamaskdir = [featdir 'masks/'];

exfunc = [featdir 'reg/example_func.nii.gz'];
highres2exfunc = [featdir 'reg/highres2example_func.mat'];
if (~exist(highres2exfunc,'file') || ~exist(exfunc,'file')) && strcmp(runtype,'main_v3_combined')
    regdir = [analysisdir 'fla/usub' num2str(us) '/' strrep(runtype,'_combined','') '_r' num2str(1+(r-1)*11) '_block_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/'];
    exfunc = [regdir 'example_func.nii.gz'];
    highres2exfunc = [regdir 'highres2example_func.mat'];
end
    
anatmask_highres = [slamaskdir roi '_highres_highres2standard_freesurfer.nii.gz'];
anatmask_func = [flamaskdir roi '_func_func2highres_' func2highres '_highres2standard_freesurfer.nii.gz'];

if (~exist(anatmask_func,'file') || optInputs(varargin, 'overwrite')) && exist([flamaskdir roi '_func_func2highres_' func2highres  '.nii.gz'],'file')
    copyfile([flamaskdir roi '_func_func2highres_' func2highres  '.nii.gz'],anatmask_func);
end

if ~exist(flamaskdir,'dir'); 
    mkdir(flamaskdir);
end

if ~exist(anatmask_func,'file') || optInputs(varargin,'overwrite')
    unix_fsl(fsl_version, ['flirt -interp trilinear -in ' anatmask_highres ' -ref ' exfunc ' -applyxfm -init ' highres2exfunc ' -out ' anatmask_func]);
end