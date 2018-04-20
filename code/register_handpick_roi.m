function register_handpick_roi(exp,roi,us,runtype,r,model,varargin)

analysisdir = [params('rootdir') exp '/analysis/'];

func2highres = 'bbreg'; % default registration from func to highres is boundary-based
if optInputs(varargin,'func2highres');
    func2highres = varargin{optInputs(varargin,'func2highres')+1}; % default registration from func to highres is boundary-based
end

fsl_version = read_fsl_version(exp,varargin{:});
fwhm = read_smooth(exp, varargin{:});

if ~isempty(strfind(roi, '2mm'))
    resolution = '2mm';
else
    resolution = '1mm';
end

%% main

roi_directory = [params('rootdir') 'handpicked_rois/us' num2str(us) '/'];
featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/'];
flamaskdir = [featdir 'masks/'];

exfunc = [featdir 'reg/example_func.nii.gz'];
highres2exfunc = [featdir 'reg/highres2example_func.mat'];

anatmask_highres = [roi_directory roi '.nii.gz'];
anatmask_func = [flamaskdir roi '_func_func2highres_' func2highres '_highres2standard_flirt.nii.gz'];

if (~exist(anatmask_func,'file') || optInputs(varargin, 'overwrite')) && exist([flamaskdir roi '_func_func2highres_' func2highres  '.nii.gz'],'file')
    copyfile([flamaskdir roi '_func_func2highres_' func2highres  '.nii.gz'],anatmask_func);
end

if ~exist(flamaskdir,'dir'); 
    mkdir(flamaskdir);
end

if ~exist(anatmask_func,'file') || optInputs(varargin,'overwrite')
    unix_fsl(fsl_version, ['flirt -interp trilinear -in ' anatmask_highres ' -ref ' exfunc ' -applyxfm -init ' highres2exfunc ' -out ' anatmask_func]);
end