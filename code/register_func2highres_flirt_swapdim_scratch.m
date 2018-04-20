function register_func2highres_flirt(exp, us, runtype, r, varargin)

% function register_func2highres_flirt(exp, us, runtype, r, varargin)
% registers a functional image to a highres anatomical using FSL's
% linear registration algorithm, flirt

%% Registration using flirt

% degrees of freedom, 6 = rigid, 12 = affine
dof = '6';

% range of search
search_range = '-180 180';
if optInputs(varargin, '90')
    search_range = '-90 90';
end

costfn = 'corratio';
if optInputs(varargin, 'costfn')
    costfn = varargin{optInputs(varargin, 'costfn')+1};
end

% weight by brain extraction image
betweight = false;
if optInputs(varargin, 'betweight')
    betweight = true;
end

if optInputs(varargin, 'search_range')
    sr = varargin{optInputs(varargin, 'search_range')+1};
    search_range = ['-'  num2str(sr) ' ' num2str(sr)];
    clear sr;
end

% addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% preprocessing directory
preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
regdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/reg_flirt/'];
if ~exist(regdir, 'dir')
  mkdir(regdir);
end

if betweight
    exfunc = [preprocdir 'example_func.nii.gz']; % brain-masked functional image used for motion correction
    exfunc_swaplr = [preprocdir 'example_func_swaplr.nii.gz']; % brain-masked functional image used for motion correction
    exfunc_betmask = [preprocdir 'bet_mask.nii.gz'];
    highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz']; % brain extracted structural image
else
    exfunc = [preprocdir 'example_func_bet.nii.gz']; % brain-masked functional image used for motion correction
    exfunc_swaplr = [preprocdir 'example_func_bet_swaplr.nii.gz']; % brain-masked functional image used for motion correction
    highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz']; % brain extracted structural image
end
exfunc2highres_regmat = [regdir 'example_func2highres.mat']; % functional to anatomical registration matrix
highres2exfunc_regmat = [regdir 'highres2example_func.mat']; % anatomical to functional registration matrix
exfunc2highres_nii = [regdir 'example_func2highres.nii.gz']; % functional image interpolated to anatomical
% exfunc2highres_regmat = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/example_func2highres' dof '.mat'];
% highres2exfunc_regmat = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/highres2example_func' dof '.mat'];
% exfunc2highres_nii = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/example_func2highres' dof '.nii.gz'];

unix_fsl(fsl_version, ['fslswapdim ' exfunc ' -x y z ' exfunc_swaplr]);

exfunc = exfunc_swaplr;

if ~exist(exfunc2highres_regmat, 'file') || ~exist(highres2exfunc_regmat, 'file') || ~exist(exfunc2highres_nii, 'file') || optInputs(varargin, 'overwrite')
  % registration
  if betweight
      %' -inweight ' exfunc_betmask ' -refweight ' highres_betmask 
      unix_fsl(fsl_version, ['flirt -in ' exfunc ' -ref ' highres ' -inweight ' exfunc_betmask ' -cost ' costfn ' -dof ' dof ' -searchrx ' search_range ' -searchry ' search_range ' -searchrz ' search_range ' -omat ' exfunc2highres_regmat]);
  else
      unix_fsl(fsl_version, ['flirt -in ' exfunc ' -ref ' highres ' -cost ' costfn ' -dof ' dof ' -searchrx ' search_range ' -searchry ' search_range ' -searchrz ' search_range ' -omat ' exfunc2highres_regmat]);
  end
  % invert registration
  unix_fsl(fsl_version, ['convert_xfm -omat ' highres2exfunc_regmat ' -inverse ' exfunc2highres_regmat]);
  % apply to example_func to evaluate success
  unix_fsl(fsl_version, ['flirt -interp trilinear -in ' exfunc ' -ref ' highres ' -applyxfm -init ' exfunc2highres_regmat ' -out ' exfunc2highres_nii]);
end

%% Create freesurfer style matrix for surface resampling

% subjec id
subjid = [exp '_us' num2str(us)];

% freesurfer volume
orig = [params('rootdir') 'freesurfer/' subjid '/mri/orig.mgz']; % orig is the anatomical volume in the freesurfer format
if ~exist(orig, 'file')
  orig = [params('rootdir') 'freesurfer/' subjid '/mri/brain.mgz']; % orig is the anatomical volume in the freesurfer format
  fprintf('Using brain.mgz instead of orig.mgz\n'); drawnow;
end

% highres to orig, based on headers
highres2orig_fsl = [regdir 'highres2orig.mat'];
highres2orig_freesurf = [regdir 'highres2orig.dat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' highres ' --targ ' orig ' --reg ' highres2orig_freesurf ' --fslregout ' highres2orig_fsl ' --regheader --noedit ']);

% catenate header-based and flirt
exfunc2orig_fsl = [regdir 'example_func2orig.mat'];
unix_fsl(fsl_version, ['convert_xfm -omat ' exfunc2orig_fsl ' -concat ' highres2orig_fsl ' ' exfunc2highres_regmat]);

% create freesurfer-style matrix
exfunc2orig_freesurf = [regdir 'example_func2orig.dat'];
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_freesurf ' --fsl ' exfunc2orig_fsl ' --noedit ']);

%% Optional viewing tools

% freeview overlay
if optInputs(varargin, 'freeview')
  unix_freesurfer_version(freesurfer_version,['freeview ' highres ':grayscale=10,150 ' exfunc2highres_nii ':grayscale=0,6000 &'])
end
  
% tksurfer overlay
if optInputs(varargin, 'tkregister2')
  tmp = ['us' num2str(us) '-r' num2str(r)];
  unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_freesurf ' --title ' tmp ' --surf &']);
end