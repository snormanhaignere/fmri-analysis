function register_func2func_flirt(exp, us, runtype, r1, r2, varargin)

% Registers a functional image from run r1 to run r2. 
% 
% 2016-02-16: Created, Sam NH
% 
% 2017-11-14: Fixed a bug (number not converted to a string in a file name)

%% Registration using flirt

% combine two runs into 2-element array
assert(isscalar(r1) && isscalar(r2));
runs = [r1, r2];
% clear r1 r2;

% degrees of freedom, 6 = rigid, 12 = affine
dof = '6';

% cost function
costfn = 'corratio';
if optInputs(varargin, 'costfn')
    costfn = varargin{optInputs(varargin, 'costfn')+1};
end

% range of search
search_range = '-180 180';
if optInputs(varargin, '90')
    search_range = '-90 90';
end

% weight by brain extraction image
betweight = false;
if optInputs(varargin, 'betweight')
    betweight = true;
end

% addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% preprocessing directory and functional for the target and reference
preprocdir = cell(1,2);
exfunc = cell(1,2);
betmask = cell(1,2);
for i = 1:2
    % preprocessing directory
    preprocdir{i} = [params('rootdir') exp '/analysis/preprocess' ...
        '/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/'];
    
    % only use brain-extracted functional image if not weighting by that image
    if betweight
        exfunc{i} = [preprocdir{i} 'example_func.nii.gz'];
        betmask{i} = [preprocdir{i} 'bet_mask.nii.gz'];
    else
        exfunc{i} = [preprocdir{i} 'example_func_bet.nii.gz'];
    end
end

% directory to save registration files
regdir = [preprocdir{1} 'reg_flirt/'];
if ~exist(regdir, 'dir')
    mkdir(regdir);
end

% registration matrix
reg_mat = [regdir 'example_func2example_func_r' num2str(runs(2)) '.mat'];

% functional image interpolated to r2
registered_image = [regdir 'example_func2example_func_r' num2str(runs(2)) '.nii.gz']; 

if ~exist(reg_mat, 'file') || ...
        ~exist(registered_image, 'file') || optInputs(varargin, 'overwrite')
    
    command = ['flirt -in ' exfunc{1} ' -ref ' exfunc{2} ...
        ' -cost ' costfn ' -dof ' dof ' -searchrx ' search_range ...
        ' -searchry ' search_range ' -searchrz ' search_range ...
        ' -omat ' reg_mat];
    
    % add weights
    if betweight
        command = [command ' -inweight ' betmask{1} ' -refweight ' betmask{2}];
    end
    
    % registration
    unix_fsl(fsl_version, command);
    
    % apply to example_func to evaluate success
    unix_fsl(fsl_version, ...
        ['flirt -interp trilinear -in ' exfunc{1} ' -ref ' exfunc{2} ...
        ' -applyxfm -init ' reg_mat ' -out ' registered_image]);
end

% create freesurfer-style matrix
subjid = [exp '_us' num2str(us)];
reg_mat_freesurf = strrep(reg_mat, '.mat', '.dat');
unix_freesurfer_version( freesurfer_version, ...
    ['tkregister2 --s ' subjid ' --mov ' exfunc{1} ' --targ ' exfunc{2} ...
    ' --reg ' reg_mat_freesurf ' --fsl ' reg_mat ' --noedit ']);

%% Optional viewing tools

% freeview overlay
if optInputs(varargin, 'freeview')
  unix_freesurfer_version(freesurfer_version,...
      ['freeview ' exfunc{2} ':grayscale=0,6000 ' ...
      registered_image ':grayscale=0,6000 &'])
end