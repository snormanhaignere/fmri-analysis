function register_func2func_handtune(exp, us, runtype, r1, r2, varargin)

% Fine-tune a registration between two functional images by hand
% 
% 2017-11-14: Created, Sam NH

% combine two runs into 2-element array
assert(isscalar(r1) && isscalar(r2));
runs = [r1, r2];

% FSL and freesurfer versions
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% preprocessing directory and functional for the target and reference
preprocdir = cell(1,2);
exfunc = cell(1,2);
for i = 1:2
    % preprocessing directory
    preprocdir{i} = [params('rootdir') exp '/analysis/preprocess' ...
        '/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/'];
    % brain-extracted functional image to register
    exfunc{i} = [preprocdir{i} 'example_func_bet.nii.gz']; 
end

% directory to save registration files
init_regdir = [preprocdir{1} 'reg_flirt/'];
init_regmat_freesurf = [init_regdir 'example_func2example_func_r' num2str(runs(2)) '.dat'];

% directory with hand-tuned files
handtune_regdir = [preprocdir{1} 'reg_handtune/'];
handtune_regmat_freesurf = [handtune_regdir 'example_func2example_func_r' num2str(runs(2)) '.dat'];
if ~exist(handtune_regdir,'dir')
  mkdir(handtune_regdir);
end

% copy over initial registration file if not already present
if ~exist(handtune_regmat_freesurf,'file') || optInputs(varargin, 'start_from_init')
  % copy registration file to new filename already in hand-tuned directory
  % if present, prevents overwriting
  if optInputs(varargin, 'start_from_init') && exist(handtune_regmat_freesurf,'file')
      date_string = strrep(strrep(datestr(clock,0),' ','-'),':','-');
      dated_file = strrep(handtune_regmat_freesurf, '.dat', ['-' date_string '.dat']);
      copyfile(handtune_regmat_freesurf, dated_file);
  end
  
  % copy over initial registration file, erasing prior file if present
  copyfile(init_regmat_freesurf, handtune_regmat_freesurf, 'f');
end

% freesurfer subject id
subjid = [exp '_us' num2str(us)];

% perform manual registration with tkregister2
tmp = ['us' num2str(us) '-r' num2str(runs(1)) '-to-r' num2str(runs(2))];
command = ['tkregister2 --s ' subjid  ' --mov ' exfunc{1} ...
    ' --targ ' exfunc{2} ' --reg ' handtune_regmat_freesurf ...
    ' --title ' tmp ' --surf'];
unix_freesurfer_version( freesurfer_version, command);
