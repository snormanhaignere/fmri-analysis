function hpfilt(exp,us,runtype,r,varargin)

% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

fwhm = read_smooth(exp, varargin{:});
smooth_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm'];
intnorm_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm_intnorm'];
filt_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm_intnorm_hpfilt' num2str(params('hpcutoff'))];

% brightness threshold
mode_file = [preprocdir 'mode.txt'];
fid = fopen(mode_file,'r');
x = textscan(fid,'%f'); fclose(fid);
factor = 1e4/x{1};

fsl_version = read_fsl_version(exp);

if ~exist([filt_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite')

    % sigma
    [~, ~, TR] = read_scanparams(exp,us,runtype,'run',r,varargin{:});
    sigma = params('hpcutoff')/(2*TR);
    
    % normalize mode
    unix_fsl(fsl_version, ['fslmaths ' smooth_file ' -mul ' num2str(factor) ' '  intnorm_file]);
    
    % call filter
    unix_fsl('4.1', ['fslmaths ' intnorm_file ' -bptf ' num2str(sigma) ' -1 ' filt_file]);
    
end