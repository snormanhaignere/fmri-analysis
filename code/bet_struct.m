function bet_struct(exp,us,varargin)

% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

% fa values
switch exp
    case 'naturalsound-nmf-localizer'
        favalue = 0.3;
    case {'tono-pitch-localizer'}
        favalue = 0.3;
    case {'tono-localizer'}
        favalue = 0.5;
    otherwise
        favalue = 0.15;
end

if optInputs(varargin, 'favalue')
    favalue = varargin{optInputs(varargin, 'favalue')+1};
end

datadir = [params('rootdir') exp '/data/'];
niidir = [datadir 'brain/nifti/usub' num2str(us) '/'];

preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

las_file = [niidir 'struct_r1_LAS'];
raw_file = [preprocdir 'raw'];
bet_file = [preprocdir 'brain'];
fsl_version = read_fsl_version(exp);

if ~exist([raw_file '.nii.gz'], 'file') || optInputs(varargin, 'overwrite');
    copyfile([las_file '.nii.gz'],[raw_file '.nii.gz'],'f');
end

if ~exist([bet_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version,['bet2 ' raw_file ' ' bet_file ' -f ' num2str(favalue) ' -m']);
end