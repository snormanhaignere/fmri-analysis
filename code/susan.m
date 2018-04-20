function susan(exp,us,runtype,r,varargin)

% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

[fwhm,actual_fwhm] = read_smooth(exp, varargin{:});
mean_func_file = [preprocdir 'mean_func'];
thresh_file = [preprocdir 'brain_thresh'];
mode_file = [preprocdir 'mode.txt'];
smooth_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm'];

fsl_version = read_fsl_version(exp);

if ~exist([smooth_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite')
    if actual_fwhm == 0
        copyfile([thresh_file '.nii.gz'],[smooth_file '.nii.gz']);
    else
        % brightness threshold
        fid = fopen(mode_file,'r');
        x = textscan(fid,'%f'); fclose(fid);
        brightness_thresh = x{1}*0.75;
        
        % sigma
        sigma = actual_fwhm/(2*sqrt(2*log(2)));
        
        % call susan
        unix_fsl(fsl_version, ['susan ' thresh_file ' ' num2str(brightness_thresh) ' ' num2str(sigma) ' 3 0 1 ' mean_func_file ' ' num2str(brightness_thresh) ' ' smooth_file]);
    end
end