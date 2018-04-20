function check_bet(exp, us, runtype, runs, varargin)

% function check_bet(exp, us, runtype, runs, varargin)
% useful function for check success of brain extraction algorithm

freesurfer_version = read_freesurfer_version(exp,varargin{:});
for i = 1:length(runs)
    preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/'];
    unix_freesurfer_version(freesurfer_version, ['freeview ' preprocdir 'example_func.nii.gz:grayscale=0,2000 ' preprocdir 'brain_mask.nii.gz:colormap=heat:opacity=0.5 &']);
end