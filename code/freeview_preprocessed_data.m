function freeview_preprocessed_data(exp, us, runtype, runs, input_fname, varargin)

% function freeview_preprocessed_data(exp, us, runtype, runs, varargin)
% wrapper function for freeview
% loads functional files from the directory:
% svnh/EXP/analysis/preprocess/usubUS/RUNTYPE_rRUNS(i)/

freesurfer_version = read_freesurfer_version(exp,varargin{:});
for i = 1:length(runs)
    fname = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/' input_fname '.nii.gz'];
    unix_freesurfer_version(freesurfer_version,['freeview ' fname ':grayscale=0,2000 &']);
end