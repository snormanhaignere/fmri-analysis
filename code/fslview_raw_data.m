function fslview_raw_data(exp, us, runtype, runs, varargin)

freesurfer_version = read_freesurfer_version(exp,varargin{:});
for i = 1:length(runs)
    fname = [params('rootdir') exp '/data/brain/nifti/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '_falseheader_reorient.nii.gz'];
    if strcmp(computer,'MACI64')
        [~,y] = unix(['/usr/local/fsl/bin/fslview.app/Contents/MacOS/fslview ' fname ' &']);
    else
        unix_freesurfer_version(freesurfer_version,['freeview ' fname ' &'])
    end
end