function smooth_surf_v2(exp,us,runtype,r,hemi,input_fname,varargin)

subjid = [exp '_us' num2str(us)];
freesurfer_version = read_freesurfer_version(exp,varargin{:});

[~,fwhm] = read_smooth(exp, varargin{:});
if optInputs(varargin, 'fwhm');
    fwhm = varargin{optInputs(varargin, 'fwhm')+1};
end

% map to fsaverage template
if optInputs(varargin, 'monkey')
    surfdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
else
    surfdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
end
unsmoothed_file = [surfdir hemi '.' input_fname '.mgz'];
smoothed_file = [surfdir hemi '.' input_fname '_smooth' num2str(100*fwhm, '%.0f') 'mm.mgz'];

if ~exist(smoothed_file,'file') || optInputs(varargin,'overwrite')
    if fwhm == 0
        copyfile(unsmoothed_file, smoothed_file, 'f');
    else
        if optInputs(varargin, 'monkey')
            unix_freesurfer_version(freesurfer_version, ['mris_fwhm --s ' subjid ' --i ' unsmoothed_file ' --o ' smoothed_file  ' --hemi ' hemi ' --fwhm ' num2str(fwhm) ' --smooth-only --sum ' strrep(smoothed_file,'.mgz','_reported_fwhm.txt')]);
        else
            unix_freesurfer_version(freesurfer_version, ['mri_surf2surf --s fsaverage --sval ' unsmoothed_file ' --tval ' smoothed_file  ' --hemi ' hemi ' --fwhm ' num2str(fwhm)]);
        end
    end
end