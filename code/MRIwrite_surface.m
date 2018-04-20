function MRIwrite_surface(data, fname, hemi)
% function MRIwrite_surface(data, fname, hemi)
% 
% Wrapper function for the function MRIwrite, designed to handle surface data.

dims = size(data);
if length(dims) > 2
    error('Data should be two dimensional');
end

if dims(2) > dims(1);
    data = data';
    dims = size(data);
end
    
s = MRIread([params('rootdir') 'freesurfer/fsaverage/surf/' hemi '.area.mgh']);
s.vol = nan([1,dims(1),1,dims(2)]);
s.vol(1,:,1,:) = data;
s.nframes = dims(2);
s.fspec = fname;
MRIwrite(s, fname);