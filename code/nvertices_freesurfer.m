function nvertices = nvertices_freesurfer(subjid, us, hemi)

% read vertices and faces
if strcmp(subjid, 'fsaverage')
    inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated'];
else
    switch us
        case 158
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated_dist3'];
        case 157
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated_dist4'];
        case 170
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated_dist2'];
        otherwise
            inflated_file = [params('rootdir') 'freesurfer/' subjid '/surf/' hemi '.inflated'];
    end
end
[vertices, ~] = freesurfer_read_surf(inflated_file);
nvertices = size(vertices,1); 