function flirt_sla2highres(exp,us,runtype,runnum,contr,model,varargin)

% function flirt_sla2highres(us,runtype,runnum,contr,model,varargin)
% 
% registers second-level contrasts to highres subject space
% 
% EDIT: should update second level gfeat to include runtype

%%
analysisdir = [params('rootdir') exp '/analysis/'];
datadir = [params('rootdir') exp '/data/'];
slasubdir = [analysisdir 'sla/usub' num2str(us) '/'];
fsl_version = read_fsl_version(exp,varargin{:});
fwhm = read_smooth(exp, varargin{:});

overwrite = false;
if optInputs(varargin,'overwrite');
    overwrite = true;
end

featdir = [slasubdir runtype '_' contr '_r' sprintf('%d',runnum) '_' model '_' num2str(100*fwhm, '%0.f') 'mm.gfeat/'];


thresh = [featdir 'cope1.feat/thresh_zstat1.nii.gz'];
zstat = [featdir 'cope1.feat/stats/zstat1.nii.gz'];
cope = [featdir 'cope1.feat/stats/cope1.nii.gz'];

threshHR = strrep(thresh, '.nii.gz', '_highres.nii.gz');
zstatHR = strrep(zstat, '.nii.gz', '_highres.nii.gz');
copeHR = strrep(cope, '.nii.gz', '_highres.nii.gz');

regmat = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(runnum(1)) '_' model  '_' num2str(100*fwhm, '%0.f') 'mm.feat/reg/standard2highres.mat'];
highres = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(runnum(1)) '_' model  '_' num2str(100*fwhm, '%0.f') 'mm.feat/reg/highres.nii.gz'];

if (overwrite || ~exist(threshHR,'file')) && exist(thresh,'file')
    unix_fsl(fsl_version, ['flirt -in ' thresh ' -ref ' highres ' -applyxfm -init ' regmat ' -out ' threshHR]);
end

if (overwrite || ~exist(zstatHR,'file')) && exist(zstat,'file')
    unix_fsl(fsl_version, ['flirt -in ' zstat ' -ref ' highres ' -applyxfm -init ' regmat ' -out ' zstatHR]);
end

if (overwrite || ~exist(copeHR,'file')) && exist(cope,'file')
    fprintf(['flirt -in ' cope ' -ref ' highres ' -applyxfm -init ' regmat ' -out ' copeHR '\n'])
    unix_fsl(fsl_version, ['flirt -in ' cope ' -ref ' highres ' -applyxfm -init ' regmat ' -out ' copeHR]);
end

threshHR_subsamp = strrep(thresh, '.nii.gz', '_highres_2mm.nii.gz');
zstatHR_subsamp = strrep(zstat, '.nii.gz', '_highres_2mm.nii.gz');
copeHR_subsamp = strrep(cope, '.nii.gz', '_highres_2mm.nii.gz');

if (overwrite || ~exist(threshHR_subsamp,'file')) && exist(threshHR,'file')
    unix_fsl(fsl_version, ['fslmaths ' threshHR ' -subsamp2 ' threshHR_subsamp]);
end

if (overwrite || ~exist(zstatHR_subsamp,'file')) && exist(zstatHR,'file')
    unix_fsl(fsl_version, ['fslmaths ' zstatHR ' -subsamp2 ' zstatHR_subsamp]);
end

if (overwrite || ~exist(copeHR_subsamp,'file')) && exist(copeHR,'file')
    unix_fsl(fsl_version, ['fslmaths ' copeHR ' -subsamp2 ' copeHR_subsamp]);
end

highrescopy = [slasubdir 'struct_r1_LAS_brain.nii.gz'];
highrescopy_subsamp = [slasubdir 'struct_r1_LAS_brain_2mm.nii.gz'];

if ~exist(highrescopy,'file')
    unix(['cp ' highres ' ' highrescopy]);
end 

if ~exist(highrescopy_subsamp,'file')
    unix_fsl(fsl_version, ['fslmaths ' highres ' -subsamp2 ' highrescopy_subsamp]);
end 