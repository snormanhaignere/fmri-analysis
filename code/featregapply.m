function featregapply(exp,us,runtype,r,model,varargin)

% function featregapply(s,runtype,r,model)
% 
% runs fsl's featregapply

fwhm = read_smooth(exp, varargin{:});
fsl_version = read_fsl_version(exp, varargin{:});
analysisdir = [params('rootdir') exp '/analysis/'];
featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/'];
unix_fsl(fsl_version,['featregapply ' featdir]);
