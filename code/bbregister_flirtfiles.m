function bbregister_flirtfiles(exp,us,runtype,r,model,varargin)

% bbregister_flirtfiles(exp,us,runtype,r,model,varargin)
% 
% generates extra registration files based on bbreg transformation
% from example_func to highres 

scriptsdir = [pwd '/'];
analysisdir = [params('rootdir') exp '/analysis/'];
fwhm = read_smooth(exp, varargin{:});

fsl_version = read_fsl_version(exp);
featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/'];
regdir = [featdir 'reg/bbreg/'];

%% files

% files based on flirt's reg
exfunc = [featdir 'reg/example_func.nii.gz'];
highres = [featdir 'reg/highres.nii.gz'];
standard = [featdir 'reg/standard.nii.gz'];
highres2standard = [featdir 'reg/highres2standard.mat'];

% files based on bbreg
exfunc2highres = [regdir 'example_func2highres.mat'];
exfunc2highres_brain = [regdir 'example_func2highres.nii.gz'];
exfunc2highres_pic1 = [regdir 'example_func2highres1.png'];
exfunc2highres_pic2 = [regdir 'example_func2highres2.png'];
exfunc2highres_pic3 = [regdir 'example_func2highres.png'];

exfunc2standard = [regdir 'example_func2standard.mat'];
standard2exfunc = [regdir 'standard2example_func.mat'];
exfunc2standard_brain = [regdir 'example_func2standard.nii.gz'];
exfunc2standard_pic1 = [regdir 'example_func2standard1.png'];
exfunc2standard_pic2 = [regdir 'example_func2standard2.png'];
exfunc2standard_pic3 = [regdir 'example_func2standard.png'];

%% make the files

fprintf('Registering example_func2highres files...');

unix_fsl(fsl_version,['flirt -in ' exfunc ' -ref ' highres ' -applyxfm -init ' exfunc2highres ' -out ' exfunc2highres_brain]);

fprintf('Registering example_func2standard files...');

unix_fsl(fsl_version,['convert_xfm -omat ' exfunc2standard ' -concat ' highres2standard ' ' exfunc2highres]);
unix_fsl(fsl_version,['convert_xfm -omat ' standard2exfunc ' -inverse ' exfunc2standard]);
unix_fsl(fsl_version,['flirt -in ' exfunc ' -ref ' standard ' -applyxfm -init ' exfunc2standard ' -out ' exfunc2standard_brain]);

%% slicefiles

fprintf('Making example_func2highres pics...');

addpath(scriptsdir);
cd(regdir);

unix_fsl(fsl_version,['slicer ' exfunc2highres_brain ' ' highres ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ']); 
unix_fsl(fsl_version,['pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png ' exfunc2highres_pic1]);
unix_fsl(fsl_version,['slicer ' highres ' ' exfunc2highres_brain ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png']);
unix_fsl(fsl_version,['pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png ' exfunc2highres_pic2]);
unix_fsl(fsl_version,['pngappend ' exfunc2highres_pic1 ' - ' exfunc2highres_pic2 ' ' exfunc2highres_pic3 ' ; rm -f sl?.png']);

fprintf('Making example_func2standard pics...');

unix_fsl(fsl_version,['slicer ' exfunc2standard_brain ' ' standard ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ']); 
unix_fsl(fsl_version,['pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png ' exfunc2standard_pic1]);
unix_fsl(fsl_version,['slicer ' standard ' ' exfunc2standard_brain ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png']);
unix_fsl(fsl_version,['pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png ' exfunc2standard_pic2]);
unix_fsl(fsl_version,['pngappend ' exfunc2standard_pic1 ' - ' exfunc2standard_pic2 ' ' exfunc2standard_pic3 ' ; rm -f sl?.png']);

cd(scriptsdir);
%% scraps



% % example func to highres, copy over from freesurfer parent directory
% tmp = [featdir 'reg/freesurfer/example_func2highres.mat'];
% unix(['cp -f ' tmp ' ' exfunc2highres]);
% 
% tmp = [featdir 'reg/freesurfer/highres2example_func.mat'];
% unix(['cp -f ' tmp ' ' highres2exfunc]);
% highres2exfunc = [regdir 'highres2example_func.mat'];


% regtype = 'bb_spm'; % 4 possible types, bb_spm, bb_flirt, nobb_spm, nobb_flirt
% 
% if optInputs(varargin,'regtype');
%     regtype = varargin{optInputs(varargin,'regtype')+1};
% end

