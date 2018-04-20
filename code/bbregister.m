function bbregister(exp,us,runtype,r,model,varargin)

% function bbregister(us,runtype,r,model,varargin)
%
% main function for running bbregister
% requires a surface reconstruction

%% directories, params
analysisdir = [params('rootdir') exp '/analysis/'];
freesurferdir = '~/freesurfer/';

fwhm = read_smooth(exp, varargin{:});
featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(fwhm*100,'%.0f') 'mm.feat'];
regdir = [featdir '/reg/bbreg/'];
subjid = [exp '_us' num2str(us)];

if ~exist(regdir,'dir');
    mkdir(regdir);
end

fsl_version = read_fsl_version(exp);
freesurfer_version = read_freesurfer_version(exp,varargin{:});

% transforms
exfunc2orig_freesurf = [regdir 'anat2exf.register.dat'];
exfunc2orig_fsl = [regdir 'exfunc2orig.mat']; %exfunc2orig
orig2highres_freesurf = [regdir 'orig2highres.dat'];
highres2orig_freesurf = [regdir 'highres2orig.dat'];
orig2highres_fsl = [regdir 'orig2highres.mat'];
highres2orig_fsl = [regdir 'highres2orig.mat'];
highres2exfunc_fsl = [regdir 'highres2example_func.mat'];
exfunc2highres_fsl = [regdir 'example_func2highres.mat'];

if exist(highres2exfunc_fsl, 'file') && exist(exfunc2highres_fsl, 'file') && ~optInputs(varargin, 'overwrite')
    return;
end

%% files

% images
orig = [freesurferdir subjid '/mri/orig.mgz']; % orig is the anatomical volume in the freesurfer space
if ~exist(orig, 'file')
    orig = [freesurferdir subjid '/mri/brain.mgz']; % orig is the anatomical volume in the freesurfer space
    fprintf('Using brain.mgz instead of orig.mgz\n'); drawnow;
end
highres = [featdir '/reg/highres.nii.gz'];

if ~optInputs(varargin,'partial_fov')
    %     exfunc_nosmooth = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/example_func.nii.gz'];
    %     if exist(exfunc_nosmooth,'file');
    %         exfunc = exfunc_nosmooth;
    %     else
    %         exfunc = [featdir '/reg/example_func.nii.gz'];
    %     end
    exfunc = [featdir '/reg/example_func.nii.gz'];
else
    exfunc = [featdir '/reg/example_func_partial_fov.nii.gz'];
    br = readmr([featdir '/reg/example_func.nii.gz'],'NOPROGRESSBAR');
    br.data(:,:,[1:17 20:31]) = 0;
    bxhwrite(br,exfunc);
    unix(['fslview ' exfunc ' &']);
end



% original flirt transforms
copyflirt(exp,us,runtype,r,model); % copy flirts for clarity and to prevent overwriting
if exist([featdir '/reg/flirt/'],'dir')
    exfunc2highres_init_fsl = [featdir '/reg/flirt/example_func2highres.mat'];
    exfunc2orig_init_fsl = [regdir 'example_func2orig_init.mat'];
    exfunc2orig_init_freesurf = [regdir 'example_func2orig_init.dat'];
else
    error('Should run swap first.')
end

%% generate initial registration matrix

% highres to orig, based on headers
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' highres ' --targ ' orig ' --reg ' highres2orig_freesurf ' --fslregout ' highres2orig_fsl ' --regheader --noedit ']);

% catenate header-based and flirt
unix_fsl(fsl_version, ['convert_xfm -omat ' exfunc2orig_init_fsl ' -concat ' highres2orig_fsl ' ' exfunc2highres_init_fsl]);

% create freesurfer-style matrix
unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_init_freesurf ' --fsl ' exfunc2orig_init_fsl ' --noedit ']);

% check initial registration
if optInputs(varargin,'check-init')
    tmp = ['us' num2str(us) '-r' num2str(r) '-init'];
    unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_init_freesurf ' --title ' tmp ' --surf &']);
end

%% registration

% main registration
if ~optInputs(varargin,'noreg')
    switch exp
        case 'pitch_overlap_v3'
            if (us == 104 && strcmp(runtype, 'overlap_v3_lowconds') && r == 6) || (us == 108 && strcmp(runtype, 'overlap_v3_highconds') && r == 5)
                unix_freesurfer_version( freesurfer_version, ['bbregister --s ' subjid ' --mov ' exfunc ' --t2 --init-reg ' exfunc2orig_init_freesurf ' --reg ' exfunc2orig_freesurf '  --brute1max 8']);
            else
                unix_freesurfer_version( freesurfer_version, ['bbregister --s ' subjid ' --mov ' exfunc ' --t2 --init-reg ' exfunc2orig_init_freesurf ' --reg ' exfunc2orig_freesurf '  --brute1max 6']);
            end
        otherwise
            unix_freesurfer_version( freesurfer_version, ['bbregister --s ' subjid ' --mov ' exfunc ' --t2 --init-reg ' exfunc2orig_init_freesurf ' --reg ' exfunc2orig_freesurf]);
    end
end

% check registration
if optInputs(varargin,'check-final')
    tmp = ['us' num2str(us) '-r' num2str(r) '-final'];
    unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' subjid ' --mov ' exfunc ' --reg ' exfunc2orig_freesurf ' --title ' tmp ' --surf &']);
end

% conversion process
% 1) output fsl file, default is no good for some reason...
unix_freesurfer_version( freesurfer_version, ['tkregister2 --mov ' exfunc ' --targ ' orig ' --reg ' exfunc2orig_freesurf ' --fslregout ' exfunc2orig_fsl ' --noedit ']);

% 2) orig to highres, based on headers
unix_freesurfer_version( freesurfer_version, ['tkregister2 --mov ' orig ' --targ ' highres ' --reg ' orig2highres_freesurf ' --fslregout ' orig2highres_fsl ' --regheader --noedit ']);

% 3) catenate header-based and bbregister-based transforms
unix_fsl(fsl_version, ['convert_xfm -omat ' exfunc2highres_fsl ' -concat ' orig2highres_fsl ' ' exfunc2orig_fsl]);

% 4) invert to get highres to exfunc transform
unix_fsl(fsl_version, ['convert_xfm -omat ' highres2exfunc_fsl ' -inverse ' exfunc2highres_fsl]);

% test it!
if optInputs(varargin,'test_highres2exfunc')
    highres_func = [regdir 'highres_func.nii.gz'];
    unix_fsl(fsl_version,['flirt -in ' highres ' -ref ' exfunc ' -applyxfm -init ' highres2exfunc_fsl ' -out ' highres_func]);
end

% test it!
if optInputs(varargin,'test_exfunc2highres')
    exfunc_highres = [featdir '/reg/freesurfer/example_func_highres.nii.gz'];
    unix_fsl(fsl_version, ['flirt -in ' exfunc ' -ref ' highres ' -applyxfm -init ' exfunc2highres_fsl ' -out ' exfunc_highres]);
end

if ~optInputs(varargin,'noswap')
    bbregister_flirtfiles(exp,us,runtype,r,model,varargin{:});
    regswap(exp,us,runtype,r,model,'bbreg',varargin{:});
    featregapply(exp,us,runtype,r,model,varargin{:});
end

