function fla_sigav_surf_v2(exp,us,runtype,r,hemi,input_fname,varargin)

%% Setup
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
subjid = [exp '_us' num2str(us)];
if optInputs(varargin, 'monkey')
    preprocdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
else
    preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
end

if optInputs(varargin, 'monkey')
    sigavdir = [params('rootdir') 'freesurfer/' subjid '/sigav/' runtype '_r' num2str(r) '/'];
else
    sigavdir = [params('rootdir') 'freesurfer/fsaverage/sigav/' subjid '/' runtype '_r' num2str(r) '/'];
end

if ~exist(sigavdir,'dir')
    mkdir(sigavdir);
end

%% Averaging
conds = read_conditions(exp,us,runtype,varargin{:});
t = read_timings(exp,us,runtype,r,varargin{:});
[blockdur, nulldur, TR, TA, stimdur, stim2scan, win] = read_scanparams(exp, us, runtype, 'run', r, varargin{:}); %#ok<ASGLU>

if optInputs(varargin, 'monkey')
    blocktp = find(win > 10 - 1e-3 & win < blockdur-TR+stim2scan + 1e-3);
else
    blocktp = find(win > 5 - 1e-3 & win < blockdur-TR+stim2scan + 1e-3);
end

sigav_file = [sigavdir hemi '.' input_fname '.mgz'];

if ~exist(sigav_file, 'file') || optInputs(varargin, 'overwrite')
    
    ts_file = [preprocdir hemi '.' input_fname '.mgz'];
    func = MRIread(ts_file);
    dims = size(func.vol);
    
    try
        nullav = zeros([dims(1:3), length(win)]);
        x = strcmp('NULL',t.conds);
        ons = repmat(win, sum(x), 1) + repmat(t.onsets(x), 1, length(win));
        ons_tr = round(ons/TR+1);
        for j = 1:length(win)
            nullav(:,:,:,j) = mean(func.vol(:,:,:,ons_tr(:,j)),4);
        end
    catch
        keyboard
    end
    
    blockav = zeros([dims(1:3), length(win), length(conds)]);
    for i = 1:length(conds)
        
        fprintf('Block average: %s\n',conds{i}); drawnow;
        x = strcmp(conds{i},t.conds);
        ons = repmat(win, sum(x), 1) + repmat(t.onsets(x), 1, length(win));
        ons_tr = round(ons/TR+1);
        
        for j = 1:length(win)
            blockav(:,:,:,j,i) = mean(func.vol(:,:,:,ons_tr(:,j)),4);
        end
    end
    
    cond_plateau = squeeze(mean(blockav(:,:,:,blocktp,:),4));
    null_plateau = repmat(squeeze(mean(nullav(:,:,:,blocktp),4))',[1 length(conds)]);
    
    sigav_plateau = func;
    sigav_plateau.vol = nan([dims(1:3),length(conds)]);
    sigav_plateau.vol(1,:,1,:) = (cond_plateau - null_plateau) ./ null_plateau;
    sigav_plateau.nframes = size(sigav_plateau.vol,4);
    sigav_plateau.fspec = sigav_file;
    MRIwrite(sigav_plateau, sigav_file);
    
end
