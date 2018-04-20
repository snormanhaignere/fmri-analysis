function fla_sigav_surf(exp,us,runtype,r,hemi,varargin)

%% Setup
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
freesurfer_version = read_freesurfer_version(exp,varargin{:});
subjid = [exp '_us' num2str(us)];
if optInputs(varargin, 'monkey')
    preprocdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
else
    preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
end
[~,fwhm] = read_smooth(exp, varargin{:});
if optInputs(varargin, 'fwhm');
    fwhm = varargin{optInputs(varargin, 'fwhm')+1};
end

cutoff = inf;
if optInputs(varargin, 'cutoff')
    cutoff = varargin{optInputs(varargin, 'cutoff')+1};
end

if optInputs(varargin, 'monkey')
    sigavdir = [params('rootdir') 'freesurfer/' subjid '/fla/' runtype '_r' num2str(r) '/'];
else
    sigavdir = [params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/' runtype '_r' num2str(r) '/'];
end

if ~exist(sigavdir,'dir')
    mkdir(sigavdir);
end

demean_flag = '';
if optInputs(varargin,'demean')
    demean_flag = '_demean';
end
if optInputs(varargin, 'nulldemean')
    demean_flag = '_nulldemean';
end

% analysisdir = [params('rootdir') exp '/analysis/'];
% featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' params('smooth') 'mm.feat/'];
% fwhm = read_smooth(exp, varargin{:});

%% Averaging
conds = read_conditions(exp,us,runtype,varargin{:});
t = read_timings(exp,us,runtype,r,varargin{:});
[blockdur, nulldur, TR, TA, stimdur, stim2scan, win] = read_scanparams(exp, us, runtype, 'run', r, varargin{:}); %#ok<ASGLU>

if optInputs(varargin, 'monkey')
    blocktp = find(win > 10 - 1e-3 & win < blockdur-TR+stim2scan + 1e-3);
else
    blocktp = find(win > 5 - 1e-3 & win < blockdur-TR+stim2scan + 1e-3);
end

% win = win(win > 5 - 1e-3 & win < blockdur-TR+stim2scan + 1e-3);
% win = win;

% sigav_file = [sigavdir hemi '.sigav' demean_flag '_' num2str(100*fwhm, '%.0f') 'mm_tp' num2str(blocktp(1)) '-' num2str(blocktp(end))  '.mgz'];
sigav_file = [sigavdir hemi '.sigav' demean_flag '_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-' num2str(cutoff) 'sec-order2_tp' num2str(blocktp(1)) '-' num2str(blocktp(end))  '.mgz'];

if ~exist(sigav_file, 'file') || optInputs(varargin, 'overwrite')
    
%     %%
    %     if optInputs(varargin, 'monkey')
    %         func = readmr([preprocdir 'smooth.nii.gz'],'NOPROGRESSBAR');
    %     else
    %         func = readmr([preprocdir 'smooth' num2str(100*fwhm, '%.0f') 'mm.nii.gz'],'NOPROGRESSBAR');
    %     end
    
    %     ts_file = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/preprocess/naturalsound_us1/main_v3_r1/rh.brain_thresh_detrend1.mgz';
    %     ts_file = '/mindhive/nklab/u/svnh/freesurfer/naturalsound_us1/preprocess/usub1/main_v3_r1/rh.brain_thresh_detrend1.mgz';
    %     ts_file = '/mindhive/nklab/u/svnh/freesurfer/fsaverage/preprocess/naturalsound_us1/main_v3_combined_r1/rh.brain_thresh_detrend1_demean.mgz';
    
    %     ts_file = [preprocdir hemi '.brain_thresh_detrend1_demean.mgz'];
    if fwhm == 0 && cutoff == inf
        ts_file = [preprocdir hemi '.brain_thresh_detrend1' demean_flag '.mgz'];
    else
        ts_file = [preprocdir hemi '.brain_thresh_detrend1' demean_flag  '_smooth' num2str(100*fwhm, '%.0f') 'mm_hpfilt-' num2str(cutoff) 'sec-order2.mgz'];
    end
    %         ts_file = [preprocdir hemi '.brain_thresh_detrend1' demean_flag  '_smooth' num2str(100*fwhm, '%.0f') 'mm.mgz'];
    func = MRIread(ts_file);
    % func = readmr([analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/motcorr.nii.gz'],'NOPROGRESSBAR');
    dims = size(func.vol);
    
    %   try
    nullav = zeros([dims(1:3), length(win)]);
    x = strcmp('NULL',t.conds);
    ons = repmat(win, sum(x), 1) + repmat(t.onsets(x), 1, length(win));
    ons_tr = round(ons/TR+1);
    for j = 1:length(win)
        nullav(:,:,:,j) = mean(func.vol(:,:,:,ons_tr(:,j)),4);
    end
    %   catch
    %     keyboard
    %   end
    
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
    
    %   sigav_plateau.data = (cond_plateau - null_plateau) ./ null_plateau;
    %   sigav_plateau.info.dimensions(4).size = size(sigav_plateau.data,4);
    %   bxhwrite(sigav_plateau,sigav_file);
    sigav_plateau = func;
    sigav_plateau.vol = nan([dims(1:3),length(conds)]);
    sigav_plateau.vol(1,:,1,:) = (cond_plateau - null_plateau) ./ null_plateau;
    sigav_plateau.nframes = size(sigav_plateau.vol,4);
    sigav_plateau.fspec = sigav_file;
    MRIwrite(sigav_plateau, sigav_file);
    
end

%
% exfunc =readmr([analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/example_func.nii.gz'],'NOPROGRESSBAR');
% showsrs2(exfunc, sigav_timecourse);


% %% ttest
%
% pbr = func;
% pbr.info.dimensions(4).size = 1;
% [h p] = ttest(squeeze(blockav(:,:,:,:,1)-blockav(:,:,:,:,2)),0,0.05,'both',4);
% pbr.data = -log10(p);
% bxhwrite(pbr,[featdir 'stats/face_obj_sigav.nii.gz']);

