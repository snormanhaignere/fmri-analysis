function fla_sigav(exp,us,runtype,r,varargin)

% keyboard;
%%
analysisdir = [params('rootdir') exp '/analysis/'];
% featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' params('smooth') 'mm.feat/'];
fwhm = read_smooth(exp, varargin{:});

%%
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

output_file = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r'   num2str(r) '_' num2str(100*fwhm, '%.0f') 'mm_sigav_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.nii.gz'];

if ~exist(output_file, 'file') || optInputs(varargin, 'overwrite')
    
    
    %%
    if optInputs(varargin, 'monkey')
        func = readmr([analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/smooth' num2str(100*fwhm, '%.0f') 'mm_intnorm_detrend1.nii.gz'],'NOPROGRESSBAR');
    else
        func = readmr([analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/smooth' num2str(100*fwhm, '%.0f') 'mm_intnorm_hpfilt' num2str(params('hpcutoff')) '.nii.gz'],'NOPROGRESSBAR');
    end
    % func = readmr([analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/motcorr.nii.gz'],'NOPROGRESSBAR');
    dims = size(func.data);
    
    try
        nullav = zeros([dims(1:3), length(win)]);
        x = strcmp('NULL',t.conds);
        ons = repmat(win, sum(x), 1) + repmat(t.onsets(x), 1, length(win));
        ons_tr = round(ons/TR+1);
        for j = 1:length(win)
            nullav(:,:,:,j) = mean(func.data(:,:,:,ons_tr(:,j)),4);
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
            blockav(:,:,:,j,i) = mean(func.data(:,:,:,ons_tr(:,j)),4);
        end
    end
    
    cond_plateau = squeeze(mean(blockav(:,:,:,blocktp,:),4));
    null_plateau = repmat(squeeze(mean(nullav(:,:,:,blocktp),4)),[1 1 1 length(conds)]);
    
    sigav_plateau = func;
    sigav_plateau.data = (cond_plateau - null_plateau) ./ null_plateau;
    sigav_plateau.info.dimensions(4).size = size(sigav_plateau.data,4);
    bxhwrite(sigav_plateau,output_file);
    
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

