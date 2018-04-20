function rfx_surf(expsub,con,fwhm,hemi,runtype,model,varargin)

%% Paramaters

projfrac = [0 1];

% whether to use weighted least squares (default)
% or ordinary least squares
glmtype = 'wls';
if optInputs(varargin,'ols')
    glmtype = 'ols';
end

fwhm_vol = read_smooth(expsub{1,1}, varargin{:});

% string to identify this particular analysis
idstring = [hemi '_' con '_usubs' sprintf('%d',sort(cat(1,expsub{:,2})))  '_fwhm' num2str(fwhm) '_' glmtype '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm' ];

% output glm directory
exp_unique = unique(expsub(:,1));
glmdir = ['~/freesurfer/fsaverage/group_rfx/' sprintf('%s_', exp_unique{:}) 'us' sprintf('%d',expsub{:,2}) '/' idstring '/'];
if ~exist(glmdir,'dir');
    mkdir(glmdir);
end

% p value to use for correction
pval = 0.05;
if optInputs(varargin,'pval')
    pval = varargin{optInputs(varargin,'pval')+1};
end

correction = 'cache';
if optInputs(varargin,'correction')
    correction = varargin{optInputs(varargin,'correction')+1};
end

con_sign = 'pos';
if optInputs(varargin,'con_sign')
    con_sign = varargin{optInputs(varargin,'con_sign')+1};
end

voxthresh = 3;
if optInputs(varargin,'voxthresh')
    voxthresh = varargin{optInputs(varargin,'voxthresh')+1};
end

nsim = 5000;
if optInputs(varargin,'nsim')
    nsim = varargin{optInputs(varargin,'nsim')+1};
end

% min and max used for plotting
if any(strcmp(con_sign,{'pos','neg'}))
    fminmax = [voxthresh,6]-log10(2);
elseif strcmp(con_sign,{'abs'})
    fminmax = [voxthresh, 6];
end

freesurfer_version = read_freesurfer_version(expsub{1,1},varargin{:});

%% Preprocesing

copestr = [];
varstr = [];
for i = 1:length(expsub)
    
    % freesurfer subject id
    subjid = [expsub{i,1} '_us' num2str(expsub{i,2})];
    
    % runs for that subject
    allruns = read_runs(expsub{i,1}, expsub{i,2}, runtype, varargin{:});
    
    % second level copes and varcopes
    x = [con '_fwhm0_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_' model  '_' num2str(fwhm_vol*100,'%.0f') 'mm'];
    copesurf = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '_' x '_r' sprintf('%d',allruns) '/beta.mgh'];
    varsurf = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '_' x '_r' sprintf('%d',allruns) '/rvar.mgh'];
    
    if ~exist(copesurf,'file') || ~exist(varsurf,'file')
        keyboard;
        error('Cannot find second level analyses for %s, usub %d', expsub{i,1},expsub{i,2});
    end
    
    copestr = [copestr ' --s fsaverage --is ' copesurf]; %#ok<AGROW>
    varstr = [varstr ' --s fsaverage --is ' varsurf]; %#ok<AGROW>
    
end

allcopes = [glmdir  'allcopes.mgh'];
allvarcopes = [glmdir  'allvarcopes.mgh'];

if fwhm > 0;
    fwhm_str = ['--fwhm ' num2str(fwhm)];
else
    fwhm_str = '';
end

if ~exist(allcopes,'file') || ~exist(allvarcopes,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version,['mris_preproc --target fsaverage --hemi ' hemi ' --out ' allcopes ' ' copestr ' ' fwhm_str]);
    unix_freesurfer_version(freesurfer_version,['mris_preproc --target fsaverage --hemi ' hemi ' --out ' allvarcopes ' ' varstr ' ' fwhm_str]);
end

%% GLM

sigmap = [glmdir 'osgm/sig.mgh'];
if ~exist(sigmap,'file')
    if strcmp(glmtype, 'wls')
        unix_freesurfer_version(freesurfer_version,['mri_glmfit --y ' allcopes ' --wls ' allvarcopes ' --osgm  --surf fsaverage ' hemi ' --glmdir ' glmdir ' --seed 0 ' ]);
    elseif strcmp(glmtype,'ols')
        unix_freesurfer_version(freesurfer_version,['mri_glmfit --y ' allcopes ' --osgm  --surf fsaverage ' hemi ' --glmdir ' glmdir ' --seed 0 ' ]);
    else
        error('No valid glm type');
    end
end

titlestring = strrep(con,'_','-');
anatarg = [' -annotation ~/freesurfer/fsaverage/label/' hemi '.stp-pm2al-5lab-recolor.annot -labels-under '];
if optInputs(varargin,'sigmap');
    if optInputs(varargin, 'tksurfer')
        unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -overlay ' sigmap ' -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' -title ' titlestring ' ' anatarg ' &']);
    else
        freeview3('fsaverage',hemi,'overlay',sigmap,'overlay_threshold',[fminmax(1), mean(fminmax), fminmax(2)],varargin{:});
    end
end

if optInputs(varargin,'nocorrection')
    return;
end

%% Cluster Correction

% base to identify the analysis
if strcmp(correction, 'cache')
    csdbase = ['cache.th' strrep(sprintf('%2.1f',voxthresh), '.', '') '.' con_sign];
else
    csdbase = [correction '_' con_sign '_thr' strrep(num2str(voxthresh),'.','') '_nsim' num2str(nsim)];
end

% cluster correction
sigmap_mask = [glmdir 'osgm/' csdbase '.sig.masked.mgh'];
cluster_annot = [glmdir 'osgm/' csdbase '.sig.ocn.annot'];

if strcmp(correction, 'cache')
    unix_freesurfer_version(freesurfer_version,['mri_glmfit-sim --glmdir ' glmdir ' --cache ' num2str(voxthresh) ' ' con_sign ' --2spaces --cwpvalthresh ' num2str(pval) ' --overwrite ']);
elseif ~exist(sigmap_mask,'file') || ~exist(cluster_annot,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version,['mri_glmfit-sim --glmdir ' glmdir ' --sim ' correction ' ' num2str(nsim) ' ' num2str(voxthresh) ' ' csdbase ' --sim-sign ' con_sign ' --2spaces --overwrite ']);
else
    % just recompute p values
    unix_freesurfer_version(freesurfer_version,['mri_glmfit-sim --glmdir ' glmdir ' --no-sim ' csdbase ' --2spaces --cwpvalthresh ' num2str(pval)]);
end

% plot the masked significance map
if optInputs(varargin,'sigmap_mask');
    if optInputs(varargin, 'tksurfer')
        unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -overlay ' sigmap_mask ' -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' -title ' titlestring ' ' anatarg ' &']);
    else
        freeview3('fsaverage',hemi,'overlay',sigmap_mask,'overlay_threshold',[fminmax(1), mean(fminmax), fminmax(2)],varargin{:});
    end
end

% plot the annotation map
if optInputs(varargin,'annotation');
    unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -gray -annotation ' cluster_annot ' &']);
end

%% Binary label

addpath('/software/Freesurfer/emf-5.1.0/matlab');

if optInputs(varargin, 'label')
    % binarize
    sigmap_mask_bin = strrep(sigmap_mask,'.mgh','.bin.mgh');
    if ~exist(sigmap_mask_bin, 'file') || optInputs(varargin,'overwrite')
        [srf, M, mr_parms] = load_mgh(sigmap_mask);
        srf(srf<fminmax(1)) = 0;
        srf(srf>=fminmax(1)) = 1;
        save_mgh(srf,sigmap_mask_bin,M,mr_parms);
    end
    
    % optionally smooth binarized map and rebinarize
    if optInputs(varargin,'smoothbin')
        
        smoothkern = varargin{optInputs(varargin,'smoothbin') + 1};
        sigmap_mask_bin_smooth = strrep(sigmap_mask_bin,'.mgh',['.smooth' num2str(smoothkern) '.mgh']);
        if ~exist(sigmap_mask_bin_smooth,'file') || optInputs(varargin,'overwrite')
            unix_freesurfer_version(freesurfer_version,['mri_surf2surf --s fsaverage --sval ' sigmap_mask_bin ' --tval ' sigmap_mask_bin_smooth  ' --hemi ' hemi ' --fwhm-src ' num2str(smoothkern)]);
        end
        
        % re-binarize
        sigmap_mask_bin_smooth_bin = strrep(sigmap_mask_bin_smooth,'.mgh','.bin.mgh');
        [srf, M, mr_parms] = load_mgh(sigmap_mask_bin_smooth);
        srf(srf<0.25) = 0;
        srf(srf>=0.25) = 1;
        save_mgh(srf,sigmap_mask_bin_smooth_bin,M,mr_parms);
    else
        sigmap_mask_bin_smooth_bin = sigmap_mask_bin;
    end
    
    % produce individual labels
    label = strrep(sigmap_mask_bin_smooth_bin, '.mgh','.label');
    if ~exist(label,'file') || optInputs(varargin,'overwrite')
        labeldir = [strrep(label,'.','_') '/'];
        if ~exist(labeldir,'dir');
            mkdir(labeldir);
        end
        
        if optInputs(varargin,'overwrite')
            unix(['rm -f ' labeldir '*']);
        end
        
        labelbase = [labeldir 'sigclusters'];
        unix_freesurfer_version(freesurfer_version,['mri_surfcluster --in ' sigmap_mask_bin_smooth_bin ' --olab ' labelbase ' --subject fsaverage --hemi ' hemi ' --thmin 0.5 ' ]);
        
        % combine labels if necessary
        labels = mydir(labeldir);
        if length(labels) > 1
            unix_freesurfer_version(freesurfer_version,['mri_mergelabels ' sprintf([' -i ' labeldir '%s'],labels{:}) ' -o ' label]);
        else
            unix(['cp -f ' labeldir labels{1} ' ' label]);
        end
        
    end
    
    if optInputs(varargin, 'tksurfer')
        unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -overlay ' sigmap_mask ' -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' -label ' label ' -label-outline &']);
    else
        try
            freeview(sigmap_mask,hemi,[fminmax(1), mean(fminmax), fminmax(2)],'label',label);
        catch
            keyboard;
        end
    end
end


%% Scraps

% unix_freesurfer_version(freesurfer_version,['mri_glmfit-sim --glmdir ' groupdir ' --cache ' num2str(voxthresh) ' ' con_sign ' --cwpvalthresh ' num2str(pval/2)]);
% voxthresh_str = strrep(sprintf('%2.1f',voxthresh), '.', '');
% fprintf(['tksurfer fsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.' con_sign '.sig.masked.mgh -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' &\n']);
%
% % raw significance map
% unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/sig.mgh -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' &']);
%
% % thresholded significance map
% unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.' con_sign '.sig.masked.mgh -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' &']);
%
% % significant clusters
% unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -annotation ' groupdir 'osgm/cache.th' voxthresh_str '.' con_sign '.sig.ocn.annot &']);

% fprintf(['tksurfer fsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.' con_sign '.sig.ocn.mgh -fminmax 0.5 5.5 &\n']);
% unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.con_sign '.sig.ocn.mgh -fminmax 0.5 5.5 &']);
% unix_freesurfer_version(freesurfer_version,['tksurfer fsaverage ' hemi ' inflated -annot ' groupdir 'osgm/cache.th' voxthresh_str '.' con_sign '.sig.ocn.annot &']);

%
%     cd(['~/' expsub{i,1} '/scripts/']);
%     allruns = read_runs(expsub{i,2}, runtype);
%     freesurfreg = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(fwhm_vol*100,'%.0f') 'mm.feat/reg/anat2std.register.dat'];
%
%     if ~exist(freesurfreg,'file') || optInputs(varargin,'overwrite')
%         highres = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(fwhm_vol*100,'%.0f') 'mm.feat/reg/highres.nii.gz'];
%         standard = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(fwhm_vol*100,'%.0f') 'mm.feat/reg/standard.nii.gz'];
%         fslreg = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(fwhm_vol*100,'%.0f') 'mm.feat/reg/highres2standard.mat'];
%         unix(['tkregister2 --mov ' highres ' --targ ' standard ' --fsl ' fslreg ' --reg ' freesurfreg ' --noedit ']);
%     end