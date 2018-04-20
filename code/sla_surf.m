function sla_surf(exp,us,con,fwhm_surf,correction,con_sign,voxthresh,nsim,hemi,runtype,model,varargin)


%% Paramaters

analysisdir = [params('rootdir') exp '/analysis/'];
subjid = [exp '_us' num2str(us)];

% analysis can either be done on either the template
% or the individual subject brain
if optInputs(varargin,'fsaverage')
    targid = 'fsaverage';
elseif optInputs(varargin,'myfsaverage')
    targid = 'myfsaverage';
else
    targid = subjid;
end

fwhm = read_smooth(exp, varargin{:});

runnum = read_runs(exp,us,runtype,varargin{:});
if optInputs(varargin, 'runs')
    runnum = varargin{optInputs(varargin, 'runs')+1};
end

projfrac = [0 1];
interpmethod = 'trilinear';

% string used to identify this analysis
idstring = [con '_fwhm' num2str(fwhm_surf) '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_' interpmethod '_' runtype '_' model  '_' num2str(100*fwhm, '%.0f') 'mm'];

% output directory
if optInputs(varargin,'fsaverage')
    glmdir = ['~/freesurfer/fsaverage/sla/' subjid '/' hemi '_' idstring '_r' sprintf('%d',runnum)];
elseif optInputs(varargin,'myfsaverage')
    glmdir = ['~/freesurfer/myfsaverage/sla/' subjid '/' hemi '_' idstring '_r' sprintf('%d',runnum)];
else
    glmdir = ['~/freesurfer/' subjid '/sla/' hemi '_' idstring '_r' sprintf('%d',runnum)];
end

if length(runnum) > 100
    glmdir = strrep(glmdir, ['_r' sprintf('%d',runnum)], ['_r' DataHash(runnum)]);
end

try
    % make output directory
    if ~exist(glmdir,'dir');
        mkdir(glmdir);
    end
catch
    keyboard;
end

% clean out directory
if optInputs(varargin,'overwrite')
    unix(['rm -rf ' glmdir '/*']);
end

% directory used to store first level surface files
fladir = ['~/freesurfer/' subjid '/fla/'];
if ~exist(fladir,'dir');
    mkdir(fladir);
end

vox_saturation_threshold = 6;
if optInputs(varargin, 'vox_saturation_threshold')
    vox_saturation_threshold = varargin{optInputs(varargin, 'vox_saturation_threshold')+1};
end

% min and max used for plotting
if any(strcmp(con_sign,{'pos','neg'}))
    fminmax = [voxthresh, vox_saturation_threshold]-log10(2);
elseif strcmp(con_sign,{'abs'})
    fminmax = [voxthresh, vox_saturation_threshold];
end

pval = 0.05;
if optInputs(varargin,'pval')
    pval = varargin{optInputs(varargin,'pval')+1};
end

figure_directory = [params('rootdir') '' exp '/figures/sla_surf/'];
if ~exist(figure_directory,'dir')
    mkdir(figure_directory);
end

freesurfer_version = read_freesurfer_version(exp,varargin{:});

%% Surface Projection

copestr = [];
varstr = [];
dof = nan(length(runnum),1);
for i = 1:length(runnum)
    
    % freesurfer subject id
    fprintf('Surface projetion: run %d\n',runnum(i)); drawnow;
    
    % registration from the functional LAS space into the anatomical RAS space
    regfile = regfile_funcLAS2anatRAS(exp,us,runtype,runnum(i),model,varargin{:});
    
    % cope and variance volumes
    featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(runnum(i)) '_' model  '_' num2str(100*fwhm, '%.0f') 'mm.feat/'];
    fid = fopen([featdir 'contrastnames.txt'],'r');
    tmp = textscan(fid,'%s\n'); fclose(fid); contrastnames = tmp{1};
    ind = strmatch(con,contrastnames,'exact');
    if length(ind) > 1
        warning('Multiple contrast matches');
        ind = ind(1);
    end
    ind = ind(1);
    copevol = [featdir 'stats/cope' num2str(ind) '.nii.gz'];
    varvol = [featdir 'stats/varcope' num2str(ind) '.nii.gz'];
    
    % degrees of freedom
    try
        fid = fopen([featdir 'stats/dof'],'r');
        tmp = textscan(fid,'%d'); fclose(fid);
        dof(i) = tmp{1};
    catch
        fprintf('Error in sla_surf: can''t find file %s\n', [featdir 'stats/dof']); drawnow;
        keyboard;
    end
    
    
    % sample cope and varcope files to the surface
    copesurf = [fladir hemi '.cope_' idstring  '_r' num2str(runnum(i)) '.mgz'];
    varsurf = [fladir hemi '.varcope_' idstring  '_r' num2str(runnum(i)) '.mgz'];
    if ~exist(copesurf,'file') || ~exist(varsurf,'file') || optInputs(varargin,'overwrite')
        refvol = ['~/freesurfer/' subjid '/mri/orig.mgz'];
        unix_freesurfer_version(freesurfer_version, ['mri_vol2surf --ref ' refvol ' --mov ' copevol ' --reg ' regfile ' --o ' copesurf ' --hemi ' hemi ' --surf white --interp ' interpmethod ' --projfrac-avg ' num2str(projfrac(1)) ' ' num2str(projfrac(2)) ' 0.05']);
        unix_freesurfer_version(freesurfer_version, ['mri_vol2surf --ref ' refvol ' --mov ' varvol ' --reg ' regfile ' --o ' varsurf ' --hemi ' hemi ' --surf white --interp ' interpmethod ' --projfrac-avg ' num2str(projfrac(1)) ' ' num2str(projfrac(2)) ' 0.05']);
    end
    
    % store subject ids
    if optInputs(varargin, 'no_sphere_reg')
        copestr = [copestr ' ' copesurf]; %#ok<AGROW>
        varstr = [varstr ' ' varsurf]; %#ok<AGROW>
    else
        copestr = [copestr ' --s ' subjid ' --is ' copesurf]; %#ok<AGROW>
        varstr = [varstr ' --s ' subjid ' --is ' varsurf]; %#ok<AGROW>
    end
end

%% Preprocesing

% combine and optionally smooth copes/varcopes
allcopes = [glmdir '/allcopes.mgz'];
allvarcopes = [glmdir '/allvarcopes.mgz'];

if fwhm_surf > 0;
    fwhm_str = ['--fwhm ' num2str(fwhm_surf)];
else
    fwhm_str = '';
end

if ~exist(allcopes,'file') || ~exist(allvarcopes,'file') || optInputs(varargin,'overwrite')
    if optInputs(varargin, 'no_sphere_reg')
        if fwhm_surf > 0;
            fprintf('Need to implement smoothing\n'); drawnow;
            keyboard;
        end
        unix_freesurfer_version(freesurfer_version, ['mri_concat ' copestr ' --o ' allcopes]);
        unix_freesurfer_version(freesurfer_version, ['mri_concat ' varstr ' --o ' allvarcopes]);
    else
        unix_freesurfer_version(freesurfer_version, ['mris_preproc --target ' targid ' --hemi ' hemi ' --out ' allcopes ' ' copestr ' ' fwhm_str]);
        unix_freesurfer_version(freesurfer_version, ['mris_preproc --target ' targid ' --hemi ' hemi ' --out ' allvarcopes ' ' varstr ' ' fwhm_str]);
    end
end

%% GLM

sigmap = [glmdir '/osgm/sig.mgh'];
if ~exist(sigmap,'file') || optInputs(varargin, 'overwrite')
    unix_freesurfer_version(freesurfer_version, ['mri_glmfit --y ' allcopes ' --yffxvar ' allvarcopes ' --ffxdof ' num2str(sum(dof)) ' --osgm  --surf ' targid ' ' hemi ' --seed 0 --glmdir ' glmdir]);
end

% anatarg = [' -annotation ~/freesurfer/fsaverage/label/' hemi '.stp-pm2al-5lab-recolor.annot -labels-under '];
% titlestring = ['us' num2str(us) '-' strrep(con,'_','-')];

if optInputs(varargin,'sigmap_screenshot');
    
    sigmap_file = [figure_directory subjid '_' hemi '_' idstring '_r' sprintf('%d',runnum)];
    if length(runnum) > 100
        sigmap_file = strrep(sigmap_file, ['_r' sprintf('%d',runnum)], ['_r' DataHash(runnum)]);
    end
    freeview3(targid,hemi,'overlay',sigmap,'overlay_threshold',[fminmax(1), mean(fminmax), fminmax(2)],'screenshot',sigmap_file);
end
    
if optInputs(varargin,'sigmap');
    %     if optInputs(varargin, 'tksurfer')
    %         unix_freesurfer_version(freesurfer_version, ['tksurfer ' targid ' ' hemi ' inflated -overlay ' sigmap ' -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' -title ' titlestring ' ' anatarg ' &']);
    %     else
    sigmap_file = [figure_directory subjid '_' hemi '_' idstring '_r' sprintf('%d',runnum)];
    if length(runnum) > 100
        sigmap_file = strrep(sigmap_file, ['_r' sprintf('%d',runnum)], ['_r' DataHash(runnum)]);
    end
    %     freeview3(targid,hemi,'overlay',sigmap,'overlay_threshold',[fminmax(1), mean(fminmax), fminmax(2)],'screenshot',sigmap_file);
    freeview3(targid,hemi,'overlay',sigmap,'overlay_threshold',[fminmax(1), mean(fminmax), fminmax(2)]);
    %     end
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

sigmap_mask = [glmdir '/osgm/' csdbase '.sig.masked.mgh'];
cluster_annot = [glmdir '/osgm/' csdbase '.sig.ocn.annot'];

if strcmp(correction, 'cache')
    unix_freesurfer_version(freesurfer_version, ['mri_glmfit-sim --glmdir ' glmdir ' --cache ' num2str(voxthresh) ' ' con_sign ' --2spaces --overwrite --cwpvalthresh ' num2str(pval) ]);
elseif ~exist(sigmap_mask,'file') || ~exist(cluster_annot,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ['mri_glmfit-sim --glmdir ' glmdir ' --sim ' correction ' ' num2str(nsim) ' ' num2str(voxthresh) ' ' csdbase ' --sim-sign ' con_sign ' --2spaces --overwrite -cwpvalthresh ' num2str(pval)]);
else
    unix_freesurfer_version(freesurfer_version, ['mri_glmfit-sim --glmdir ' glmdir ' --no-sim ' csdbase ' --2spaces --cwpvalthresh ' num2str(pval)]);
end


if optInputs(varargin,'sigmap_mask');
    if optInputs(varargin, 'tksurfer')
        unix_freesurfer_version(freesurfer_version, ['tksurfer ' targid ' ' hemi ' inflated -overlay ' sigmap_mask ' -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' -title ' titlestring ' ' anatarg ' &']);
    else
        freeview(sigmap_mask,hemi,[fminmax(1), mean(fminmax), fminmax(2)],varargin{:});
    end
end

if optInputs(varargin,'annotation');
    unix_freesurfer_version(freesurfer_version, ['tksurfer ' targid ' ' hemi ' inflated -gray -annotation ' cluster_annot ' -title ' titlestring ' &']);
end

lb = read_label(['~/freesurfer/fsaverage/label/' hemi '.stp.label']);
sigmap_labelmask = [glmdir '/osgm/sig_stp_v3.mgh'];
sigmapB = MRIread(sigmap);
sigmapB.fspec = sigmap_labelmask;
xi = setdiff(1:length(sigmapB.vol), lb.vnums+1);
sigmapB.vol(xi) = 0;
sigmapB.vol(sigmapB.vol<fminmax(1)) = 0;
MRIwrite(sigmapB, sigmap_labelmask);


%%

if optInputs(varargin, 'label')
    
    
    addpath('/software/Freesurfer/emf-5.1.0/matlab');
    
    % binarize
    sigmap_mask_bin = strrep(sigmap_labelmask,'.mgh','.bin.mgh');
    if ~exist(sigmap_mask_bin, 'file') || optInputs(varargin,'overwrite')
        [srf, M, mr_parms] = load_mgh(sigmap_labelmask);
        srf(srf<fminmax(1)) = 0;
        srf(srf>=fminmax(1)) = 1;
        save_mgh(srf,sigmap_mask_bin,M,mr_parms);
    end
    
    % optionally smooth binarized map and rebinarize
    if optInputs(varargin,'smoothbin')
        
        smoothkern = varargin{optInputs(varargin,'smoothbin') + 1};
        sigmap_mask_bin_smooth = strrep(sigmap_mask_bin,'.mgh',['.smooth' num2str(smoothkern) '.mgh']);
        if ~exist(sigmap_mask_bin_smooth,'file') || optInputs(varargin,'overwrite')
            unix_freesurfer_version(freesurfer_version, ['mri_surf2surf --s fsaverage --sval ' sigmap_mask_bin ' --tval ' sigmap_mask_bin_smooth ' --hemi ' hemi ' --fwhm-src ' num2str(smoothkern)]);
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
    if true || ~exist(label,'file') || optInputs(varargin,'overwrite')
        labeldir = [strrep(label,'.','_') '/'];
        if ~exist(labeldir,'dir');
            mkdir(labeldir);
        end
        
        unix(['rm -f ' labeldir '*']);
        
        %         if optInputs(varargin,'overwrite')
        %             unix(['rm -f ' labeldir '*']);
        %         end
        
        labelbase = [labeldir 'sigclusters'];
        unix_freesurfer_version(freesurfer_version, ['mri_surfcluster --in ' sigmap_mask_bin_smooth_bin ' --olab ' labelbase ' --subject fsaverage --hemi ' hemi ' --thmin 0.5 ' ]);
        
        % combine labels if necessary
        labels = mydir(labeldir);
        if length(labels) > 1
            unix_freesurfer_version(freesurfer_version, ['mri_mergelabels ' sprintf([' -i ' labeldir '%s'],labels{:}) ' -o ' label]);
        else
            unix(['cp -f ' labeldir labels{1} ' ' label]);
        end
        
    end
    
    
    if optInputs(varargin, 'tksurfer')
        unix_freesurfer_version(freesurfer_version, ['tksurfer fsaverage ' hemi ' inflated -overlay ' sigmap_mask ' -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' -label ' label ' -label-outline &']);
    else
        try
            freeview(sigmap_labelmask,hemi,[fminmax(1), mean(fminmax), fminmax(2)],'label',label);
        catch
            keyboard;
        end
    end
end

%%

% voxthresh_str = strrep(sprintf('%2.1f',voxthresh), '.', '');
% fprintf(['tksurfer myfsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.' contrastsign '.sig.masked.mgh -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' &\n']);
%
% % raw significance map
% unix_freesurfer_version(freesurfer_version, ['tksurfer myfsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/sig.mgh -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' &']);
%
% % thresholded significance map
% unix_freesurfer_version(freesurfer_version, ['tksurfer myfsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.' contrastsign '.sig.masked.mgh -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' &']);
%
% % significant clusters
% unix_freesurfer_version(freesurfer_version, ['tksurfer myfsaverage ' hemi ' inflated -annotation ' groupdir 'osgm/cache.th' voxthresh_str '.' contrastsign '.sig.ocn.annot &']);

% fprintf(['tksurfer myfsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.' contrastsign '.sig.ocn.mgh -fminmax 0.5 5.5 &\n']);
% unix_freesurfer_version(freesurfer_version, ['tksurfer myfsaverage ' hemi ' inflated -overlay ' groupdir 'osgm/cache.th' voxthresh_str '.contrastsign '.sig.ocn.mgh -fminmax 0.5 5.5 &']);
% unix_freesurfer_version(freesurfer_version, ['tksurfer myfsaverage ' hemi ' inflated -annot ' groupdir 'osgm/cache.th' voxthresh_str '.' contrastsign '.sig.ocn.annot &']);

%
%     cd(['~/' expsub{i,1} '/scripts/']);
%     allruns = read_runs(expsub{i,2}, runtype);
%     freesurfreg = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/anat2std.register.dat'];
%
%     if ~exist(freesurfreg,'file') || optInputs(varargin,'overwrite')
%         highres = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/highres.nii.gz'];
%         standard = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/standard.nii.gz'];
%         fslreg = ['~/' expsub{i,1} '/analysis/fla/sub' num2str(expsub{i,2}) '/main_r' num2str(allruns(1)) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/reg/highres2standard.mat'];
%         unix(['tkregister2 --mov ' highres ' --targ ' standard ' --fsl ' fslreg ' --reg ' freesurfreg ' --noedit ']);
%     end