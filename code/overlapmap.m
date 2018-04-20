function overlapmap(usubs,roi,loc,surface_fwhm,pthresh,varargin)

%% Paramaters

% string to identify this particular analysis
% idstring = [hemi '_' con '_usubs' sprintf('%d',sort(cat(1,expsub{:,2})))  '_fwhm' num2str(fwhm) '_' glmtype '_projfrac' num2str(projfrac(1)) '-' num2str(projfrac(2)) '_trilinear_' runtype '_' model  '_' num2str(fwhm*100,'%.0f') 'mm' ];

% output glm directory
outputdir = ['~/freesurfer/fsaverage/overlapmap/usubs' sprintf('-%d',usubs) '/'];
if ~exist(outputdir,'dir');
    mkdir(outputdir);
end

% smoothing
fwhm = read_smooth(loc(1).exp, usubs(1), loc(1).runtype, varargin{:});

% roi id string
roi_idstring = [roi.name '_' cat(2,roi.hemis{:}) ];

% id string for localizer analysis
if ~isempty(loc)
    loc_idstring = [];
    for i = 1:length(loc)
        loc_idstring = [loc_idstring 'loc' num2str(i) '_' loc(i).con  '_' DataHash([loc(i).exp loc(i).runtype loc(i).sign loc(i).op sprintf('%g',loc(i).opval) loc(i).model num2str(surface_fwhm) 'mm'])]; %#ok<AGROW>
    end
else
    loc_idstring = 'anatomical';
end

%%

nverts_hemi = 163842;
pmap = zeros(nverts_hemi*length(roi.hemis),1);
for q = 1:length(usubs)
    
    vertarea = nan(nverts_hemi*length(roi.hemis),1);
    for i = 1:length(roi.hemis)
        fid = fopen(['~/freesurfer/fsaverage/surf/' roi.hemis{i} '.area.asc'],'r');
        vertarea((1:nverts_hemi) + (i-1)*nverts_hemi) = fscanf(fid,'%f');
        fclose(fid);
    end
    
    mask = zeros(nverts_hemi*length(roi.hemis),1);
    for i = 1:length(roi.hemis)
        roi_brain = read_label(['~/freesurfer/fsaverage/label/' roi.hemis{i} '.' roi.name '.label']);
        mask(roi_brain.vnums + (i-1)*nverts_hemi + 1) = 1;
    end
    fprintf('Mask size: %d vertices, %.2f mm^2\n',sum(mask),sum(vertarea(logical(mask))));
        
    for i = 1:length(loc)
        
        loc(i).runs = read_runs(loc(i).exp, usubs(q), loc(i).runtype, varargin{:});
        idstring = [loc(i).con '_fwhm' num2str(surface_fwhm) '_projfrac0-1_trilinear_' loc(i).runtype '_' loc(i).model  '_' num2str(fwhm*100,'%.0f') 'mm'];
        
        locmask = zeros(nverts_hemi*length(roi.hemis),1);
        for j = 1:length(roi.hemis)
            if strcmp(loc(i).con(1:5), 'anova')
                sigmap = ['~/freesurfer/fsaverage/sla/' loc(i).exp '_us' num2str(usubs(q)) '/' roi.hemis{j} '.' loc(i).con '_fwhm' num2str(surface_fwhm) '_projfrac0-1_trilinear_' loc(i).runtype '_r' sprintf('%d',loc(i).runs) '_' loc(i).model  '_' num2str(fwhm*100,'%.0f') 'mm.mgz'];
                if ~exist(sigmap,'file') || optInputs(varargin,'overwrite')
                    anova_vol = [params('rootdir') loc(i).exp '/analysis/sla/usub' num2str(usubs(q)) '/' loc(i).runtype '_' loc(i).con '_r' sprintf('%d',loc(i).runs) '_' loc(i).model '_' num2str(fwhm*100,'%.0f') 'mm.gfeat/cope1.feat/stats/zfstat1.nii.gz'];
                    regfile = regfile_stdLAS2anatRAS(loc(i).exp, usubs(q), loc(i).model, varargin{:});
                    refvol = ['~/freesurfer/' loc(i).exp '_us' num2str(usubs(q)) '/mri/orig.mgz'];
                    unix_freesurfer(['mri_vol2surf --ref ' refvol ' --mov ' anova_vol ' --reg ' regfile ' --o ' sigmap ' --hemi ' roi.hemis{j} ' --trgsubject fsaverage --surf white --interp trilinear --projfrac-avg 0 1 0.05 '  ifelse(surface_fwhm > 0, [' --surf-fwhm ' num2str(surface_fwhm) ' '], '')]);
                end
            else
                sigmap = ['~/freesurfer/fsaverage/sla/' loc(i).exp '_us' num2str(usubs(q)) '/' roi.hemis{j} '_' idstring '_r' sprintf('%d',loc(i).runs) '/osgm/sig.mgh'];
                if ~exist(sigmap,'file')
                    sla_surf(loc(i).exp,usubs(q),loc(i).con,surface_fwhm,[],[],[],[],roi.hemis{j},loc(i).runtype,loc(i).model,'nocorrection',varargin{:},'fsaverage','runs', loc(i).runs);
                end
            end
            
            % read files
            x = MRIread(sigmap);
            locmask((1:nverts_hemi) + (j-1)*nverts_hemi) = x.vol;
        end
        
        % flip contrast if necessary
        switch loc(i).sign
            case 'pos'
            case 'neg'
                locmask = -locmask; %#ok<*AGROW>
            case 'abs'
                locmask = abs(locmask);
            otherwise
                error('Bad sign');
        end
        
        locmask(~logical(mask)) = 0;
        mask = brainops(locmask,loc(i).op,loc(i).opval,sum(mask));
        fprintf('Mask size: %d vertices, %.2f mm^2\n',sum(mask),sum(vertarea(logical(mask))));
    end
    
    pmap(logical(mask)) = pmap(logical(mask)) + 1/length(usubs);
    
end

%%

for i = 1:length(roi.hemis)
    ovmap_file = [outputdir '/' roi.hemis{i} '.ovmap_' roi_idstring '_' loc_idstring '_pthresh' sprintf('%.2f',pthresh) '.mgh'];
    ovmap = MRIread(['~/freesurfer/fsaverage/surf/' roi.hemis{i} '.area.mgh']);
    ovmap.vol = pmap((1:nverts_hemi)+(i-1)*nverts_hemi);
    ovmap.vol(ovmap.vol < pthresh - 1e-3) = 0;
    ovmap.fspec = ovmap_file;
    MRIwrite(ovmap, ovmap_file);
    
    pmax = 0.5;
    %     inputrange = [0, pthresh, pmax, 1];
    %     colormapping = [0, pthresh + 0.3*(pmax-pthresh), pmax, 1];
    %
    %     ovmap_file_colorwarp = strrep(ovmap_file,'.mgh','_colorwarp.mgh');
    %     ovmap.fspec = ovmap_file_colorwarp;
    %     xi = ovmap.vol > inputrange(1) & ovmap.vol < inputrange(end);
    %     ovmap.vol(xi) = interp1(inputrange,colormapping,ovmap.vol(xi));
    %     MRIwrite(ovmap, ovmap_file_colorwarp);
    
    if optInputs(varargin, 'freeview')
        freeview(ovmap_file, roi.hemis{i}, [0.01, pthresh, pmax]);
    end
end


%%

addpath('/software/Freesurfer/emf-5.1.0/matlab');

% freeview(sigmap_mask_bin, roi.hemis{i}, [0.1 0.25 0.5]);


if optInputs(varargin, 'label')
    % binarize
    
    for i = 1:length(roi.hemis)
        sigmap = [outputdir '/' roi.hemis{i} '.ovmap_' roi_idstring '_' loc_idstring '_pthresh' sprintf('%.2f',pthresh) '.mgh'];
        sigmap_mask_bin = strrep(sigmap,'.mgh','.bin.mgh');
        if ~exist(sigmap_mask_bin, 'file') || optInputs(varargin,'overwrite')
            [srf, M, mr_parms] = load_mgh(sigmap);
            srf(srf < pthresh - 1e-3) = 0;
            srf(srf > pthresh - 1e-3) = 1;
            save_mgh(srf,sigmap_mask_bin,M,mr_parms);
        end
        
        % optionally smooth binarized map and rebinarize
        if optInputs(varargin,'smoothbin')
            
            smoothkern = varargin{optInputs(varargin,'smoothbin') + 1};
            sigmap_mask_bin_smooth = strrep(sigmap_mask_bin,'.mgh',['.smooth' num2str(smoothkern) '.mgh']);
            if ~exist(sigmap_mask_bin_smooth,'file') || optInputs(varargin,'overwrite')
                unix_freesurfer(['mri_surf2surf --s fsaverage --sval ' sigmap_mask_bin ' --tval ' sigmap_mask_bin_smooth  ' --hemi ' roi.hemis{i} ' --fwhm-src ' num2str(smoothkern)]);
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
            unix_freesurfer(['mri_surfcluster --in ' sigmap_mask_bin_smooth_bin ' --olab ' labelbase ' --subject fsaverage --hemi ' roi.hemis{i} ' --thmin 0.5 ' ]);
            
            % combine labels if necessary
            labels = mydir(labeldir);
            if length(labels) > 1
                unix_freesurfer(['mri_mergelabels ' sprintf([' -i ' labeldir '%s'],labels{:}) ' -o ' label]);
            else
                unix(['cp -f ' labeldir labels{1} ' ' label]);
            end
            
        end
        
        if optInputs(varargin, 'tksurfer')
            unix_freesurfer(['tksurfer fsaverage ' roi.hemis{i} ' inflated -overlay ' sigmap_mask ' -fminmax ' num2str(fminmax(1)) ' ' num2str(fminmax(2)) ' -label ' label ' -label-outline &']);
        else
            try
                freeview(sigmap,roi.hemis{i},[0.01, pthresh, pmax],'label',label);
            catch
                keyboard;
            end
        end
    end
end


% %% Preprocesing
%
%
% for i = 1:length(usubs)
%
%     % freesurfer subject id
%     subjid = [expsub{i,1} '_us' num2str(expsub{i,2})];
%
%     % runs for that subject
%     allruns = read_runs(expsub{i,1}, expsub{i,2}, runtype, varargin{:});
%
%     % second level copes and varcopes
%     x = [con '_fwhm' num2str(fwhm) '_projfrac0-1_trilinear_' runtype '_' model  '_' num2str(fwhm*100,'%.0f') 'mm'];
%     pmap = MRIread(['~/freesurfer/fsaverage/sla/' subjid '/' hemi '_' x '_r' sprintf('%d',allruns) '/osgm/sig.mgh']);
%
%     switch con_sign
%         case 'pos'
%         case 'neg'
%             pmap.vol = -pmap.vol; %#ok<*AGROW>
%         case 'abs'
%             pmap.vol = abs(pmap.vol);
%         otherwise
%             error('Bad sign');
%     end
%
%     try
%     ovmap.vol(pmap.vol > voxthresh) = ovmap.vol(pmap.vol > voxthresh) + 1;
%     catch
%         keyboard
%     end
% end
%
% ovmap.vol = ovmap.vol/length(expsub);
%
% %%



%
% allcopes = [outputdir  'allcopes.mgh'];
% allvarcopes = [outputdir  'allvarcopes.mgh'];
%
% if fwhm > 0;
%     fwhm_str = ['--fwhm ' num2str(fwhm)];
% else
%     fwhm_str = '';
% end
%
% if ~exist(allcopes,'file') || ~exist(allvarcopes,'file') || optInputs(varargin,'overwrite')
%     unix_freesurfer(['mris_preproc --target fsaverage --hemi ' hemi ' --out ' allcopes ' ' copestr ' ' fwhm_str]);
%     unix_freesurfer(['mris_preproc --target fsaverage --hemi ' hemi ' --out ' allvarcopes ' ' varstr ' ' fwhm_str]);
% end
