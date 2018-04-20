function psc = roi_psc(us, roi, test, loc, varargin)

% [psc testruns_used main_idstring] = roi_psc(s, testruntype, testmodel, masktype, roi, varargin)
% For anatomical + contrast: roi_psc(s, testruntype, testmodel, masktype, roi, locruntype, locmodel, contrast)
%
% Main function for computing the psc value for an roi. If loc and test runs are overlapping, uses leave-one-out style analysis.
% Specificy a subject, a test runtype (e.g. 'main', 'loc'), a test model (e.g. 'block' or 'event'), and then a type of mask: anat, anat_contrast, and cluster
%
% For anat you only need to specify an roi (needs to be updated)
% For cluster and anat_contrast, you need to specify an anatomical roi, a localizer contrast, a localizer runtype, & a localizer model (next two arguments respectively)
%
% Sample
% [psc testruns_used main_idstring] = roi_psc(1, 'main', 'block', 'anat_contrast', 'destrieux_pp', 'harm-highfreq_vs_noise-highfreq', 'localizer', 'block')
% roi_psc(1, 'main', 'block', 'anat', 'hvprob_pp', 'harm_vs_noise', 'main', 'block')

%% Directories, Parameters

roidir = [params('rootdir') test.exp '/analysis/roi/'];
if ~exist(roidir,'dir');
    mkdir(roidir);
end

% id string for this roi
x = sprintf('%g-',roi.opval);
roi_idstring = [roi.name '_' roi.op strrep(x(1:end-1),'.','p')];

allconds = read_conditions(test.exp, us, test.runtype, varargin{:});
% id string for the test run
test_idstring = [test.exp '_' test.runtype '_r' num2str(test.run) '_' test.model '_'  num2str(100*test.fwhm, '%.0f') 'mm' '_conds' sprintf('%d',find(ismember(allconds,test.conds)))];

% id string for localizer analysis
if ~isempty(loc)
    loc_idstring = [];
    for i = 1:length(loc)
        x = sprintf('%g-',loc(i).opval);
        loc_idstring = [loc_idstring 'loc' num2str(i) '_' loc(i).exp '_' loc(i).runtype '_' loc(i).con '_' loc(i).sign  '_' loc(i).op '_' strrep(x(1:end-1),'.','p') '_r' sprintf('%d',loc(i).runs) '_' loc(i).model '_' num2str(100*loc(i).fwhm, '%.0f') 'mm']; %#ok<AGROW>
    end
    if length(loc)>1
        loc_idstring = DataHash(loc_idstring);
    end
    loc_idstring = DataHash(loc_idstring);
else
    loc_idstring = 'anatomical';
end

if optInputs(varargin, 'sort_surface_axis');
    percrange = varargin{optInputs(varargin, 'sort_surface_axis')+1};
    loc_idstring = [loc_idstring '_sort_surface_axis' num2str(percrange(1)) '-' num2str(percrange(2))];
end

% second ROI to remove from the first roi
% roidisjoin = '';
% if optInputs(varargin,'roidisjoin');
%     roidisjoin = varargin{optInputs(varargin,'roidisjoin')+1};
%     roi_idstring = [roi_idstring '_disjoin_' roidisjoin];
% end

%% Registration

% Registration Algorithm Used

% types of registration to move between
% functional, highres, and standard spaces.
% only used for registering anatomical masks.

% defaults are bbregister to go from functional
% to highres (flirt is alternative)
% and flirt to go from highres to standard (fnirt is alternative).

func2highres = 'bbreg';
if optInputs(varargin,'func2highres');
    func2highres = varargin{optInputs(varargin,'func2highres')+1};
end

highres2standard = 'flirt';
if optInputs(varargin,'highres2standard');
    highres2standard = varargin{optInputs(varargin,'highres2standard')+1};
end

% if freesurfer roi is being used
% then there is no highres to standard transformation
% and so this is just set to 'freesurfer'
if ~isempty(findstr('destrieux_',roi.name)) || ~isempty(findstr('mylabel_',roi.name)) % destrieux is from freesurfer and so doesnt really have a standard to highres transformation
    highres2standard = 'freesurfer';
end

%% Signal Averaging Parameters if Used instead of Betas

sigav = false;
win = [NaN, NaN];
stat_idstring = 'betas';
if optInputs(varargin,'sigav_timecourse') || optInputs(varargin,'sigav_plateau') || optInputs(varargin,'sigav_peaktime');
    
    sigav = true;
    [blockdur, nulldur, TR, TA, stimdur, stim2scan, win, ~,~,nullwin] = read_scanparams(test.exp, us, test.runtype); %#ok<ASGLU>
    
    % window
    stat_idstring = ['sigav_win' strrep(num2str(win(1)),'.','p') '-' strrep(num2str(win(end)),'.','p')];
    
    if isempty(nullwin) && optInputs(varargin,'nulltc')
        nullwin = win;
        stat_idstring = [stat_idstring '_nulltc'];
    elseif isempty(nullwin)
        nullwin = (TR:TR:nulldur-TR) + stim2scan;
        %         nullwin = (TR:TR:nulldur-TR) + stim2scan;
    end
end

%% Misc Setup

% string used to remember this specific analysis
main_idstring = ['psc_us' num2str(us) '_' stat_idstring '_' roi_idstring '_' test_idstring '_' loc_idstring];
if length(main_idstring) > 400
    main_idstring = DataHash(main_idstring);
end

% initialize psc matrix
if sigav
    psc = nan(length(win), length(test.conds));
else
    psc = nan(length(test.conds),1);
end

if ~optInputs(varargin,'nowrite')
    fprintf('\nS%d, %s\n',us, roi.name);
    for i = 1:length(loc)
        fprintf('Localizer: %s, runs %s\n',loc(i).con, sprintf('%d ',loc(i).runs))
        fprintf('Test: run %d\n', test.run);
    end
end

if exist([roidir main_idstring '.mat'],'file') && ~optInputs(varargin,'overwrite')
    load([roidir main_idstring '.mat']);
    return;
end

%% Read in ROI and Localizer

testfeat = [params('rootdir') test.exp '/analysis/fla/usub' num2str(us) '/' test.runtype '_r' num2str(test.run) '_' test.model  '_' num2str(100*test.fwhm, '%.0f') 'mm.feat/'];

% keyboard;
if length(roi.name) >= 8 && strcmp(roi.name(1:8),'handpick')
    
    error('Needs to be updated');
    % cluster file, register to testrun if not already
    register_cluster(roi,us,testruntype,r,testmodel,loc,locruntype,locruns,locmodel,varargin{:});
    
    clustfile = [testfeat 'masks/' roi '_' locruntype '_' loc '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm_z' strrep(num2str(zthresh),'.','') '_func_func2highres_' func2highres '.nii.gz'];
    if ~exist(clustfile,'file');
        psc = NaN;
        fprintf('No cluster for this run.\n'); drawnow;
        return;
    end
    roi_brain = readmr(clustfile,'NOPROGRESSBAR');
    
elseif ~isempty(strfind(roi.name, 'func_cluster'))
    
    anatfile = register_func_cluster(roi.name, us, test, varargin{:});
    roi_brain = readmr(anatfile,'NOPROGRESSBAR');
    
else
    
    anatfile = [testfeat 'masks/' roi.name  '_func_func2highres_' func2highres  '_highres2standard_' highres2standard  '.nii.gz'];
    roi_brain = readmr(anatfile,'NOPROGRESSBAR');
    
end

if ~isempty(loc)
    
    clear loc_brain;
    
    for i = 1:length(loc)
        if optInputs(varargin, 'sort_surface_axis')
            switch roi.name
                case 'mylabel_stp'
                    roitemp.name = 'stp';
                    roitemp.hemis = {'rh','lh'};
                otherwise
                    error('No valid ROI');
            end
            zstatfile = contrast_sort_surface_axis(us, roitemp, test, loc, percrange, varargin{:});
            loc_brain(i) = readmr(zstatfile,'NOPROGRESSBAR');
            

        else
            
            if strcmp(loc(i).con(1:10), 'sigav-corr')
                register_sigav_corr(us,test,loc(i),varargin{:});
            else
                % register contrast, will not overwrite unless overwrite is specified in optional arguments
                register_contrast(us,test,loc(i),varargin{:});
            end
            
            % read filesif length(runnum) > 100
            zstatfile = [testfeat 'slacontrasts/' loc(i).exp '_' loc(i).runtype '_' loc(i).con '_r' sprintf('%d',loc(i).runs) '_' loc(i).model '_' num2str(100*loc(i).fwhm, '%.0f') 'mm.nii.gz'];
            if length(loc(i).runs) > 100
                zstatfile = strrep(zstatfile, ['_r' sprintf('%d',loc(i).runs)], ['_r' DataHash(loc(i).runs)]);
            end
            loc_brain(i) = readmr(zstatfile,'NOPROGRESSBAR');
            
            % flip contrast if necessary
            switch loc(i).sign
                case 'pos'
                case 'neg'
                    loc_brain(i).data = -loc_brain(i).data; %#ok<*AGROW>
                case 'abs'
                    loc_brain(i).data = abs(loc_brain(i).data);
                otherwise
                    error('Bad sign');
            end
        end
    end
end

%% Create Mask

% roi with operation applied, should probably be a binary operation, such as thresholding
mask = roi_brain;
mask.data = brainops(mask.data, roi.op, roi.opval);
masksize = sum(mask.data(:));
fprintf('Mask size: %d voxels, %.2f cubic mm\n',masksize,masksize*2.0833*2.0833*4.4);

% remove a second roi from the initial roi, should recheck code before using
% if ~isempty(roidisjoin)
%     error('Check this functionality before using.');
%     anatfile = [testfeat 'masks/' roidisjoin  '_func_func2highres_' func2highres  '_highres2standard_' highres2standard  '.nii.gz']; %#ok<UNRCH>
%     x = readmr(anatfile,'NOPROGRESSBAR');
%     y = brainops(x.data, roi_op, roi_opval);
%     mask.data(y==1) = 0;
% end

if ~isempty(loc)
    for i = 1:length(loc)
        if optInputs(varargin, 'sort_surface_axis')
            mask.data(loc_brain(i).data < 0.5) = 0;
        else
            locmask = loc_brain(i).data;
            locmask(~logical(mask.data)) = 0;
            mask.data = brainops(locmask,loc(i).op,loc(i).opval,sum(mask.data(:)));
        end
        fprintf('Mask size: %d voxels, %.2f cubic mm\n',sum(mask.data(:)),sum(mask.data(:))*2.0833*2.0833*4.4);
    end
end

% if ~isempty(loc2)
%     roimask = logical(mask.data);
%     mask.data = loc2_brain.data;
%     mask.data(~roimask) = -inf;
%     mask.data = brainops(mask.data,loc2_op,loc2_opval,masksize);
%     masksize = sum(mask.data(:));
%     fprintf('Mask size: %d voxels, %.2f cubic mm\n',masksize,masksize*2.0833*2.0833*4.4);
% end

%% Calculate PSC Values

% normalization factor
masknorm = sumdims(mask.data,1:3);

% 3) Calculate PSC values, either by signal averaging or using betas
if sigav
        
    fprintf('Reading functional data...\n');
    
    func_file = [testfeat 'filtered_func_data.nii.gz'];
    if ~exist(func_file,'file');
        func_file = [params('rootdir') test.exp '/analysis/preprocess/usub' num2str(us) '/' test.runtype '_r' num2str(test.run) '/smooth' num2str(100*test.fwhm, '%.0f') 'mm.nii.gz'];
    end
    if ~exist(func_file, 'file')
        fprintf('Error in roi_psc: cannot find file %s\n', func_file);
        drawnow;
        keyboard;
    end
    funcbr = readmr( func_file , 'NOPROGRESSBAR' );
    
    %     funcbr = readmr( preproc_file, 'NOPROGRESSBAR' );
    funcmask = funcbr.data .* repmat(mask.data, [1 1 1 size(funcbr.data,4)]);
    
    % fixation based on response during null periods
    behavdat = read_timings(test.exp, us, test.runtype, test.run, test.model);
    x = behavdat.onsets(strcmp('NULL',behavdat.conds));
    y = repmat(x, 1, length(nullwin)) + repmat(nullwin, length(x), 1);
    z = round(y/TR) + 1;
    if any(abs(z(:) - (y(:)/TR + 1)) > 1e-3)
        fprintf('Error in roi_psc.m: index should be an integer\n');
        drawnow;
        keyboard;
    end
    
    if optInputs(varargin,'nulltc')
        onsets = behavdat.onsets(strcmp('NULL',behavdat.conds));
        % loop through each stimulus presentation/repetition
        nullmean_allreps = nan( length(nullwin), length(onsets) );
        for j = 1:length(onsets)
            x = round((onsets(j) + nullwin)/TR) + 1;%(onsets(j)/TR + 1) + (win(1)/TR : win(2)/TR);
            try
                y = squeeze(sumdims(funcmask(:,:,:,x), 1:3 )) / masknorm;
            catch
                keyboard
            end
            if any(isnan(y));
                error('signal average is nan for some reason');
            end
            nullmean_allreps(:,j) = y;
        end
        nullmean = nanmean(nullmean_allreps,2);
    else
        nullmean = squeeze(sumdims(funcmask(:,:,:,z(:)),1:4)) / (numel(z)*masknorm);
    end
    
    for i = 1:length(test.conds)
        
        onsets = behavdat.onsets(strcmp(test.conds{i},behavdat.conds));
        
        % loop through each stimulus presentation/repetition
        psc_allreps = nan( length(win), length(onsets) );
        for j = 1:length(onsets)
            x = round((onsets(j) + win)/TR) + 1;%(onsets(j)/TR + 1) + (win(1)/TR : win(2)/TR);
            if abs(x - ((onsets(j) + win)/TR + 1)) > 1e-3
                fprintf('Error in roi_psc.m: index should be an integer\n');
                drawnow;
                keyboar;
            end
            y = squeeze(sumdims(funcmask(:,:,:,x), 1:3 )) / masknorm;
            if any(isnan(y));
                error('signal average is nan for some reason');
            end
            psc_allreps(:,j) = 100*((y-nullmean)./nullmean);
        end
        psc(:,i) = nanmean(psc_allreps,2);
        
    end
    
else
    
    
    % 2) Read in Mean Func
    nullbr = readmr([testfeat 'mean_func.nii.gz'],'NOPROGRESSBAR');
    nullmask = mask.data.*nullbr.data;
    nullmean = squeeze(sumdims(nullmask,1:3)/masknorm);
    if isnan(nullmean);
        fprintf(['fixation is nan.\n']); drawnow;
        keyboard;
    end
    
    
    % read in names of evs in the proper order
    fid = fopen([testfeat 'evname.txt'],'r');
    tmp = textscan(fid,'%s'); fclose(fid);
    evname = tmp{1};
    evs = test.conds;
    
    for i = 1:length(evs)
        
        x = find(strcmp(evs{i},evname));
        
        % pe is basically multiplied by 2 if temporal derivatives are used
        if strcmp(params('tempderiv',test.model,test.exp), '1');
            pe = (x-1)*2+1;
        else
            pe = x;
        end
        
        % weighted pe for mask
        pefile = [testfeat 'stats/pe' num2str(pe) '.nii.gz'];
        pebr = readmr(pefile, 'NOPROGRESSBAR');
        pemask = mask.data.*pebr.data;
        pestat = squeeze(sumdims(pemask,1:3))/masknorm;
        
        if isnan(pestat);
            error('pe is nan for some reason');
        end
        
        % normalize to psc
        psc(i) = 100*(pestat/nullmean);
        
    end
    
end

save([roidir main_idstring '.mat'],'psc');

%% Scraps
%
%
%     switch masktype
%
%         case {'anat'}
%
%             m2 = m1;
%
%         case {'anat_contrast', 'cluster'}
%
%             zstatfile = [testfeat 'slacontrasts/' locruntype '_' loc '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm.nii.gz'];
%             % register contrast, will not overwrite unless overwrite is specified in optional arguments
%             register_contrast(s,testruntype,r,testmodel,loc,locruntype,locruns,locmodel,varargin{:});
%
%             % read files
%             anatbr = readmr(anatfile,'NOPROGRESSBAR');
%             zstatbr = readmr(zstatfile,'NOPROGRESSBAR');
%
%             % flip contrast if necessary
%             if negcontrast;
%                 zstatbr.data = -zstatbr.data;
%             end
%
%             if ~isempty(ztopvoxel)
%
%                 m1inds = find(m1.data);
%                 [~,ztopinds] = zstatbr.data(m1inds);
%                 m2 = m1;
%                 m2.data = zeros(size(m1.data));
%                 m2.data(m1inds(ztopinds(1:ztopvoxels))) = 1;
%
%             else
%
%                 % weight zstatistic and truncate extreme values
%                 zpow = zstatbr.data.^zpower;
%                 zpow(zpow > 1e10) = 1e10;
%                 m2 = m1;
%                 m2.data = m1.data .* (zstatbr.data > zthresh) .* zpow;
%             end
%
%     end
%
%             end
%
%
%
%             if strcmp(masktype, 'anat')
%
%                 m2 = anatmaskbr;
%
%             elseif strcmp(masktype, 'anat_contrast');
%
%                 zstatfile = [testfeat 'slacontrasts/' locruntype '_' loc '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm.nii.gz'];
%
%                 % register contrast, will not overwrite unless overwrite is specified in optional arguments
%                 register_contrast(s,testruntype,r,testmodel,loc,locruntype,locruns,locmodel,varargin{:});
%
%                 % read files
%                 anatbr = readmr(anatfile,'NOPROGRESSBAR');
%                 zstatbr = readmr(zstatfile,'NOPROGRESSBAR');
%
%                 % flip contrast if necessary
%                 if negcontrast;
%                     zstatbr.data = -zstatbr.data;
%                 end
%
%                 if ~
%
%
%         case 'cluster'
%
%             % if test and localizer runtypes are the same, use left-over runs to localizer
%             % else use all of the localizer runs
%             if strcmp( testruntype, locruntype );
%                 locruns = setdiff(testruns,r);
%             else
%                 locruns = read_runs(s,locruntype);
%             end
%
%             fprintf('Localizer: %s, runs %s\n', locruntype, sprintf('%d ',locruns));
%
%             % contrast and anatomical file
%             clustfile = [testfeat 'masks/' roi '_' locruntype '_' loc '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm_z' strrep(num2str(zthresh),'.','') '_func_func2highres_' func2highres '.nii.gz'];
%             zstatfile = [testfeat 'slacontrasts/' locruntype '_' loc '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm.nii.gz'];
%
%             % register contrast, will not overwrite unless overwrite is specified in optional arguments
%             register_contrast(s,testruntype,r,testmodel,localizer,locruntype,locruns,locmodel,varargin{:});
%
%             % register contrast, will not overwrite unless overwrite is specified in optional arguments
%             register_cluster(roi,s,testruntype,r,testmodel,localizer,locruntype,locruns,locmodel,varargin{:});
%
%             if ~exist(clustfile,'file');
%                 fprintf('No cluster for this run.\n'); drawnow;
%                 continue;
%             end
%
%             % read files
%             clustbr = readmr(clustfile,'NOPROGRESSBAR');
%             zstatbr = readmr(zstatfile,'NOPROGRESSBAR');
%
%             % flip contrast if necessary
%             if negcontrast;
%                 zstatbr.data = -zstatbr.data;
%             end
%
%             % weight zstatistic and truncate extreme values
%             zpow = zstatbr.data.^zpower;
%             zpow(zpow > 1e10) = 1e10;
%
%             % weight anatomical mask and truncate extreme values
%             mpow = clustbr.data.^maskpower;
%             mpow(mpow > 1e10) = 1e10;
%
%             % threshold and weight
%             m2 = clustbr;
%             maskbr.data = mpow .* (zstatbr.data>zthresh) .* zpow;
%
%     end
%
%
%
%
%
%         case 'anat_contrast'  % anatomical ANDed with functional contrast
%
%             % if test and localizer runtypes are the same, use left-over runs to localizer
%             % else use all of the localizer runs
%             if strcmp( testruntype, locruntype );
%                 locruns = setdiff(testruns,r);
%             else
%                 locruns = read_runs(s,locruntype);
%             end
%
%             fprintf('Localizer: %s, runs %s\n', locruntype, sprintf('%d ',locruns));
%
%             % contrast and anatomical file
%             anatfile = [testfeat 'masks/' roi  '_func_func2highres_' func2highres  '_highres2standard_' highres2standard  '.nii.gz'];
%             zstatfile = [testfeat 'slacontrasts/' locruntype '_' localizer '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm.nii.gz'];
%
%             % register contrast, will not overwrite unless overwrite is specified in optional arguments
%             register_contrast(s,testruntype,r,testmodel,localizer,locruntype,locruns,locmodel,varargin{:});
%
%             % read files
%             anatbr = readmr(anatfile,'NOPROGRESSBAR');
%             zstatbr = readmr(zstatfile,'NOPROGRESSBAR');
%
%             % flip contrast if necessary
%             if negcontrast;
%                 zstatbr.data = -zstatbr.data;
%             end
%
%
%             % weight zstatistic and truncate extreme values
%             zpow = zstatbr.data.^zpower;
%             zpow(zpow > 1e10) = 1e10;
%
%             % weight anatomical mask and truncate extreme values
%             mpow = anatbr.data.^maskpower;
%             mpow(mpow > 1e10) = 1e10;
%
%             if ~isempty(masksize) % threshold anatomical by taking the best "masksize" number of voxels
%                 nvox = round(masksize/prod(funcdims));
%                 x = mpow(:);
%                 [maskvals, maskinds] = sort(x,'descend');
%                 if maskvals(nvox) == 0;
%                     error('Mask nvox is Too Big');
%                 end
%                 threshinds = zeros(size(maskinds));
%                 threshinds(maskinds(1:nvox)) = 1;
%                 mthresh = reshape(threshinds, size(anatbr.data));
%
%                 maskbr = anatbr;
%                 maskbr.data = mthresh .* (zstatbr.data>zthresh) .* zpow;
%             elseif ~isempty(maskthresh) % threshold based on some absolute value
%                 maskbr = anatbr;
%                 maskbr.data = (mpow > maskthresh).* (zstatbr.data>zthresh) .* zpow;
%             else % weighted mask
%                 maskbr = anatbr;
%                 maskbr.data = mpow .* (zstatbr.data>zthresh) .* zpow;
%             end
%
%         case 'cluster'
%
%             % if test and localizer runtypes are the same, use left-over runs to localizer
%             % else use all of the localizer runs
%             if strcmp( testruntype, locruntype );
%                 locruns = setdiff(testruns,r);
%             else
%                 locruns = read_runs(s,locruntype);
%             end
%
%             fprintf('Localizer: %s, runs %s\n', locruntype, sprintf('%d ',locruns));
%
%             % contrast and anatomical file
%             clustfile = [testfeat 'masks/' roi '_' locruntype '_' localizer '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm_z' strrep(num2str(zthresh),'.','') '_func_func2highres_' func2highres '.nii.gz'];
%             zstatfile = [testfeat 'slacontrasts/' locruntype '_' localizer '_r' sprintf('%d',locruns) '_' locmodel '_' params('smooth') 'mm.nii.gz'];
%
%             % register contrast, will not overwrite unless overwrite is specified in optional arguments
%             register_contrast(s,testruntype,r,testmodel,localizer,locruntype,locruns,locmodel,varargin{:});
%
%             % register contrast, will not overwrite unless overwrite is specified in optional arguments
%             register_cluster(roi,s,testruntype,r,testmodel,localizer,locruntype,locruns,locmodel,varargin{:});
%
%             if ~exist(clustfile,'file');
%                 fprintf('No cluster for this run.\n'); drawnow;
%                 continue;
%             end
%
%             % read files
%             clustbr = readmr(clustfile,'NOPROGRESSBAR');
%             zstatbr = readmr(zstatfile,'NOPROGRESSBAR');
%
%             % flip contrast if necessary
%             if negcontrast;
%                 zstatbr.data = -zstatbr.data;
%             end
%
%             % weight zstatistic and truncate extreme values
%             zpow = zstatbr.data.^zpower;
%             zpow(zpow > 1e10) = 1e10;
%
%             % weight anatomical mask and truncate extreme values
%             mpow = clustbr.data.^maskpower;
%             mpow(mpow > 1e10) = 1e10;
%
%             % threshold and weight
%             maskbr = clustbr;
%             maskbr.data = mpow .* (zstatbr.data>zthresh) .* zpow;
%
%     end

% The default type of mask is just a binary mask with a set threshold (default = 0.5)
% For non-binary masks, masks can be binarized by picking the top.
% You can also just weight voxels by their mask weight, raised to a power.
% The power can be used to emphasize voxels that better fit the roi.
% You can only use 1 of the 3 methods, not combinations.
% % threshold mask with a specific number, default
% roithresh = 0.5;
% if optInputs(varargin,'roithresh');
%     roithresh = varargin{optInputs(varargin,'roithresh')+1};
% end
% roi_idstring = ['roithresh' strrep(num2str(roithresh), '.', 'p')];
%
% % threshold mask by picking top masksize voxels
% roisize = [];
% if optInputs(varargin,'roisize');
%     roisize = varargin{optInputs(varargin,'roisize')+1};
%     roithresh = [];
%     roi_idstring = ['roisize' num2str(roisize)];
% end
%
% % default is to linearly weight localizer voxels by mask weight
% roipower = [];
% if optInputs(varargin,'roipower');
%     roipower = varargin{optInputs(varargin,'roipower')+1};
%     roisize = [];
%     roithresh = [];
%     roi_idstring = ['roipower' num2str(roipower)];
% end
%
% roidisjoin = '';
% if optInputs('roidisjoin');
%     roidisjoin = varargin{ptInputs('roidisjoin')+1};
%     roi_idstring = [roi_idstring '_disjoin_' roidisjoin];
% end

%
%
% % voxels above threshold from the localizer
% zthresh = params('zthresh');
% if optInputs(varargin,'zthresh');
%     zthresh = varargin{optInputs(varargin,'zthresh')+1};
% end
% zstr = ['zthresh' strrep(num2str(zthresh), '.', 'p')];
%
% % weight voxels z-stat of the localizer
% zpower = [];
% if optInputs(varargin,'zpower');
%     zpower = varargin{optInputs(varargin,'zpower')+1};
%     zthresh = [];
%     zstr = ['zpower' num2str(zpower)];
% end
%
% % select ztopvoxels with the highest z-stat from the localizer
% % requires a set mask threshold
% zvoxsort = [];
% if optInputs(varargin,'zvoxsort');
%     zvoxsort = varargin{optInputs(varargin,'zvoxsort')+1};
%     zthresh = [];
%     zpower = [];
%     zstr = ['zvoxsort' num2str(zvoxsort(1)) '-' num2str(zvoxsort(2))];
%     if ~isempty(roipower)
%         error('Selecting Top Localizer Voxels Requires a Binary Mask');
%     end
% end

% dimensions of functional image
% funcdims = [2.083, 2.083, 4.400];

% string used to remember contrast, not relevant if a contrast isn't specified
% if isempty(localizer)
%     loc_idstring = '';
% else
%     loc_idstring = ['loc_' ifelse(negcontrast,'neg-','') localizer '_' locruntype '_r' sprintf('%d',locruns) '_' locmodel];
% end


% % 4) Threshold or Weight the Localizer Data
% if isempty(localizer)
%
%     m2 = m1;
%
% elseif ~isempty(zthresh)
%
%     % binarize initial mask
%     m2 = m1;
%     m2.data = m1.data .* (loc.data > zthresh);
%
% elseif ~isempty(zvoxsort)
%
%     % threshold anatomical by taking the best "masksize" number of voxels
%     m1inds = find(m1.data > 0.5);
%     [~, zsortinds] = sort(loc.data(m1inds),'descend');
%
%     m2 = m1;
%     m2.data = zeros(size(m1.data));
%     m2.data(m1inds(zsortinds(zvoxsort(1):zvoxsort(2)))) = 1;
%
% elseif ~isempty(zpower)
%
%     m2 = m1;
%     m2.data = m1.data .* (loc.data .^ zpower);
%     m2.data(m2.data > 1e10) = 1e10;
%
% end
%
% if sumdims(m2.data,1:3)==0;
%     psc = NaN;
%     fprintf('Mask is empty.\n');
%     return;
% end




% 2) Threshold or Weight the Initial ROI
% if ~isempty(roithresh)
%
%
%     % binarize initial mask
%     m1 = roibr;
%     m1.data = roibr.data > roithresh;
%
% elseif ~isempty(roisize)
%
%     % threshold anatomical by taking the best "masksize" number of voxels
%     n = round(roisize/prod(funcdims));
%     [sortvals, sortinds] = sort(roibr.data(:),'descend');
%
%     if sortvals(n) == 0;
%         error('Mask nvox is Too Big');
%     end
%
%     m1 = roibr;
%     m1.data = zeros(size(roibr.data));
%     m1.data(sortinds(1:n)) = 1;
%
% else
%
%     m1 = roibr;
%     m1.data = roibr.data.^roipower;
%     m1.data(m1.data > 1e10) = 1e10;
%
% end





