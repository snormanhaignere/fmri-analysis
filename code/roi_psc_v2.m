function [psc_func_roi] = roi_psc_v2(us, roi, test, loc, varargin)

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

% keyboard;

roidir = [params('rootdir') test.exp '/analysis/roi/'];
if ~exist(roidir,'dir');
    mkdir(roidir);
end

% id string for this roi
x = sprintf('%g-',roi.opval);
roi_idstring = [roi.name '_' roi.op strrep(x(1:end-1),'.','p')];

allconds = read_conditions(test.exp, us, test.runtype, varargin{:});

% id string for the test run
test_idstring = [...
    test.exp '_' test.runtype '_r' sprintf('%d',test.runs) '_' test.model ...
    '_' num2str(100*test.fwhm, '%.0f') 'mm' '_conds' sprintf('%d',find(ismember(allconds,test.conds)))];

% id string for localizer analysis
if ~isempty(loc)
    loc_idstring = [];
    for i = 1:length(loc)
        x = sprintf('%g-',loc(i).opval);
        loc_idstring = [...
            loc_idstring 'loc' num2str(i) '_' loc(i).exp '_' loc(i).runtype ...
            '_' loc(i).con '_' loc(i).sign  '_' loc(i).op '_' strrep(x(1:end-1),'.','p') ...
            '_r' sprintf('%d',loc(i).runs) '_' loc(i).model '_' num2str(100*loc(i).fwhm, '%.0f') 'mm']; 
    end
    if length(loc)>1
        loc_idstring = DataHash(loc_idstring);
    end
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

[blockdur, nulldur, TR, TA, stimdur, stim2scan, win] = read_scanparams(test.exp, us, test.runtype); %#ok<ASGLU>

% plateau samples
plat_win = win(win > 4 & win < blockdur + stim2scan + 1e-3);

% window
stat_idstring = [...
    'sigav_win' strrep(num2str(win(1)),'.','p') '-' strrep(num2str(win(end)),'.','p') ...
    '_plateau_win' strrep(num2str(plat_win(1)),'.','p') '-' strrep(num2str(plat_win(end)),'.','p')];

if optInputs(varargin,'nulltc')
    nullwin = win;
    stat_idstring = [stat_idstring '_nulltc'];
else
    nullwin = (TR:TR:nulldur-TR) + stim2scan;
end
    
%% Misc Setup

% string used to remember this specific analysis
main_idstring = ['psc_v2_us' num2str(us) '_' stat_idstring '_' roi_idstring '_' test_idstring '_' loc_idstring];

% if ~optInputs(varargin,'nowrite')
%     fprintf('\nS%d, %s\n',us, roi.name);
%     for i = 1:length(loc)
%         fprintf('Localizer: %s, runs %s\n',loc(i).con, sprintf('%d ',loc(i).runs))
%         fprintf('Test: run %d\n', test.run);
%     end
% end

fprintf('Localizer: %s, runs %s\n',loc(i).con, sprintf('%d ',loc(i).runs))
fprintf('Test: runs %s\n',sprintf('%d ',test.runs));

if exist([roidir main_idstring '.mat'],'file') && ~optInputs(varargin,'overwrite')
    load([roidir main_idstring '.mat']);
    return;
end

%% Read in ROI

testfeat = [...
    params('rootdir') test.exp '/analysis/fla/usub' num2str(us) '/'...
    test.runtype '_r' num2str(test.runs) '_' test.model  '_' num2str(100*test.fwhm, '%.0f') 'mm.feat/'];
    
anatfile = [...
    testfeat 'masks/' roi.name  '_func_func2highres_' func2highres  ...
    '_highres2standard_' highres2standard  '.nii.gz'];

roi_brain = readmr(anatfile,'NOPROGRESSBAR');

%% Compute signal average

funcbr = readmr([testfeat 'filtered_func_data.nii.gz'], 'NOPROGRESSBAR');

dims = size(funcbr.data);
func_matrix = reshape(shiftdim(funcbr.data, 3), dims(4), prod(dims(1:3)));
xi = logical(brainops(roi_brain.data, roi.op, roi.opval));
func_matrix_roi = func_matrix(:,xi(:));

% timings data
behavdat = read_timings(test.exp, us, test.runtype, test.runs, test.model, varargin{:});

% BOLD signal at null time points
onsets = behavdat.onsets(strcmp('NULL',behavdat.conds));
null_tps = nan( length(nullwin), size(func_matrix_roi,2), length(onsets));
for j = 1:length(onsets)
    null_tps(:,:,j) = func_matrix_roi(round((onsets(j) + nullwin)/TR) + 1,:);
end

if optInputs(varargin, 'nulltc')
    nullmean = squeeze(mean(null_tps, 3));
else
    nullmean = ones(length(win),1)*squeeze(mean(mean(null_tps, 3),1));
end

nblocks = length(behavdat.onsets(strcmp(test.conds{1},behavdat.conds)));
psc = nan(length(win), size(func_matrix_roi,2), nblocks, length(test.conds));
for i = 1:length(test.conds)
    onsets = behavdat.onsets(strcmp(test.conds{i},behavdat.conds));
    for j = 1:length(onsets)
        x = func_matrix_roi(round((onsets(j) + nullwin)/TR) + 1,:);
        psc(:,:,j,i) = 100*(x-nullmean)./ nullmean;
    end
end

%% Localizer difference matrix

groups = regexp(loc.con,'_vs_','split');
psc_groups = nan(length(win), size(func_matrix_roi,2), nblocks, length(groups));
for j = 1:length(groups)
    evn = evgroup2evname(test.exp,us,groups{j},test.model);
    inds = nan(1,length(evn));
    for k = 1:length(evn)
        inds(k) = find(strcmp(evn{k},test.conds));
    end
    psc_groups(:,:,:,j) = mean(psc(:,:,:,inds),4);
end

if length(groups) == 1
    psc_diff = psc_groups;
elseif length(groups) == 2;
    psc_diff = psc_groups(:,:,:,1) - psc_groups(:,:,:,2);
else
    error('Contrast should be a difference between two condition groups');
end

%% Leave-one-out analysis

voxrange = round((size(func_matrix_roi,2)-1) * loc.opval) + 1;
psc_diff_plat_win = squeeze(mean(psc_diff(ismember(win, plat_win),:,:),1))';
psc_func_roi = nan(nblocks, length(win), length(test.conds));
for i = 1:nblocks
    locblocks = setdiff(1:nblocks,i);
    [~,p,~,t] = ttest(psc_diff_plat_win(locblocks, :));
    logP = -10*log10(p) .* sign(t.tstat);
    [~,xi] = sort(t.tstat,'descend');
    psc_func_roi(i,:,:) = squeeze(mean(psc(:,xi(voxrange(1):voxrange(2)),i,:),2)); 
end

psc_func_roi_mean = squeeze(mean(psc_func_roi,1));

%%

save([roidir main_idstring '.mat'],'psc_func_roi', 'psc_func_roi_mean');

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





