function naturalsound_linear_analysis_surf_v2

% subjects
varargin = {};
usubs = [1 11 106 121 122 123 143 160 161 164]; % excluded 156 (abs motion), excluded 37 (behavior), 159 (behavior)

% preprocessing parameters
fwhm = 3;
normalization_type = 'sumnorm_demeanSubjects';
roi = 'hand-audctx';

% directory
soundZ_thresh = 0;
soundP_thresh = 3; % p < 0.001
corrP_thresh = 0;
regressP_thresh = 0;
regressExVar_thresh = 30;
% percPos_thresh = 0;
analysis_name = ['us' sprintf('-%d', usubs) '_sZ-' num2str(soundZ_thresh) '_sP-' num2str(soundP_thresh) '_cP-' strrep(num2str(corrP_thresh),'.','p') '_exV-' strrep(num2str(regressExVar_thresh),'.','p') '_fwhm-' num2str(fwhm) '_' roi '_' normalization_type];
figure_directory = [params('rootdir') 'naturalsound/figures/linear_surf_analysis_downsamp_v2/' analysis_name '/'];
if ~exist(figure_directory, 'dir')
    mkdir(figure_directory);
end
addpath(genpath('FastICA_25'));
addpath(genpath('functional_systems'));
addpath(genpath('export_fig'));
freesurfer_version = read_freesurfer_version('naturalsound');

%% Read in data

maxruns = 3;
data_directory = [figure_directory 'data/'];
if ~exist(data_directory,'dir')
    mkdir(data_directory);
end
matfile = [data_directory 'voxel_matrix_unfilt.mat'];

if ~exist(matfile,'file')
    
    % initialize
    hemis = {'rh','lh'};
    patch = cell(1,2); vras = cell(1,2); vi = cell(1,2); roi_mask = cell(1,2);
    grid_x = cell(1,2); grid_y = cell(1,2); bounds = cell(1,2); vgrid = cell(1,2);
    
    for i = 1:2
        
        % patch coordinates and surface voxel indices
        if strcmp(hemis{i},'lh')
            patch{i} = read_patch(['~/freesurfer/fsaverage/surf/' hemis{i} '.cortex2.patch.flat']);
        elseif strcmp(hemis{i}, 'rh')
            patch{i} = read_patch(['~/freesurfer/fsaverage/surf/' hemis{i} '.cortex.patch.flat']);
        end
        
        % roi_mask
        roi_mask{i} = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
        
        % vertices
        % sign of patch vertices indicates whether face is on border of
        % patch, which is not relevant for selecting the vertices
        [~,ai,ib] = intersect(abs(patch{i}.vnums),roi_mask{i}.vnums);
        vras{i} = patch{i}.vras(ai,1:2); % 2-D coordinates
        vi{i} = roi_mask{i}.vnums(ib); % indices within roi
        
        % bounds for regridding, 2 mm resolution
        bounds{i} = [floor(min(vras{i})); ceil(max(vras{i}))]; % min and max of 2-d coordinates
        [grid_x{i},grid_y{i}] = meshgrid(bounds{i}(1,1):2:bounds{i}(2,1),(bounds{i}(1,2):2:bounds{i}(2,2))); % regrid in units of mm
        
        % initialize voxel grid
        vgrid{i} = nan(size(grid_x{i},1), size(grid_y{i},2), length(usubs), 165, maxruns);
        
    end
    
    % subjects
    for j = 1:length(usubs)
        
        % subjid
        subjid = ['naturalsound_us' num2str(usubs(j))];
        
        % number of runs for current subject
        runs = read_runs('naturalsound',usubs(j),'main_v3_combined');
        nruns = length(runs);
        
        % voxel grid
        for k = 1:nruns
            
            % hemispheres
            for i = 1:2
                
                % grid each sound
                subject_grid_file = [data_directory 'vgrid_us' num2str(usubs(j)) '_r' num2str(k) '_' hemis{i} '.mat'];
                if ~exist(subject_grid_file,'file')
                    
                    fprintf('\nhemi %s, usub %d, scan %d\n', hemis{i}, usubs(j), runs(k)); drawnow;
                    
                    % initialize grid
                    vgrid_subject = nan(size(grid_x{i},1),size(grid_y{i},2),165);
                    
                    % read in signal-averaged data
                    %                     fname = [params('rootdir') 'freesurfer/' subjid '/fla/main_combined_r' num2str(runs(k)) '/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp5-11.mgz'];
                    fname = [params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r' num2str(runs(k)) '/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp3-6.mgz'];
                    sigav_vol = MRIread(fname);
                    vr = squeeze(sigav_vol.vol(:,vi{i}+1,:,:))';
                    
                    % add boundary points as NaN values, useful for interpolation
                    vras_boundary_pts = [vras{i}; bounds{i}(1,:); bounds{i}(2,:); bounds{i}(1,1), bounds{i}(2,2); bounds{i}(2,1), bounds{i}(1,2)];
                    vr_boundary_pts = [vr,nan(165,4)];
                    
                    % grided nan mask
                    mask = all(~isnan(vr_boundary_pts));
                    grid_mask = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),double(mask),grid_x{i},grid_y{i},'natural');
                    
                    % regridding
                    for q = 1:165
                        fprintf('%d ', q); drawnow;
                        grid_response = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),vr_boundary_pts(q,:),grid_x{i},grid_y{i},'natural');
                        grid_response(grid_mask<0.99) = NaN;
                        vgrid_subject(:,:,q) = grid_response;
                        %                     figure;
                        %                     surf(grid_x{j}{i}, grid_y{j}{i}, 100*grid_response,'EdgeColor','none');
                        %                     caxis([0 3]);
                        %                     view(0, 90);
                    end
                    save(subject_grid_file, 'vgrid_subject');
                else
                    load(subject_grid_file);
                end
                
                % save grid
                vgrid{i}(:,:,j,:,k) = vgrid_subject;
                
                %                 vgrid{i}(:,:,j,q,k) = griddata(vras{i}(:,1),vras{i}(:,2),vr(q,:),grid_x{i},grid_y{i});
                %                 vgrid{i}(:,:,j,q,k) = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),vr_boundary_pts(q,:),grid_x{i},grid_y{i});
                
                %                 g = griddata(vras_patch{i}(:,1),vras_patch{i}(:,2),vr_patch(1,:),grid_x{i},grid_y{i});
                %                 g = griddata(vras{i}(:,1),vras{i}(:,2),vr(1,:),grid_x{i},grid_y{i});
                
                %                 figure;
                %                 q = 1;
                %                 scatter(vras_boundary_pts(:,1),vras_boundary_pts(:,2),3,100*vr_boundary_pts(q,:));
                %                 caxis([0 3]);
                %
                %                 figure;
                %                 q = 1;
                %                 scatter(vras_boundary_pts(:,1),vras_boundary_pts(:,2),3,double(mask));
                %                 caxis([0 1]);
                %
                %                 figure;
                %                 gx_mask = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),double(mask),grid_x{i},grid_y{i},'natural');
                %                 surf(grid_x{1}, grid_y{1}, gx_mask,'EdgeColor','none');
                %                 caxis([0 1]);
                %                 view(0, 90);
                %
                %                 figure;
                %                 gx = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),vr_boundary_pts(q,:),grid_x{i},grid_y{i},'natural');
                %                 surf(grid_x{1}, grid_y{1}, 100*gx.*double(gx_mask>0.99),'EdgeColor','none');
                %                 view(0, 90);
                %                 caxis([0 3]);
                
                %                 figure;
                %                 q = 1;
                %                 scatter(vras{i}(:,1),vras{i}(:,2),3,100*vr(q,:));
                %                 caxis([0 3]);
                %
                %                 figure;
                %                 q = 1;
                %                 scatter(vras{i}(:,1),vras{i}(:,2),3,double(~isnan(vr(q,:))));
                %                 caxis([0 3]);
                %
                %                 figure;
                %                 mask = all(~isnan(vr));
                %                 gx_mask = griddata(vras{i}(:,1),vras{i}(:,2),double(mask),grid_x{i},grid_y{i},'natural');
                %                 surf(grid_x{1}, grid_y{1}, gx_mask,'EdgeColor','none');
                %                 caxis([0 1])
                %                 view(0, 90);
                %                 figure;
                %                 gx = griddata(vras{i}(:,1),vras{i}(:,2),vr(q,:),grid_x{i},grid_y{i},'natural');
                %                 surf(grid_x{1}, grid_y{1}, 100*gx.*double(gx_mask>0.99),'EdgeColor','none');
                %                 view(0, 90);
                %                 caxis([0 3]);
                
                %                 figure;
                %                 scatter(vras{i}(:,1),vras{i}(:,2),3,vr(1,:));
                %                 figure;
                %                 scatter(vras_patch{i}(:,1),vras_patch{i}(:,2),3,vr_patch(1,:));
                %                 figure;
                %                 surf(grid_x{i},grid_y{i},g);
                %                 view(0, 90);
                %                 figure;
                %                 mesh(grid_x{i},grid_y{i},m);
                %                 view(0, 90);
                
            end
        end
    end
    save(matfile, 'vras','vi','mask','grid_x','grid_y','vgrid','-v7.3');
else
    load(matfile);
end

% roll-out the grid of voxels for each sound/subject/scan to a vector
n_voxels_rh = size(vgrid{1},1)*size(vgrid{1},2);
vgrid_reshape_rh = reshape(vgrid{1}, [n_voxels_rh, length(usubs), 165, maxruns]);
n_voxels_lh = size(vgrid{2},1)*size(vgrid{2},2);
vgrid_reshape_lh = reshape(vgrid{2}, [n_voxels_lh, length(usubs), 165, maxruns]);

% combine right and left hemisphere voxels
n_voxels = n_voxels_lh + n_voxels_rh;
v_hemi_combined = cat(1, vgrid_reshape_rh, vgrid_reshape_lh);

% combine across subjects
% each row: [rh_vox_s1, lh_vox_s1, rh_vox_s2, lh_vox_s2, ...]
v_unfilt_allscans = nan([165, n_voxels*length(usubs), maxruns]);
for i = 1:maxruns
    v_unfilt_allscans(:,:,i) = reshape(v_hemi_combined(:,:,:,i), [n_voxels*length(usubs), 165])';
end
v_unfilt_mean = nanmean(v_unfilt_allscans,3);

% v_unfilt_s1 = reshape(v_hemi_combined(:,:,:,1), [n_voxels*length(usubs), 165])'; % stimuli x voxels
% v_unfilt_s2 = reshape(v_hemi_combined(:,:,:,2), [n_voxels*length(usubs), 165])'; % stimuli x voxels
% v_unfilt_s3 = reshape(v_hemi_combined(:,:,:,3), [n_voxels*length(usubs), 165])'; % stimuli x voxels
% v_unfilt_s12 = (v_unfilt_s1+v_unfilt_s2)/2;
% v_unfilt_s123 = (v_unfilt_s1+v_unfilt_s2+v_unfilt_s3)/3;

% subject indices
si_unfilt = ones(n_voxels,1)*usubs;
si_unfilt = si_unfilt(:);

% grid indices
gridi_unfilt = repmat([(1:n_voxels_rh)';(1:n_voxels_lh)'],length(usubs),1);

% hemisphere indices, 1 = left, 2 = right
hi_unfilt = repmat([2*ones(n_voxels_rh,1); ones(n_voxels_lh,1)],length(usubs),1);

%% Stats

stats_file = [data_directory 'stats.mat'];
if ~exist(stats_file,'file') || optInputs(varargin, 'overwrite')
    
    % p-value sound threshold
    [~,x] = ttest(v_unfilt_mean, 0, 0.05, 'right');
    soundP = -log10(x);
    soundP_individual_scans = nan(length(soundP), maxruns);
    for i = 1:maxruns
        [~,x] = ttest(v_unfilt_allscans(:,:,i), 0, 0.05, 'right');
        soundP_individual_scans(:,i) = -log10(x);
    end
    
    % sound-responsivitiy
    soundZ = mean(v_unfilt_mean) ./ std(v_unfilt_mean);
    soundZ_individual_scans = nan(length(soundZ), maxruns);
    for i = 1:maxruns
        soundZ_individual_scans(:,i) = mean(v_unfilt_allscans(:,:,i)) ./ std(v_unfilt_allscans(:,:,i));
    end
    
    % percent of positive sounds
    percPos = 100*mean(v_unfilt_mean>0);
    
    % cross-scan correlation
    corrR = fastcorr3(v_unfilt_allscans(:,:,1), v_unfilt_allscans(:,:,2));
    corrT = corrR .* sqrt((165-2) ./ (1 - corrR.^2));
    corrP = -log10(tpvalue_copy(-abs(corrT), 165-2));
    
    % cross-scan regression
    regressP = nan(size(corrP));
    regressSE = nan(size(corrP));
    regressB = nan(size(corrP));
    regressE = nan(165,length(corrP));
    regressExVar = nan(size(corrP));
    corr_norm_reduction = nan(size(corrP));
    for i = 1:size(v_unfilt_allscans,2)
        %     if mod(i,1000)==0
        %         fprintf('%d\n',i);
        %     end
        x = v_unfilt_allscans(:,i,1);
        y = v_unfilt_allscans(:,i,2);
        [regressB(i),regressSE(i),regressP(i),regressE(:,i)] = regress_SNH(x,y);
        regressExVar(i) = 100*(sum(y.^2) - sum(regressE(:,i).^2)) / sum(y.^2);
        corr_norm_reduction(i) = norm_reduction(x,y);
    end
    save(stats_file, 'soundP','soundP_individual_scans','corrR','corrP','soundZ','soundZ_individual_scans','percPos','regressP','regressExVar','corr_norm_reduction');
else
    load(stats_file);
end

%% Voxel Selection

% regressP histogram
[y,x] = hist(regressP(~isnan(regressP)), 100);
figure;
set(gcf, 'Position', [0 0 1200 400]);
subplot(1,2,1);
plot(x,y,'LineWidth', 2);
hold on;
h = plot(regressP_thresh*[1 1], [0 max(y)], '--k', 'LineWidth', 2);
ylim([0 max(y)]);
xlabel('-log10(P)');
ylabel('Number of Voxels');
subplot(1,2,2);
plot(x,sum(y)-cumsum(y),'LineWidth', 2);
hold on;
h = plot(regressP_thresh*[1 1], [0 sum(y)], '--k', 'LineWidth', 2);
ylim([0 sum(y)]);
xlabel('-log10(P)');
ylabel('1-cumsum(NumVox)');
export_fig([figure_directory 'regressP_histogram.pdf'],'-pdf','-nocrop');

% regress explained varaince histogram
[y,x] = hist(regressExVar(~isnan(regressExVar)), 100);
figure;
set(gcf, 'Position', [0 0 1200 400]);
subplot(1,2,1);
plot(x,y,'LineWidth', 2);
hold on;
h = plot(regressExVar_thresh*[1 1], [0 max(y)], '--k', 'LineWidth', 2);
ylim([0 max(y)]);
xlabel('% Explained Variance');
ylabel('Number of Voxels');
subplot(1,2,2);
plot(x,sum(y)-cumsum(y),'LineWidth', 2);
hold on;
h = plot(regressExVar_thresh*[1 1], [0 sum(y)], '--k', 'LineWidth', 2);
ylim([0 sum(y)]);
xlabel('% Explained Variance');
ylabel('1-cumsum(NumVox)');
export_fig([figure_directory 'regressExVar_histogram.pdf'],'-pdf','-nocrop');

% corrR histogram
[y,x] = hist(corrR(~isnan(corrR)), 100);
figure;
set(gcf, 'Position', [0 0 1200 400]);
subplot(1,2,1);
plot(x,y,'LineWidth', 2);
hold on;
h = plot([0 0], [0 max(y)], '--k', 'LineWidth', 2);
ylim([0 max(y)]);
xlabel('Correlation (r)');
ylabel('Number of Voxels');
subplot(1,2,2);
plot(x,sum(y)-cumsum(y),'LineWidth', 2);
hold on;
h = plot([0 0], [0 sum(y)], '--k', 'LineWidth', 2);
ylim([0 sum(y)]);
xlabel('Correlation (r)');
ylabel('1-cumsum(NumVox)');
export_fig([figure_directory 'corrR_histogram.pdf'],'-pdf','-nocrop');

% corrP histogram
[y,x] = hist(corrP(~isnan(corrP)), 100);
figure;
set(gcf, 'Position', [0 0 1200 400]);
subplot(1,2,1);
plot(x,y,'LineWidth', 2);
hold on;
h = plot(corrP_thresh*[1 1], [0 max(y)], '--k', 'LineWidth', 2);
ylim([0 max(y)]);
% legend(h, ['p < ' num2str(pthresh)]);
xlabel('-log10(p)');
ylabel('Number of Voxels');
legend(h, ['-log10(p) < ' num2str(corrP_thresh)]);
subplot(1,2,2);
plot(x,sum(y)-cumsum(y),'LineWidth', 2);
hold on;
h = plot(corrP_thresh*[1 1], [0 sum(y)], '--k', 'LineWidth', 2);
ylim([0 sum(y)]);
xlabel('-log10(p)');
ylabel('1-cumsum(NumVox)');
legend(h, ['-log10(p) < ' num2str(corrP_thresh)]);
export_fig([figure_directory 'corrP_histogram.pdf'],'-pdf','-nocrop');

% soundZ histogram
[y,x] = hist(min(soundZ_individual_scans,[],2), 100);
figure;
set(gcf, 'Position', [0 0 1200 400]);
subplot(1,2,1);
plot(x,y,'LineWidth', 2);
hold on;
h = plot(soundZ_thresh*[1 1], [0 max(y)], '--r', 'LineWidth', 2);
ylim([0 max(y)]);
legend(h, ['Z > ' num2str(soundZ_thresh)]);
xlabel('Z values');
ylabel('Number of Voxels');
subplot(1,2,2);
plot(x,sum(y)-cumsum(y),'LineWidth', 2);
hold on;
h = plot(soundZ_thresh*[1 1], [0 sum(y)], '--r', 'LineWidth', 2);
legend(h, ['Z > ' num2str(soundZ_thresh)]);
ylim([0 sum(y)]);
xlabel('Z values');
ylabel('1-cumsum(NumVox)');
export_fig([figure_directory 'soundZ_histogram.pdf'],'-pdf','-nocrop');

% soundP histogram
[y,x] = hist(min(soundP_individual_scans,[],2), 100);
figure;
set(gcf, 'Position', [0 0 1200 400]);
subplot(1,2,1);
plot(x,y,'LineWidth', 2);
hold on;
h = plot(soundP_thresh*[1 1], [0 max(y)], '--r', 'LineWidth', 2);
ylim([0 max(y)]);
legend(h, ['-log10(p) < ' num2str(soundP_thresh)]);
xlabel('-log10(p)');
ylabel('Number of Voxels');
subplot(1,2,2);
plot(x,sum(y)-cumsum(y),'LineWidth', 2);
hold on;
h = plot(soundP_thresh*[1 1], [0 sum(y)], '--r', 'LineWidth', 2);
legend(h, ['-log10(p) < ' num2str(soundP_thresh)]);
ylim([0 sum(y)]);
xlabel('-log10(p)');
ylabel('1-cumsum(NumVox)');
export_fig([figure_directory 'soundP_histogram.pdf'],'-pdf','-nocrop');

% percPos histogram
% [y,x] = hist(percPos, 100);
% figure;
% set(gcf, 'Position', [0 0 1200 400]);
% subplot(1,2,1);
% plot(x,y,'LineWidth', 2);
% hold on;
% h = plot(percPos_thresh*[1 1], [0 max(y)], '--r', 'LineWidth', 2);
% ylim([0 max(y)]);
% legend(h, ['percPos > ' num2str(percPos_thresh)]);
% xlabel('-log10(p)');
% ylabel('Number of Voxels');
% subplot(1,2,2);
% plot(x,sum(y)-cumsum(y),'LineWidth', 2);
% hold on;
% h = plot(percPos_thresh*[1 1], [0 sum(y)], '--r', 'LineWidth', 2);
% legend(h, ['percPos > ' num2str(percPos_thresh)]);
% ylim([0 sum(y)]);
% xlabel('-log10(p)');
% ylabel('1-cumsum(NumVox)');
% export_fig([figure_directory 'percPos_histogram.pdf']   ,'-pdf','-nocrop');
% fprintf('Mean soundP for fraction index is %.1f\n',mean(soundP(percPos>percPos_thresh)));


% select voxels
% ResetRandStream(0);
% voxel_selection1 = percPos > percPos_thresh & soundP > soundP_thresh & soundZ > soundZ_thresh & corrP > corrP_thresh & regressP > regressP_thresh;

usubs_two_scans = [122 161];
usubs_three_scans = [1 106 11 121 123 143 160 164];

voxel_selection = false(size(soundP));
xi = ismember(si_unfilt, usubs_two_scans);
voxel_selection(xi) = all(soundP_individual_scans(xi,1:2)' > soundP_thresh) & regressExVar(xi) > regressExVar_thresh;
% voxel_selection(xi) = percPos(xi) > percPos_thresh & all(soundP_individual_scans(xi,1:2)' > soundP_thresh) & all(soundZ_individual_scans(xi,1:2)' > soundZ_thresh) & corrP(xi) > corrP_thresh & regressP(xi) > regressP_thresh & regressExVar(xi) > regressExVar_thresh;
xi = ismember(si_unfilt, usubs_three_scans);
voxel_selection(xi) = all(soundP_individual_scans(xi,1:3)' > soundP_thresh) & regressExVar(xi) > regressExVar_thresh;
% voxel_selection(xi) = percPos(xi) > percPos_thresh & all(soundP_individual_scans(xi,1:3)' > soundP_thresh) & all(soundZ_individual_scans(xi,1:3)' > soundZ_thresh) & corrP(xi) > corrP_thresh & regressP(xi) > regressP_thresh & regressExVar(xi) > regressExVar_thresh;

% voxel_selection_justscan1 = soundP_scan1 > soundP_thresh & soundZ_scan1 > soundZ_thresh;
% voxel_selection_justscan2 = soundP_scan2 > soundP_thresh & soundZ_scan2 > soundZ_thresh;

fprintf('%d total voxels (%.1f%%)\n', sum(voxel_selection), 100*sum(voxel_selection)/sum(all(~isnan(v_unfilt_allscans(:,:,1)))));

% soundZ_selected = soundZ(xi);
% vox_response_matrix_scan1 = vox_response_matrix_scan1_raw(:,xi);
% vox_response_matrix_scan2 = vox_response_matrix_scan2_raw(:,xi);
% fprintf('Selecting %.1f%% of voxels\n',100*mean(xi));
% if cluster_scan1
%     vox_response_matrix = vox_response_matrix_scan1;
% else
%     vox_response_matrix = vox_response_matrix_scan1 + vox_response_matrix_scan2;
% end

si = si_unfilt(voxel_selection);
hi = hi_unfilt(voxel_selection);
% vox_response_matrix_scan1_raw(vox_response_matrix_scan1_raw < 0) = 0;
% vox_response_matrix_scan2_raw(vox_response_matrix_scan2_raw < 0) = 0;

switch normalization_type
    
    %     case 'sumnorm_demean'
    %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
    %         [v, grandmean] = sumnorm_demean_subjects(v_unfilt_s12(:,voxel_selection),ones(size(si)));
    %         v1 = sumnorm_demean_subjects(v_unfilt_s1(:,voxel_selection),ones(size(si)));
    %         v2 = sumnorm_demean_subjects(v_unfilt_s2(:,voxel_selection),ones(size(si)));
    %         v1_selected_just_scan1 = sumnorm_demean_subjects(v_unfilt_s1(:,voxel_selection_justscan1),ones(size(si)));
    %         v2_selected_just_scan2 = sumnorm_demean_subjects(v_unfilt_s2(:,voxel_selection_justscan2),ones(size(si)));
    %
    case 'demeanSubjects'
        %             v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %             [v, normed_data] = demean_subjects(v_unfilt_s12(:,voxel_selection),si);
        %             v1 = demean_subjects(v_unfilt_s1(:,voxel_selection),si);
        %             v2 = demean_subjects(v_unfilt_s2(:,voxel_selection),si);
        %             grandmean = mean(v,2);
        
        [v, grandmean] = demean_subjects(v_unfilt_mean(:,voxel_selection),si);
        v1 = demean_subjects(v_unfilt_allscans(:,voxel_selection,1),si);
        v2 = demean_subjects(v_unfilt_allscans(:,voxel_selection,2),si);
        v3 = demean_subjects(v_unfilt_allscans(:,voxel_selection,3),si);
        v12 = demean_subjects(mean(v_unfilt_allscans(:,voxel_selection,1:2),3),si);
        % %         naturalsound_acoustic_category_figures(grandmean, figure_directory, 'grandmean-norm-');
        %
        %     case 'demeanSubjects_unitvar'
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         [v, normed_data] = demean_subjects_unitvar(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = demean_subjects_unitvar(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = demean_subjects_unitvar(v_unfilt_s2(:,voxel_selection),si);
        %         grandmean = mean(normed_data,2);
        %         %         naturalsound_acoustic_category_figures(grandmean, figure_directory, 'grandmean-norm-');
        
    case 'sumnorm_demeanSubjects'
        v_unfilt_s12 = mean(v_unfilt_allscans(:,:,1:2),3);
        [v, grandmean] = sumnorm_demean_subjects(v_unfilt_mean(:,voxel_selection),si);
        [v1, v1_grandmean] = sumnorm_demean_subjects(v_unfilt_allscans(:,voxel_selection,1),si);
        [v2, v2_grandmean] = sumnorm_demean_subjects(v_unfilt_allscans(:,voxel_selection,2),si);
        v3 = sumnorm_demean_subjects(v_unfilt_allscans(:,voxel_selection,3),si);
        v12 = sumnorm_demean_subjects(mean(v_unfilt_allscans(:,voxel_selection,1:2),3),si);
        
        %     case 'sumnorm_demeanSubjects_addgrandmean'
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         [v, grandmean] = sumnorm_demean_subjects_addgrandmean(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = sumnorm_demean_subjects_addgrandmean(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = sumnorm_demean_subjects_addgrandmean(v_unfilt_s2(:,voxel_selection),si);
        %         v1_selected_just_scan1 = sumnorm_demean_subjects_addgrandmean(v_unfilt_s1(:,voxel_selection_justscan1),si);
        %         v2_selected_just_scan2 = sumnorm_demean_subjects_addgrandmean(v_unfilt_s2(:,voxel_selection_justscan2),si);
        %
        %     case 'sumnorm_unitvar'
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         [v, grandmean] = sumnorm_demean_unitvar_subjects(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = sumnorm_demean_unitvar_subjects(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = sumnorm_demean_unitvar_subjects(v_unfilt_s2(:,voxel_selection),si);
        %     case 'sumnorm_unitcov_within_subjects'
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         [v, grandmean] = sumnorm_demean_unitcov_within_subjects(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = sumnorm_demean_unitcov_within_subjects(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = sumnorm_demean_unitcov_within_subjects(v_unfilt_s2(:,voxel_selection),si);
        %     case 'L1norm'
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         v = L1norm_demean_subjects(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = L1norm_demean_subjects(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = L1norm_demean_subjects(v_unfilt_s2(:,voxel_selection),si);
        %     case 'L1norm_unitvar'
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         v = L1norm_demean_unitvar_subjects(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = L1norm_demean_unitvar_subjects(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = L1norm_demean_unitvar_subjects(v_unfilt_s2(:,voxel_selection),si);
        %     case 'truncate_L1norm'
        %         v_unfilt_s1(v_unfilt_s1<0) = 0;
        %         v_unfilt_s2(v_unfilt_s2<0) = 0;
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         v = L1norm_demean_subjects(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = L1norm_demean_subjects(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = L1norm_demean_subjects(v_unfilt_s2(:,voxel_selection),si);
        %     case 'L2norm'
        %         v_unfilt_s12 = (v_unfilt_s1 + v_unfilt_s2)/2;
        %         [v, grandmean] = L2norm_demean_subjects(v_unfilt_s12(:,voxel_selection),si);
        %         v1 = L2norm_demean_subjects(v_unfilt_s1(:,voxel_selection),si);
        %         v2 = L2norm_demean_subjects(v_unfilt_s2(:,voxel_selection),si);
    otherwise
        error('Normalization type not found.');
end

% v = L1norm_demean_subjects(v_unfilt_s12(:,voxel_selection),si);
% v1 = L1norm_demean_subjects(v_unfilt_s1(:,voxel_selection),si);
% v2 = L1norm_demean_subjects(v_unfilt_s2(:,voxel_selection),si);

% v12 = L1norm_demean(v_unfilt_s12(:,voxel_selection),si);
% x1 = v_unfilt_s1(:,voxel_selection);
% x2 = v_unfilt_s2(:,voxel_selection);
% xi = Shuffle((1:165)'*ones(1,size(x1,2)) + ones(165,1)*(0:size(x1,2)-1)*165);
% % v_shuffle = L1norm_demean_subjects((x1(xi) + x2)/2,si);

% n_stims = size(v,1);
% n_voxels = size(v,2);

save([figure_directory 'vox_matrix.mat'], 'v', 'v1', 'v2', 'v3', 'v12', 'si');
close all;

v1_threescans = v1(:,ismember(si,usubs_three_scans));
v2_threescans = v2(:,ismember(si,usubs_three_scans));
v12_threescans = v12(:,ismember(si,usubs_three_scans));
v3_threescans = v3(:,ismember(si,usubs_three_scans));

%% Voxel matrix

ResetRandStream(1);
figure;
imagesc(v(Shuffle(1:165),  Shuffle(1:size(v,2))),std(v(:))*[-1.5 1.5]);
print('-dpng',[figure_directory 'voxel_matrix.png'],'-r200');

%% 
v_mode_centered = nan(size(v));
v12_mode_centered = nan(size(v));
v3_mode_centered = nan(size(v));
for i = 1:165
    bins = linspace(-5,5,100);
    x = v(i,:)';
    x = x/std(x);
    yh = hist(x,bins);
    % figure;
    % plot(bins,yh/sum(yh), 'r'); hold on;
    h = normpdf(bins,0,0.5);
    yh = conv(yh,h/sum(h),'same');
    % plot(bins,yh/sum(yh), 'g');
    [~,xi] = max(yh);
    x = x - bins(xi);
    % yh = hist(x,bins);
    % plot(bins,yh/sum(yh), 'b');
    % yL = ylim;
    % plot([0 0], yL,'k--');
    v_mode_centered(i,:) = x;
    
    bins = linspace(-5,5,100);
    x = v12(i,:)';
    x = x/std(x);
    yh = hist(x,bins);
    % figure;
    % plot(bins,yh/sum(yh), 'r'); hold on;
    h = normpdf(bins,0,0.5);
    yh = conv(yh,h/sum(h),'same');
    % plot(bins,yh/sum(yh), 'g');
    [~,xi] = max(yh);
    x = x - bins(xi);
    % yh = hist(x,bins);
    % plot(bins,yh/sum(yh), 'b');
    % yL = ylim;
    % plot([0 0], yL,'k--');
    v12_mode_centered(i,:) = x;
    
    bins = linspace(-5,5,100);
    x = v3(i,:)';
    nonan = ~isnan(x);
    x = x(nonan);
    x = x/std(x);
    yh = hist(x,bins);
    % figure;
    % plot(bins,yh/sum(yh), 'r'); hold on;
    h = normpdf(bins,0,0.5);
    yh = conv(yh,h/sum(h),'same');
    % plot(bins,yh/sum(yh), 'g');
    [~,xi] = max(yh);
    x = x - bins(xi);
    % yh = hist(x,bins);
    % plot(bins,yh/sum(yh), 'b');
    % yL = ylim;
    % plot([0 0], yL,'k--');
    v3_mode_centered(i,nonan) = x;
    
end

%% PCA Analyses

% cross-validation analysis
pca_reliability(v, v12, v3, usubs_three_scans, si, [figure_directory 'crossval-PCA-']);
pca_reliability_mediancorr(v, v12, v3, usubs_three_scans, si, [figure_directory 'crossval-PCA-medcorr-']);
pca_reliability_mediancorr(v1, v12, v3, usubs_three_scans, si, [figure_directory 'crossval-PCA-medcorr-scan1-']);
pca_reliability_mediancorr(v12, v12, v3, usubs_three_scans, si, [figure_directory 'crossval-PCA-medcorr-scan12-']);

% subject specificity of pca components
[U,S,V] = svd(v,'econ');
subject_specific_variance = nan(length(usubs),165);
for i = 1:length(usubs)
    subject_specific_variance(i,:) = mean(V(si == usubs(i),:).^2,1);
end
subject_specific_variance = subject_specific_variance ./ (ones(length(usubs),1)*sum(subject_specific_variance));
figure;
imagesc(subject_specific_variance(:,1:20),[0 0.3]);
xlabel('PCs'); ylabel('Subjects'); colorbar;
export_fig([figure_directory 'PCA_variance_across_subjects.png'],'-png','-nocrop','-r100');

% skewness of pca compoennts
skew = sum(subject_specific_variance.^3);
figure;
plot(skew(1:20),'-o')
ylim([0 max(skew(1:20))*1.2]); xlabel('Skew');
box off;
export_fig([figure_directory 'PCA_subject_skew.pdf'],'-pdf','-transparent');

% histogram
plothist_bysubject(V(:,1:10)'*sqrt(size(V,1)-1), si, [figure_directory 'PCA-10-hist-bysubject']);

% explained variance
figure;
set(gca,'FontSize',20);
exvar = diag(S).^2/sum(diag(S).^2);
bar(1:21,100*exvar([1:20,165]),'FaceColor',[1 1 1]);
ylim([0 max(exvar)*1.2]*100); xlim([0 22]);
set(gca,'XTick',1:2:21,'XTickLabel',[1:2:20,165]);
xlabel('Principle Components'); ylabel('% of Explained Variance');
box off;
export_fig([figure_directory 'PCA_exvar.pdf'],'-pdf','-nocrop','-transparent');

% cumulative explained variance
figure;
set(gca,'FontSize',20);
exvar = diag(S).^2/sum(diag(S).^2);
bar(1:21,100*[cumsum(exvar(1:20))',1],'FaceColor',[1 1 1]);
ylim([0 1]*100); xlim([0 22]);
set(gca,'XTick',[1:2:20],'XTickLabel',[1:2:20,165],'YTick',[20 40 60]);
xlabel('Number of PCs'); ylabel('Cumulative Explained Variance');
box off;
export_fig([figure_directory 'PCA_cumexvar.pdf'],'-pdf','-nocrop','-transparent');

Vnorm = nan(size(V,1),165);
for i = 1:165
    Vnorm(:,i) = V(:,i)/std(V(:,i));
end

% scatter plot of pca components
myscatter2(Vnorm(:,1:5),si,'allpairs',3,[figure_directory 'PCA_scatter']);
mycontour(Vnorm(:,1:5),'allpairs',[figure_directory 'PCA_contour']);

% kurtosis
figure;
bar(1:20,kurtosis(Vnorm(:,1:20))-3,'FaceColor',[1 1 1]);
xlabel('PCA Components'); ylabel('Kurtosis-3');
xlim([0 21]);
box off;
export_fig([figure_directory 'PCA_kurtosis.pdf'],'-pdf','-nocrop','-transparent');

% skewness
figure;
bar(1:20,abs(skewness(Vnorm(:,1:20))),'FaceColor',[1 1 1]);
xlabel('PCA Components'); ylabel('abs(skewness)');
xlim([0 21]);
export_fig([figure_directory 'PCA_skewness.pdf'],'-pdf','-nocrop','-transparent');

% entropy
x = nan(1,20);
for i = 1:20
    x(i) = log(sqrt(2*pi*exp(1)))-entropy(Vnorm(:,i)');
end
figure;
bar(1:20,x,'FaceColor',[1 1 1]);
xlabel('PCA Components'); ylabel('Negentropy (Nats)');
xlim([0 21]);
export_fig([figure_directory 'PCA_entropy.pdf'],'-pdf','-nocrop','-transparent');

% acoustic and category figures
Ur = U(:,1:20);
UrNorm = nan(size(Ur));
for i = 1:20
    UrNorm(:,i) = (Ur(:,i) - min(Ur(:,i)))/(max(Ur(:,i))-min(Ur(:,i)));
end
naturalsound_acoustic_category_figures(UrNorm, figure_directory, ['PCA-' num2str(nPCs) '-']);

%%

%%

% cross-validation analysis
% pca_reliability(v, v12, v3, usubs_three_scans, si, [figure_directory 'crossval-PCA-']);
pca_reliability_mediancorr(v_mode_centered, v12_mode_centered, v3_mode_centered, usubs_three_scans, si, [figure_directory 'crossval-PCA-medcorr-mode-centered-']);
% pca_reliability_mediancorr(v1, v12, v3, usubs_three_scans, si, [figure_directory 'crossval-PCA-medcorr-scan1-']);
% pca_reliability_mediancorr(v12, v12, v3, usubs_three_scans, si, [figure_directory 'crossval-PCA-medcorr-scan12-']);

% Select first N PCs, rescale V matrix to have unit variance
nPCs = 6
[U,S,V] = svd(v_mode_centered,'econ');
Vr = V(:,1:nPCs)'*sqrt(size(V,1)-1);
Ur = U(:,1:nPCs)/sqrt(size(V,1)-1);
Sr = S(1:nPCs,1:nPCs);

% ICA optimization
[Vopt_best, h_pairs_best, negEnt_best, Vopt_all, h_pairs_all, negEnt_all] = ICA_NegEnt(Vr, 10, seed, 1);
% save(NegEntOpt_matfile, 'Vopt_best', 'h_pairs_best', 'negEnt_best', 'Vopt_all', 'h_pairs_all', 'negEnt_all');


% orient with maximally-skewed dimension
Vne = Vopt_best;
Une = Ur*Sr*Vr*pinv(Vne);
xi = sign(skewness(Vne'))==-1;
Une(:,xi) = -Une(:,xi);
Vne(xi,:) = -Vne(xi,:);


% plot figures
UneNorm = nan(size(Une));
for i = 1:nPCs
    UneNorm(:,i) = (Une(:,i) - min(Une(:,i)))/(max(Une(:,i))-min(Une(:,i)));
end

component_names = {'1','2','3','4','5','6'};

% category and response profiles
[R, C, conds] = naturalsound_acoustic_category_regressors('gammatone_Spec_SpecTempMod');
for i = 1:6
    response_profile_figure(UneNorm(:,i), C, conds,[0 1]); 
    export_fig([figure_directory 'ICA-mode-centered-response-profile-' component_names{i} '.pdf'],'-pdf');
    category_figures(UneNorm(:,i), C, R, conds, [figure_directory 'ICA-mode-centered-category-' component_names{i}]);
end

% plot histogram by subject
plothist_bysubject(Vne, si, [figure_directory 'ICA-mode-centered-' num2str(nPCs) '-NegEntOpt-hist-bysubject'], component_names);



%% Entropy-based ICA with Correlation Constraint

% ICA analysis
nPCs = 5;
nreps = 1000;
seed = 1;
NegEntOpt_matfile = [figure_directory 'NegEntOpt-' num2str(nPCs) '-' num2str(nreps) '_seed' num2str(seed) '.mat'];
if ~exist(NegEntOpt_matfile,'file')
    % Select first N PCs, rescale V matrix to have unit variance
    [U,S,V] = svd(v,'econ');
    Vr = V(:,1:nPCs)'*sqrt(size(V,1)-1);
    Ur = U(:,1:nPCs)/sqrt(size(V,1)-1);
    Sr = S(1:nPCs,1:nPCs);
    
    % ICA optimization
    [Vopt_best, h_pairs_best, negEnt_best, Vopt_all, h_pairs_all, negEnt_all] = ICA_NegEnt(Vr, nreps, seed, 1);
    save(NegEntOpt_matfile, 'Vopt_best', 'h_pairs_best', 'negEnt_best', 'Vopt_all', 'h_pairs_all', 'negEnt_all');
else
    [U,S,V] = svd(v,'econ');
    Vr = V(:,1:nPCs)'*sqrt(size(V,1)-1);
    Ur = U(:,1:nPCs)/sqrt(size(V,1)-1);
    Sr = S(1:nPCs,1:nPCs);
    load(NegEntOpt_matfile);
end

% orient with maximally-skewed dimension
Vne = Vopt_best;
Une = Ur*Sr*Vr*pinv(Vne);
xi = sign(skewness(Vne'))==-1;
Une(:,xi) = -Une(:,xi);
Vne(xi,:) = -Vne(xi,:);

% specific to particular analysis
Une(:,1) = -Une(:,1); % flipping pitch
Vne(1,:) = -Vne(1,:);
component_names = {'pitch','speech','tonotopy','music','envsounds'};
% component_names = {'envsounds', 'pitch', 'music', 'hfreq', 'lfreq', 'speech'};

% plot figures
UneNorm = nan(size(Une));
for i = 1:nPCs
    UneNorm(:,i) = (Une(:,i) - min(Une(:,i)))/(max(Une(:,i))-min(Une(:,i)));
end

% category and response profiles
[R, C, conds] = naturalsound_acoustic_category_regressors('gammatone_Spec_SpecTempMod');
for i = 1:5
    response_profile_figure(UneNorm(:,i), C, conds,[0 1]); 
    export_fig([figure_directory 'ICA-response-profile-' component_names{i} '.pdf'],'-pdf');
    category_figures(UneNorm(:,i), C, R, conds, [figure_directory 'ICA-category-' component_names{i}]);
end

naturalsound_acoustic_category_figures(UneNorm, figure_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-NormProfile-ExtraAcoustics-'],'corr_bounds',[-1 1]);

% negentropy ICA
plot_rotation_negentropy_v2(Vne,[-0.02 0.4]);
export_fig([figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');

% negentropy PCA
plot_rotation_negentropy_v2(Vr,[-0.02 0.4]);
export_fig([figure_directory 'PCA-' num2str(nPCs) '-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');

% negentropy gaussian
Vgauss = normrnd(0,1,size(Vr));
[Vopt_best_gauss] = ICA_NegEnt(Vgauss, 10, seed, 1);
plot_rotation_negentropy_v2(Vgauss,[-0.02 0.3]);
export_fig([figure_directory 'Gauss-' num2str(nPCs) '-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');
plot_rotation_negentropy_v2(Vopt_best_gauss,[-0.02 0.3]);
export_fig([figure_directory 'Gauss-NegEntOpt' num2str(nPCs) '-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');

% rotation figure for factorized distribution
ResetRandStream(1);
Vfac = nan(size(Vne));
for i = 1:size(Vne,1)
    Vfac(i,:) = Shuffle(Vne(i,:));
end
plot_rotation_negentropy_v2(Vfac,[-0.02 0.4]);
export_fig([figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt-NegEntPairs-Factorized.pdf'],'-nocrop','-transparent','-pdf');

% kurtosis
figure;
x = kurtosis(Vne')-3;
[~,xi] = sort(x,'descend');
mybar(x,component_names, component_names(xi), 'FaceColor',[1 1 1]);
ylabel('Kurtosis-3');
box off;
export_fig([figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt_kurtosis.pdf'],'-pdf','-nocrop','-transparent');

% scatter
myscatter2(Vne',si,'allpairs',3,[figure_directory 'ICA-' num2str(size(Vne,1)) '-NegEntOpt-scatter'],component_names);
mycontour(Vne','allpairs',[figure_directory 'ICA_contour'],component_names);

% plot histogram by subject
plothist_bysubject(Vne, si, [figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt-hist-bysubject'], component_names);

% entropy
figure;
x = nan(1,size(Vne,1));
for i = 1:size(Vne,1)
    x(i) = log(sqrt(2*pi*exp(1)))-entropy(Vne(i,:));
end
[~,xi] = sort(x,'descend');
mybar(x, component_names, component_names(xi));
ylabel('Negentropy');
export_fig([figure_directory 'ICA-' num2str(nPCs) '_entropy.pdf'],'-pdf','-nocrop','-transparent');

% variance
UneUnitNorm = nan(size(Une));
for i = 1:nPCs
    UneUnitNorm(:,i) = Une(:,i)/std(Une(:,i));
end
VneUnitNorm = pinv(UneUnitNorm)*v;
[~,xi] = sort(var(VneUnitNorm'),'descend');
mybar(var(VneUnitNorm'), component_names, component_names(xi));
ylabel('Variance');
export_fig([figure_directory 'ICA-' num2str(nPCs) '_variance.pdf'],'-pdf','-nocrop','-transparent');

% skewness
nsmps = 1000;
skew_smps = nan(nsmps,nPCs);
for i = 1:nsmps
    xi = randi(size(Vne,2), [1,size(Vne,2)]);
    skew_smps(i,:) = skewness(Vne(:,xi)');
end
skew_data = abs(skewness(Vne'));
[~,xi] = sort(skew_data,'descend');
mybar(skew_data, component_names, component_names(xi), 'errorbar', std(skew_smps));
ylabel('Skewness');
export_fig([figure_directory 'ICA-' num2str(nPCs) '_skewness.pdf'],'-pdf','-nocrop','-transparent');

% laterality
laterality = nan(length(usubs), nPCs); 

for i = 1:length(usubs)
    
    % compute voxel weights for different splits of the data
    if any(usubs(i) == usubs_three_scans)
        VneS1 = pinv(Une)*v1(:,si==usubs(i));
        VneS2 = pinv(Une)*v2(:,si==usubs(i));
        VneS3 = pinv(Une)*v2(:,si==usubs(i));
    else
        VneS1 = pinv(Une)*v1(:,si==usubs(i));
        VneS2 = pinv(Une)*v2(:,si==usubs(i));
    end
    
    for j = 1:nPCs
        
        % compute difference in laterality
        %         laterality_direction = mean(Vne(j,si == usubs(i) & hi==2)') - mean(Vne(j,si == usubs(i) & hi==1)');
        
        % calculate percent of explainable variance captured by laterality
        if any(usubs(i) == usubs_three_scans)
            y = [VneS1(j,:)', VneS2(j,:)', VneS3(j,:)'];
        else
            y = [VneS1(j,:)', VneS2(j,:)'];
        end
        [normalized_r, normalized_r2] = normalized_correlation(hi(si==usubs(i)), y);
        laterality(i,j) = normalized_r2 * sign(normalized_r) * 100;
        
    end
end

figure;
cols = distinguishable_colors(length(usubs),1);
bar(1:nPCs,median(laterality), 'FaceColor', [1 1 1]);
hold on;
for i = 1:length(usubs)
    plot(1:nPCs, laterality(i,:), 'o', 'LineWidth',2, 'Color', cols(i,:));
end
ylabel('Right - Left (% Explainable Variance)');
set(gca, 'XTickLabel', component_names);
box off;
export_fig([figure_directory 'laterality.pdf'],'-pdf','-nocrop','-transparent');

kurtosis_positive = nan(1,nPCs);
kurtosis_negative = nan(1,nPCs);
for i = 1:nPCs
    [N,bins] = hist(Vne(i,:)',100);
    [~,xi] = max(N);
    x = Vne(i,:)'-bins(xi);
    xpos = x(x>0);
    xneg = x(x<0);
    %     xpos = xpos/sqrt(mean(xpos.^2));
    %     xneg = xneg/sqrt(mean(xneg.^2));
    
    kurtosis_positive(i) = mean(xpos.^4)./mean(x.^2).^2;
    kurtosis_negative(i) = mean(xneg.^4)./mean(x.^2).^2;
    
    %     x = [xpos;xneg];
    %     figure;
    %     [y,bins] = hist(x,50);
    %     yn = normpdf(bins,0,1);
    %     plot(bins, y/sum(y),'r','LineWidth',2); hold on;
    %     plot(bins, yn/sum(yn),'k--','LineWidth',2);
    %     xlim([-6 6]);
end

for i = 1:nPCs
    
    
    bins = linspace(-5,5,75);
    x = Vne(i,:)';
    x = x/std(x);
    yh = hist(x,bins);
    figure;
%     plot(bins,yh/sum(yh), 'r'); hold on;
    h = normpdf(bins,0,0.5);
    yh = conv(yh,h/sum(h),'same');
%     plot(bins,yh/sum(yh), 'g');
    [~,xi] = max(yh);
    x = x - bins(xi);
    yh = hist(x,bins);
    plot(bins,yh/sum(yh), 'b'); hold on;
    yL = ylim;
    plot([0 0], yL,'k--');
    
    %     yh = hist(x,bins);
    %
    %     bins = linspace(-4,4,100)';
    %     figure;
    %     y = hist(Vne(i,:)' - median(Vne(i,:)),bins);
    %     y([1,end]) = [];
    %     bins([1,end]) = [];
    %     yn = normpdf(bins,0,1);
    %     plot(bins, y/sum(y),'r','LineWidth',2); hold on;
    %     plot(bins, yn/sum(yn),'k--','LineWidth',2);
%     xlim([- 5]);
end



%%

% ICA analysis
nPCs = 5;
nreps = 10;
seed = 1;

medcorr = nan(length(usubs_three_scans),nPCs);
for i = 1:length(usubs_three_scans)
    
    NegEntOpt_matfile = [figure_directory 'NegEntOpt-' num2str(nPCs) '-' num2str(nreps) '_seed' num2str(seed) '_leftout' num2str(usubs_three_scans(i)) '.mat'];
    if ~exist(NegEntOpt_matfile,'file')
        % Select first N PCs, rescale V matrix to have unit variance
        [U,S,V] = svd(v(:,~ismember(si,usubs_three_scans(i))),'econ');
        Vr = V(:,1:nPCs)'*sqrt(size(V,1)-1);
        Ur = U(:,1:nPCs)/sqrt(size(V,1)-1);
        Sr = S(1:nPCs,1:nPCs);
        
        % ICA optimization
        [Vopt_best, h_pairs_best, negEnt_best, Vopt_all, h_pairs_all, negEnt_all] = ICA_NegEnt(Vr, nreps, seed, 0);
        save(NegEntOpt_matfile, 'Vopt_best', 'h_pairs_best', 'negEnt_best', 'Vopt_all', 'h_pairs_all', 'negEnt_all');
    else
        [U,S,V] = svd(v(:,~ismember(si,usubs_three_scans(i))),'econ');
        Vr = V(:,1:nPCs)'*sqrt(size(V,1)-1);
        Ur = U(:,1:nPCs)/sqrt(size(V,1)-1);
        Sr = S(1:nPCs,1:nPCs);
        load(NegEntOpt_matfile);
    end
    
    % orient with maximally-skewed dimension
    VneX = Vopt_best;
    UneX = Ur*Sr*Vr*pinv(VneX);
    
    [~,xi] = max(abs(corr(UneX, Une)));
    UneX = UneX(:,xi);
    VneX = VneX(xi,:);
    
    for j = 1:nPCs
        scan12 = v12(:,ismember(si,usubs_three_scans(i)));
        scan3 = v3(:,ismember(si,usubs_three_scans(i)));
        X = UneX(:,~ismember(1:nPCs,j));
        prediction_IC_leftout = X*pinv(X)*scan12;
        prediction_IC_leftin = UneX*pinv(UneX)*scan12;
        
        medcorr(i,j) = median(fastcorr3(prediction_IC_leftin, scan3)) - median(fastcorr3(prediction_IC_leftout, scan3));
    end
end



%% Acoustic analysis

UneE1 = v2_threescans*pinv(pinv(Une)*v1_threescans);
UneE2 = v3_threescans*pinv(pinv(Une)*v1_threescans);
% [R2, C] = naturalsound_acoustic_category_regressors;
% acoustic_category_figures_v2([UneE1(:,1), UneE2(:,1)], C, R2, [figure_directory 'acoustic_category_ICA-']);

%% Hand-Tuned ICA based on maximally neg-entropic directions

HandTuneRotMat_file = [figure_directory 'HandTuneRotMat.mat'];
if ~exist(HandTuneRotMat_file,'file')
    HandTuneRotMat = ICA_NegEnt_HandTune(Vne);
    save(HandTuneRotMat_file, 'HandTuneRotMat');
else
    load(HandTuneRotMat_file);
end

n = size(HandTuneRotMat,1);
T = nan(n,n);
for i = 1:n
    T(i,:) = HandTuneRotMat(i,:,i);
end

Vhand = T*Vne;
Uhand = Une*pinv(T);
UhandNorm = nan(size(Uhand));
for i = 1:size(Uhand,2)
    UhandNorm(:,i) = (Uhand(:,i) - min(Uhand(:,i)))/(max(Uhand(:,i))-min(Uhand(:,i)));
end
naturalsound_acoustic_category_figures(UhandNorm, figure_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-NormProfile-HandTune-'],'corr_bounds',[-1 1]);

% negentropy ICA
plot_rotation_negentropy_v2(Vhand,[-0.02 0.3]);
export_fig([figure_directory 'ICA-HandTune-' num2str(nPCs) '-NegEntOpt-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');

% scatter
myscatter2(Vhand',si,'allpairs',3,[figure_directory 'ICA-HandTune-' num2str(size(Vne,1)) '-NegEntOpt-scatter'],component_names);

%% Variance explained by speech, music and language
conds = read_conditions('naturalsound',1,'main_v3_combined');
english = {'stim332_angry_shouting','stim174_girl_speaking','stim232_man_speaking','stim38_boy_speaking','stim414_woman_speaking'};
foreign = {'stim502_french', 'stim502_french', 'stim504_german', 'stim508_russian', 'stim503_italian'};
nonspeech = {'stim318_school_bell', 'stim51_car_accelerating','stim25_bike_bell','stim31_blender','stim488_busy_signal'};
music = {'stim268_piano','stim394_violin','stim298_contemporary_rock_song','tim254_orchestra_music','stim72_cello'};

english_response = mean(v(ismember(conds, english),:)); 
foreign_response = mean(v(ismember(conds, foreign),:)); 
nonspeech_response = mean(v(ismember(conds, nonspeech),:));
music_response = mean(v(ismember(conds, music),:));

contrast_names = {'Foreign-EnvSounds', 'Music-EnvSounds', 'English-Foreign'};
contrast_variance = [var(foreign_response-nonspeech_response), var(music_response-nonspeech_response), var(english_response-foreign_response)];
mybar(contrast_variance,contrast_names,contrast_names);
box off;
export_fig([figure_directory 'language_speech_music_explained_variance.pdf'],'-nocrop','-transparent','-pdf');

%%



%% Music experience

% NEED 122 and 123!!
yrs_practice = [1 10; 11 2; 106 15; 121 12; 143 10; 160 3; 161 0; 164 0];

%% Acoustic analyses using Shamma model

addpath(genpath([pwd '/nsltools']));
conds = read_conditions('naturalsound', 1, 'main_v3_combined');
stimulus_directory = [params('rootdir') 'naturalsound/stims/final_stims165/'];
rv = [0.5 1 2 4 8 16 32 64 128 256];
sv = [0.125 0.25 0.5 1 2 4 8];
paras(1) = 1; % 500 frames per second
paras(2) = 8; % leaky integrator with 8 ms time constant
paras(3) = 0.1; % compression factor
paras(4) = 0; % specifies sampling rate as 16*2^(paras(4))

nr = length(rv);
ns = length(sv);
loadload; close all;

cortical_model_file = [figure_directory 'cortical_model-Hash' DataHash([rv,sv,paras]) '.mat'];
if ~exist(cortical_model_file, 'file')
    rs = nan(ns,nr*2,165);
    for i = 1:165
        [wav,sr] = wavread([stimulus_directory conds{i} '.wav']);
        wav_16kHz = resample(wav, 16000, sr);
        wav_unitseq = unitseq(wav_16kHz);
        y = wav2aud(wav_unitseq, paras);
        
        fcorname = [params('rootdir') 'naturalsound/stims/final_stims165-cor/' conds{i} '.cor'];
        cr = aud2cor(y, paras, rv, sv, fcorname);
        rs(:,:,i) = mean(mean(abs(cr), 4),3);
    end
    save(cortical_model_file, 'rs');
else
    load(cortical_model_file);
end

% correlation with all rate/scale variables
rs_collapsed = (rs(:,1:nr,:) + rs(:,(1:nr)+nr,:))/2;
rs_corr = nan(ns,nr,nPCs);
for i = 1:nPCs
    for j = 1:ns
        for k = 1:nr
            x = squeeze(rs_collapsed(j,k,:));
            rs_corr(j,k,i) = mean([corr(x, UneE1(:,i)),corr(x, UneE2(:,i))])/sqrt(corr(UneE1(:,i),UneE2(:,i)));
        end
    end
end

for i = 1:nPCs
    figure;
    imagesc(flipud(rs_corr(:,:,i)),[-1 1]);
    set(gca,'XTick',1:nr,'XTickLabel',rv,'YTick',1:ns,'YTickLabel',fliplr(sv));
    xlabel('Temporal Rate (Hz)','FontSize',20); ylabel('Spectral Scale (cyc / oct)');
    colorbar;
    saveas(gcf,[figure_directory 'ICA_cortical_model_' component_names{i} '.png']);
end

% Acoustics Just-Based on Cortical Model
Rcortical.stat_names = cell(1,ns*nr);
Rcortical.stats = nan(165,ns*nr);
Rcortical.ids = R2.ids;
Rcortical.groups = nan(1,ns*nr);
Rcortical.regressors = 1:ns*nr;
for j = 1:nr
    for i = 1:ns
        xi = i + (j-1)*ns;
        Rcortical.stat_names{xi} = ['scale' num2str(sv(i)) '-rate' num2str(rv(j))];
        Rcortical.stats(:,xi) = squeeze(rs_collapsed(i,j,:));
        Rcortical.groups(xi) = j;
    end
end

R6 = R2;
R6.stat_names = [R6.stat_names, Rcortical.stat_names];
R6.stats = [R6.stats, Rcortical.stats];
R6.groups = [R6.groups, Rcortical.groups+6];
R6.regressors = [R6.regressors, Rcortical.regressors+23];

% acoustic correlations
for i = 1:nPCs
    acoustic_figures_normalized_correlation([UneE1(:,i), UneE2(:,i)], R6, [figure_directory 'ICA-acoustic-kitchensink-' component_names{i}],[-1 1]);
end

% acoustic correlations
for i = 1:5
    acoustic_figures_normalized_correlation([UneE1(:,i), UneE2(:,i)], Rcortical, [figure_directory 'ICA-cortical-' component_names{i}],[-1 1]);
end
acoustic_figures_normalized_correlation([v1_grandmean, v2_grandmean], Rcortical, [figure_directory 'grandmean-cortical-' component_names{i}],[-1 1]);

% marginals
scale_marginals = squeeze(mean(rs_collapsed,2))';
rate_marginals = squeeze(mean(rs_collapsed,1))';

scale_names = cell(1,ns);
for i = 1:ns
    scale_names{i} = ['scale' num2str(sv(i))];
end

rate_names = cell(1,nr);
for i = 1:nr
    rate_names{i} = ['rate' num2str(rv(i))];
end

% combined structure that includes statistics from Josh and Shamma model
freqi = 2:5;
resi = 20;
cochi = 22:23;
envi = 1;
ampmodi = 6:11;
specmodi = 12:19;
R4.stat_names = [R2.stat_names(envi), R2.stat_names(freqi), rate_names, scale_names, R2.stat_names(resi), R2.stat_names(cochi)]; %, R2.stat_names(ampmodi), R2.stat_names(specmodi)
R4.stats = [R2.stats(:,envi), R2.stats(:,freqi), rate_marginals, scale_marginals, R2.stats(:,resi), R2.stats(:,cochi)]; % , R2.stats(:,ampmodi), R2.stats(:,specmodi)
R4.ids = R2.ids;
R4.groups = [1,2*ones(1,4),3*ones(1,nr),4*ones(1,ns), 5, 6*ones(1,2)]; %, 7*ones(1,6), 8*ones(1,8)
R4.regressors = find(~ismember(R4.stat_names,{'3200','scale8'})); % ,'AmpVar','SpecVar'

% acoustic correlations across stimuli
figure;
x = corr(R4.stats);
imagesc(x.*sign(x),[-1 1]);
set(gca,'XTick',1:length(R4.stat_names),'XTickLabel',R4.stat_names,'YTick',1:length(R4.stat_names),'YTickLabel',R4.stat_names);
rotateticklabel(gca,45,-20.25,14,'Demi');
colorbar;
saveas(gcf,[figure_directory 'acoustic-correlations.png']);

R2 = naturalsound_acoustic_category_regressors('shamma');

%% Correlation partialling out spectrum response

R = naturalsound_acoustic_category_regressors('gammatone_Spec');
[U,S] = svd(zscore(R.stats),'econ');
exvar = cumsum(diag(S).^2)/sum(diag(S).^2);
xi = find(exvar>0.95,1);
spectrum_regressors = U(:,1:xi);

UneE1_SpectrumResidual = nan(size(UneE1));
UneE2_SpectrumResidual = nan(size(UneE1));
for i = 1:nPCs
    response = zscore(UneE1(:,i));
    prediction = nan(165,1);
    for k = 1:165;
        xi = setdiff(1:165,k);
        prediction(k) = spectrum_regressors(k,:) * pinv(spectrum_regressors(xi,:)) * response(xi);
    end
    UneE1_SpectrumResidual(:,i) = response-prediction;
    
    response = zscore(UneE2(:,i));
    prediction = nan(165,1);
    for k = 1:165;
        xi = setdiff(1:165,k);
        prediction(k) = spectrum_regressors(k,:) * pinv(spectrum_regressors(xi,:)) * response(xi);
    end
    UneE2_SpectrumResidual(:,i) = response-prediction;
end

R = naturalsound_acoustic_category_regressors('gammatone_Spec_SpecTempModMarg');
% acoustic correlations with 
for i = 1:nPCs
    acoustic_figures_normalized_correlation([UneE1_SpectrumResidual(:,i), UneE2_SpectrumResidual(:,i)], R, [figure_directory 'ICA-acoustics-spectrum-partialled-' component_names{i}],[-1 1]);
end

%% Partialling out spectrum from regressors

clear R;
R(1) = naturalsound_acoustic_category_regressors('gammatone_Spec');
R(2) = naturalsound_acoustic_category_regressors('gammatone_SpecTempModMarg');
X = zscore(R(1).stats(:,1:6));
for i = 1:length(R(2).stat_names)
    y = zscore(R(2).stats(:,i));
    R(2).stats(:,i) = y - X*pinv(X)*y;
end

% acoustic correlations with 
for i = 1:nPCs
    acoustic_figures_normalized_correlation([UneE1(:,i), UneE2(:,i)], R(2), [figure_directory 'ICA-acoustics-spectrum-partialled-' component_names{i}],[-1 1]);
end

%% Spectrotemporal modulation partialling out spectrum

R(1) = naturalsound_acoustic_category_regressors('gammatone_Spec');
R(2) = naturalsound_acoustic_category_regressors('gammatone_SpecTempMod');
X = zscore(R(1).stats(:,1:6));
for i = 1:length(R(2).stat_names)
    y = zscore(R(2).stats(:,i));
    R(2).stats(:,i) = y - X*pinv(X)*y;
end
    
for i = 1:nPCs
    denominator = sqrt(corr(UneE1(:,i), UneE2(:,i)));
    r = mean(corr([UneE1(:,i), UneE2(:,i)], R(2).stats)) / denominator;
    rmatrix = reshape(r, [7,9]);
    figure;
    imagesc(flipud(rmatrix),[-0.7 0.7]);
    box off;
    print('-dpng',[figure_directory 'spectrotemporal-modulation-matrix-' component_names{i}],'-r200');
end

%% Acoustic analysis

% acoustic_type = 'gammatone_Spec';
% acoustic_type = 'gammatone_Spec_AmpMod';
acoustic_type = 'gammatone_Spec_SpecTempMod';
% acoustic_type = 'gammatone_Spec_AmpMod_CochCorr';
[R, C] = naturalsound_acoustic_category_regressors(acoustic_type);

% acoustic correlations with 
for i = 1:nPCs
    acoustic_figures_normalized_correlation([UneE1(:,i), UneE2(:,i)], R, [figure_directory 'ICA-acoustics-' acoustic_type '-' component_names{i}],[-1 1]);
end
acoustic_figures_normalized_correlation([v1_grandmean, v2_grandmean], R, [figure_directory 'ICA-acoustics-' acoustic_type '-grandmean'],[-1 1]);

% comparing variance captured by acoustic and category models
r2 = nan(3,nPCs);
for i = 1:nPCs
    [r2(:,i)] = acoustic_vs_category_explained_variance([UneE1(:,i), UneE2(:,i)], C, R);
end

figure;
hold on;
order = [3 1 5 4 2];
cols = lines(3);
for i = 1:3
    y = r2(i,order);
    err = r2rerr(sqrt(y),165).^2;
    errorbar(1:nPCs, 100*y, 100*(y-err(1,:)), 100*(err(2,:)-y), '-o', 'Color',cols(i,:),'LineWidth',2);
end
legend('acoustic','category','acoustic+category','Location','Best');
ylim([0 100]); xlim([0 nPCs+1]);
set(gca,'XTick',1:nPCs,'XTickLabel',component_names(order));
ylabel('% Explainable Variance');
box off;
export_fig([figure_directory 'acoustic_vs_category_exvar_lineplot_' acoustic_type '.pdf'],'-pdf','-transparent','-nocrop');

% figure comparison of acoustic and category explained variance
figure;
set(gcf,'Position',[0 0 800 600]);
set(gca,'FontSize',14);
cols = lines(3);
order = [3 1 5 4 2];
model_names = {'acoustic','category','combined'};
bar_names = cell(1,nPCs*4);
ticks = false(size(bar_names));
for i = 1:nPCs
    for j = 1:3
        xi = j+(i-1)*4;
        y = r2(j,order(i));
        bar(xi,100*y,'FaceColor',cols(j,:));
        if i == 1 && j == 1
            hold on;
        end
        err = r2rerr(sqrt(y),165).^2;
        errorbar(xi,100*y, 100*(y-err(1)), 100*(err(2)-y), 'k')
        bar_names{xi} = [component_names{order(i)} '-' model_names{j}];
        ticks(xi) = true;
    end
end
set(gca,'Position',[0.15 0.25 0.7 0.65]);
ylim([0 100]); xlim([0,length(bar_names)+1]);
set(gca,'XTick',find(ticks),'XTickLabel',bar_names(ticks));
rotateticklabel(gca,45,0.5,14,'Normal');
ylabel('% Explainable Variance');
box off;
export_fig([figure_directory 'acoustic_vs_category_exvar_' acoustic_type '.pdf'],'-pdf','-transparent','-nocrop');

acoustic_models = {'gammatone_Spec','gammatone_Spec_SpecMod','gammatone_Spec_AmpMod','gammatone_Spec_AmpMod_CochCorr','gammatone_Spec_SpecTempMod','gammatone_Spec_SpecTempMod_CochCorr'};
acoustic_models = {'gammatone_Spec','gammatone_Spec_AmpMod','gammatone_Spec_SpecMod','gammatone_Spec_SpecTempModMarg','gammatone_Spec_SpecTempMod','gammatone_Spec_SpecTempModMarg_CochCorr'};
acoustic_models = {'gammatone_Spec','gammatone_Spec_AmpMod','gammatone_Spec_SpecMod','gammatone_Spec_SpecTempMod'};

clear R;
for i = 1:length(acoustic_models)
    R(i) = naturalsound_acoustic_category_regressors(acoustic_models{i});
end

% comparing variance captured by acoustic and category models
r2 = nan(length(acoustic_models),nPCs);
for j = 1:length(R)
    for i = 1:nPCs
        %         [r2(j,i)] = linear_explained_variance_3fold_crossvalidation([UneE1(:,i), UneE2(:,i)], R(j));
        [r2(j,i)] = linear_explained_variance([UneE1(:,i), UneE2(:,i)], R(j));
    end
end

models = {'gammatone_Spec','gammatone_Spec_AmpMod','gammatone_Spec_SpecMod','gammatone_Spec_SpecTempMod','category','category+SpecTempMod'};

clear R;
for i = 1:length(models)
    switch models{i}
        case 'category'
            R(i).stats = C.continuous_scores;
        case 'category+SpecTempMod'
            X = naturalsound_acoustic_category_regressors('gammatone_Spec_SpecTempMod');
            R(i).stats = [X.stats, C.continuous_scores];
        otherwise
            X = naturalsound_acoustic_category_regressors(models{i});
            R(i).stats = X.stats;
    end
end

% comparing variance captured by acoustic and category models
r2 = nan(length(R),nPCs);
for j = 1:length(R)
    for i = 1:nPCs
        %         [r2(j,i)] = linear_explained_variance_3fold_crossvalidation([UneE1(:,i), UneE2(:,i)], R(j));
        [r2(j,i)] = linear_explained_variance([UneE1(:,i), UneE2(:,i)], R(j));
    end
end

% figure comparison of acoustic and category explained variance
figure;
set(gcf,'Position',[0 0 800 600]);
set(gca,'FontSize',10);
hold on;
cols = lines(length(models));
order = [3 1 5 4 2];
nspaces_percomponent = (length(models)+1);
bar_names = cell(1,nPCs*nspaces_percomponent);
ticks = false(size(bar_names));
for i = 1:nPCs
    for j = 1:length(models)
        xi = j+(i-1)*nspaces_percomponent;
        y = r2(j,order(i));
        bar(xi,100*y,'FaceColor',cols(j,:));
        if i == 1 && j == 1
            hold on;
        end
        err = r2rerr(sqrt(y),165).^2;
        errorbar(xi,100*y, 100*(y-err(1)), 100*(err(2)-y), 'k')
        bar_names{xi} = [component_names{order(i)} '-' strrep(strrep(models{j},'_','-'),'gammatone-','')];
        ticks(xi) = true;
    end
end
set(gca,'Position',[0.15 0.35 0.7 0.55]);
ylim([0 100]); xlim([0,length(bar_names)+1]);
set(gca,'XTick',find(ticks),'XTickLabel',bar_names(ticks));
rotateticklabel(gca,45,0.5,10,'Normal');
ylabel('% Explainable Variance');
box off;
export_fig([figure_directory 'model_comparison_spectrotemporal_category.pdf'],'-pdf','-transparent','-nocrop');

% % acoustic correlations with 
% for i = 1:nPCs
%     acoustic_figures_normalized_correlation([UneE1(:,i), UneE2(:,i)], R2, [figure_directory 'ICA-acoustics-v2-' component_names{i}],[-1 1]);
% end

% regression analysis
xi = find(~ismember(R4.stat_names,{'400','800','scale8'}));
R5.stat_names = R4.stat_names(xi);
R5.stats = R4.stats(:,xi);
R5.ids = R4.ids;
R5.groups = R4.groups(:,xi);
R5.regressors = R5.regressors(1:length(R5.stat_names));

% acoustic regression
% for i = 1:5
%     acoustic_figures_regression(Une(:,i), R5, [figure_directory 'ICA-acoustic-regression-' component_names{i}]);
% end

[R2, C, conds, R3] = naturalsound_acoustic_category_regressors;



% residuals of ICA components
Uresid = nan(size(Une));
for i = 1:5
    y = Une(:,i);
    X = Une(:,setdiff(1:5,i));
    Uresid(:,i) = y - X*pinv(X)*y;
end

% acoustic correlations with 
for i = 1:5
    acoustic_figures_normalized_correlation([Uresid(:,i),Uresid(:,i)], R4, [figure_directory 'ICA-resid-acoustics-' component_names{i}],[-1 1]);
end

% acoustic correlations with 
for i = 1:5
    category_figures(UneNorm(:,i), C, R4, [figure_directory 'ICA-category-' component_names{i}],[-1 1]);
end
acoustic_figures_normalized_correlation(UneNorm, R4, [figure_directory 'grandmean-acoustics-' component_names{i}],[-1 1]);

%% Subcortical nuclei
register_handpick_roiAll({'IC_2mm','MGN_2mm'},'naturalsound',usubs,'runtypes',{'main_v3_combined'});

usubs_IC = [1 11 121 123 143 160];
psc_scans = nan(3,165,length(usubs_IC));
for i = 1:length(usubs_IC)
    [~,~,psc_scans(:,:,i)] = roi_pscAll(usubs_IC(i), 'IC_2mm', 'naturalsound', 'main_v3_combined',  'sigav_plateau', 5, 'withsub_sem','noplot');
end
psc_scans_mean = squeeze(mean(psc_scans,3))';
acoustic_figures_normalized_correlation(mean(psc_scans_mean,2)*[1 1], R4, [figure_directory 'Colliculus-acoustics'],[-1 1]);
% acoustic_figures_normalized_correlation(psc_scans_mean, R2, [figure_directory 'Colliculus-acoustics-v2'],[-1 1]);

usubs_MGN = [1 11 121 123 143];
psc_scans = nan(3,165,length(usubs_MGN));
for i = 1:length(usubs_MGN)
    [~,~,psc_scans(:,:,i)] = roi_pscAll(usubs_MGN(i), 'MGN_2mm', 'naturalsound', 'main_v3_combined',  'sigav_plateau', 5, 'withsub_sem','noplot');
end
psc_scans_mean = squeeze(mean(psc_scans,3))';
acoustic_figures_normalized_correlation(mean(psc_scans_mean,2)*[1 1], R4, [figure_directory 'MGN-acoustics'],[-1 1]);

naturalsound_acoustic_category_figures(psc', [params('rootdir') 'naturalsound/summary_figures/'], ['IC-'],'corr_bounds',[-1 1]);



%% Project ICs back to grid

% parameters
min_subjects = 3;
fwhm_surf = 5; % additional smoothing
nsurfpts = 163842;

% directory and file information
fname_prefix = 'ICA-voxel_weights-';
surface_projections_directory = [figure_directory 'surface_maps/'];
if ~exist(surface_projections_directory,'file')
    mkdir(surface_projections_directory);
end

% initialize variables
voxel_weights_surface = nan(nPCs, nsurfpts*2, length(usubs));
norm_coeff_surface = nan(nsurfpts*2, length(usubs));
NaN_surface = nan(nsurfpts*2, length(usubs));

% ICA voxel weights
voxel_weights = Vne;

% loop through components
for k = 1:size(voxel_weights,1)
    
    % n_voxels x subject voxel weight matrix
    x = nan(n_voxels*length(usubs),1);
    x(voxel_selection) = voxel_weights(k,:)';
    voxel_weights_reshaped = reshape(x, [n_voxels, length(usubs)]);
    
    x = nan(n_voxels*length(usubs),1);
    x(voxel_selection) = norm_coeff;
    norm_coeff_reshaped = reshape(x, [n_voxels, length(usubs)]);
    
    hemis = {'rh','lh'};
    for i = 1:length(hemis)
        
        % convert to grid: xgrid x ygrid x n_subjects
        switch hemis{i}
            case 'rh'
                voxel_weights_gridded = reshape(voxel_weights_reshaped(1:n_voxels_rh, :), [size(grid_x{1}), length(usubs)]);
                norm_coeff_gridded = reshape(norm_coeff_reshaped(1:n_voxels_rh, :), [size(grid_x{1}), length(usubs)]);
                
            case 'lh'
                voxel_weights_gridded = reshape(voxel_weights_reshaped((1:n_voxels_lh) + n_voxels_rh, :), [size(grid_x{2}), length(usubs)]);
                norm_coeff_gridded = reshape(norm_coeff_reshaped((1:n_voxels_lh) + n_voxels_rh, :), [size(grid_x{2}), length(usubs)]);
        end
        
        % gaussian kernel for smoothing
        x = round(3*fwhm_surf/2);
        x = x-mod(x,2)+1; % make odd
        gaussian_kernel = fspecial('gaussian', [x,x], fwhm2std(fwhm_surf/2));
        gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
        
        
        % interpolate to surface
        for j = 1:length(usubs)
            
            %             if fwhm_surf > 0
            voxel_weights_gridded_smoothed = conv2_setNaNs_tozero(voxel_weights_gridded(:,:,j), gaussian_kernel);
            norm_coeff_gridded_smoothed = conv2_setNaNs_tozero(norm_coeff_gridded(:,:,j), gaussian_kernel);
            %             else
            %                 voxel_weights_gridded_smoothed = interpNaNs_by_convolution(:,:,j);
            %                 norm_coeff_gridded_smoothed = interpNaNs_by_convolution(:,:,j);
            %             end
            NaN_smooth = double(~isnan(voxel_weights_gridded(:,:,j)));

            %             figure;
            %             subplot(1,2,1);
            %             imagesc(voxel_weights_gridded(:,:,j));
            %             subplot(1,2,2);
            %             imagesc(voxel_weights_gridded_smoothed);
            
            voxel_weights_surface(k,vi{i} + 1 + (i-1)*nsurfpts,j) = interp2(grid_x{i},grid_y{i},voxel_weights_gridded_smoothed,vras{i}(:,1),vras{i}(:,2));
            norm_coeff_surface(vi{i} + 1 + (i-1)*nsurfpts,j) = interp2(grid_x{i},grid_y{i},norm_coeff_gridded_smoothed,vras{i}(:,1),vras{i}(:,2));
            NaN_surface(vi{i} + 1 + (i-1)*nsurfpts,j) = interp2(grid_x{i},grid_y{i},NaN_smooth,vras{i}(:,1),vras{i}(:,2));
            
            %             voxel_selection_surface(vi{i} + 1 + (i-1)*nsurfpts,j) = interp2(grid_x{i},grid_y{i},voxel_selection_gridded(:,:,j),vras{i}(:,1),vras{i}(:,2));
            %             x = interp2(grid_x{i},grid_y{i},voxel_weights_gridded(:,:,j),vras{i}(:,1),vras{i}(:,2));
            %             figure;
            %             scatter(vras{i}(:,1),vras{i}(:,2),3,x);
            %             figure;
            %             surf(grid_x{i},grid_y{i},voxel_weights_gridded(:,:,j));
            %             view(0,90);
        end
    end
end

% average across subjects
% voxel_weights_group = nanmean(voxel_weights_surface,3);
% voxel_weights_group(isnan(voxel_weights_group)) = 0;
% voxel_weights_group = mean(voxel_weights_surface,3);
% voxel_weights_group_norm = voxel_weights_group ./ (ones(size(voxel_weights_group,1),1)*sum(voxel_weights_group,1));
% voxel_weights_group_norm(isnan(voxel_weights_group_norm)) = 0;

% group_voxelthresh = sum(NaN_surface,2) > min_subjects - 1e-3;
% voxel_weights_group_norm(:,~group_voxelthresh) = 0;
% weight_std = std(voxel_weights_group_norm(:,group_voxelthresh)');
% 
% % histogram
% figure;
% [y,bins] = hist(voxel_weights_group_norm(1:nPCs,group_voxelthresh)',100);
% plot(bins,y/sum(y(:,1)),'LineWidth',2);
% legend(component_names{:});
% export_fig([surface_projections_directory 'histogram.pdf'],'-pdf','-nocrop','-transparent');

% weighted average across subjects
voxel_weights_group_weighted = nan(size(voxel_weights_surface(:,:,1)));
for i = 1:nPCs
    x = squeeze(voxel_weights_surface(i,:,:)) .* norm_coeff_surface;
    voxel_weights_group_weighted(i,:) =  nansum(x,2) ./ nansum(norm_coeff_surface,2);
end
voxel_weights_group_weighted(isnan(voxel_weights_group_weighted)) = 0;

% threshold by number of voxels included at that vertex
group_voxelthresh = sum(NaN_surface,2) > min_subjects - 1e-3;
voxel_weights_group_weighted(:,~group_voxelthresh) = 0;
weight_std = std(voxel_weights_group_weighted(:,group_voxelthresh)');

% histogram
figure;
[y,bins] = hist(voxel_weights_group_weighted(1:nPCs,group_voxelthresh)',100);
plot(bins,y/sum(y(:,1)),'LineWidth',2);
legend(component_names{:});
export_fig([surface_projections_directory 'histogram.pdf'],'-pdf','-nocrop','-transparent');

% surface directory
surf_directory = [params('rootdir') 'freesurfer/fsaverage/linear_analysis_surf/naturalsound_us' sprintf('-%d',usubs) '/'];
if ~exist(surf_directory,'dir')
    mkdir(surf_directory);
end

% group results
for k = 1
    
    name = component_names{k};

    for i = 1
        % group beta figure
        fname = [surf_directory fname_prefix '-beta-minsubjects' num2str(min_subjects) '-' name '-' hemis{i} '.mgz'];
        delete(fname);
        x = voxel_weights_group_weighted(k,(1:nsurfpts) + (i-1)*nsurfpts);
        xi = x~=0;
        xn = x;
        bounds = [-2 2];
        xn(x<bounds(1)*weight_std(k)) = 1;
        xn(x>bounds(2)*weight_std(k)) = 2;
        xn(x>bounds(1)*weight_std(k) & x<bounds(2)*weight_std(k) & xi) = interp1(weight_std(k)*bounds,[1 2], x(x>bounds(1)*weight_std(k) & x<bounds(2)*weight_std(k) & xi));
        MRIwrite_SNH(xn, fname, hemis{i});
        
        %         % for plotting
        %         t = 3;
        %         x = ica_weights_group(k,:,i)/weight_std(k);
        %         x(x<-t) = -t;
        %         x(x>t) = t;
        %         ica_normalized_forplotting = interp1([-2,-0.4,1.2,2]*(t/2),[1 2.66, 4.33, 6], x);
        %         ica_normalized_forplotting(~group_voxelthresh) = 0;
        %
        %         % group beta figure
        %         fname = [surf_directory idstring 'beta-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-component' num2str(k) '-' hemis{i} '-forplotting.mgz'];
        %         delete(fname);
        %         MRIwrite_SNH(ica_normalized_forplotting, fname, hemis{i});
        
        % save screenshot
        surface_image = [surface_projections_directory idstring 'beta-minsubjects' num2str(min_subjects) '-' name '-group-' hemis{i} '.png'];
        if true || ~exist(surface_image, 'file')
            %             freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.00001 0.00002 weight_std(k)*2],'piecewise');
            %             freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[0.00001 0.00002 weight_std(k)*2],'piecewise');
            %             freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.25 1.25 2]);
            freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[1 1.5 2]);
        end
    end
end

nsurfpts = 163842;
vnorm_recon = interp2(grid_x{1},grid_y{1},ic_group{i,1},vras{1}(:,1),vras{1}(:,2));
vnorm_recon(isnan(vnorm_recon)) = 0;
surf_recon = zeros(nsurfpts,1);
surf_recon(vi{1}+1) = vnorm_recon;
MRIwrite_SNH(surf_recon, [pwd '/test.mgz'], 'rh');
freeview3('fsaverage','rh','overlay',[pwd '/test.mgz'],'overlay_threshold',[0.25 1.25 2]);

figure;
scatter(vras{1}(:,1),vras{1}(:,2),3,vnorm_recon);



% grid directory
grid_directory = [figure_directory 'grids/'];
if ~exist(grid_directory, 'dir')
    mkdir(grid_directory);
end

% plot ic grids
close all;
for i = 1:size(Vica,1);
    %     figure(1);clf(1);
    %     set(gcf,'Position',[0 0 1000, 500])
    %     subplot(1,2,1);
    %     imagesc(fliplr(ic_group{i,2}),[-2 2]);
    %     %     surf(grid_x{2}, grid_y{2}, ic_group{i,2});
    %     %     view(0,90);
    %     title('Left Hemisphere');
    %     subplot(1,2,2);
    %     imagesc(fliplr(ic_group{i,1}),[-2 2]);
    %     %     surf(grid_x{1}, grid_y{1}, ic_group{i,1});
    %     %     view(0,90);
    %     title('Right Hemisphere');
    %     drawnow;
    %     export_fig([grid_directory 'ic-' num2str(i) '-group-vras-coordinates.png'],'-png');
    %     %     print(gcf, '-dpng', [grid_directory 'ic-' num2str(i) '-group-vras-coordinates.png'], '-r100');
    %
    close all;
    figure(1);clf(1);
    set(gcf,'Position',[0 0 1000, 500]);
    subplot(1,2,1);
    imagesc(rot90(ic_group{i,2},3),[-2 2]);
    title('Left Hemisphere');
    subplot(1,2,2);
    imagesc(rot90(ic_group{i,1},-1),[-2 2]);
    title('Right Hemisphere');
    %     export_fig([grid_directory 'ic-' num2str(i) '-group.png'],'-png');
    %     print(gcf, '-dpng', [grid_directory 'ic-' num2str(i) '-group.png'], '-r100');
    
    %     for j = 1:length(usubs)
    %         close all;
    %         figure(1);clf(1);
    %         set(gcf,'Position',[0 0 1000, 500]);
    %         subplot(1,2,1);
    %         imagesc(rot90(ic{i,2}(:,:,j),3),[-2 2]);
    %         %     surf(grid_x{2}, grid_y{2}, ic_group{i,2});
    %         title('Left Hemisphere');
    %         subplot(1,2,2);
    %         imagesc(rot90(ic{i,1}(:,:,j),-1),[-2 2]);
    %         title('Right Hemisphere');
    %         drawnow;
    %         %         print(gcf, '-dpng', [grid_directory 'ic-' num2str(i) '-us' num2str(usubs(j)) '.png'], '-r100');
    %         export_fig([grid_directory 'ic-' num2str(i) '-us' num2str(usubs(j)) '.png'],'-png');
    %     end
end


% % voxel selection, can optionally change criteria for voxels used to plot data, currently using the orignal 
% voxel_selection_forplotting = false(size(soundP));
% xi = ismember(si_unfilt, usubs_two_scans);
% voxel_selection_forplotting(xi) = all(soundP_individual_scans(xi,1:2)' > soundP_thresh) & regressExVar(xi) > regress;
% xi = ismember(si_unfilt, usubs_three_scans);
% voxel_selection_forplotting(xi) = all(soundP_individual_scans(xi,1:3)' > soundP_thresh) & regressExVar(xi) > 30;
% fprintf('%d total voxels (%.1f%%)\n', sum(voxel_selection_forplotting), 100*sum(voxel_selection_forplotting)/sum(all(~isnan(v_unfilt_allscans(:,:,1)))));

% ICA voxel weights
% v_forplotting = v_unfilt_mean(:,voxel_selection_forplotting);
% norm_coeff = sum(v_forplotting,1);
% v_norm = v_forplotting ./ (ones(size(v_forplotting,1),1)*norm_coeff);
% v_norm_demean = v_norm - grandmean*ones(1,size(v_forplotting,2));
% voxel_weights = pinv(Une)*v_norm_demean;

%% Regress IC components in whole brain

idstring = ['Une-' num2str(size(Une,2)) '-grandmean-'];
surface_projections_directory = [figure_directory 'surface_projections/'];
if ~exist(surface_projections_directory,'file')
    mkdir(surface_projections_directory);
end

nsurfpts = 163842;

% initialize variables
sound_responseP = zeros(nsurfpts*2, length(usubs));
ica_weights_soundthresh = zeros(nPCs+1, nsurfpts*2, length(usubs));
ica_logP = zeros(nPCs+1, nsurfpts*2, length(usubs));
normalization_coefficients = zeros(nsurfpts*2, length(usubs)); 

% sound threshold
soundthresh = 3;
min_subjects = 4;

% subjects
usubs = [1 11 106 121 122 123 143 160 161 164];
hemis = {'rh','lh'};

whole_brain_projection_file = [figure_directory 'whole_brain_projection-' DataHash({usubs, soundthresh, 5}) '.mat'];
if ~exist(whole_brain_projection_file, 'file')
    for j = 1:length(usubs)
            
        % subject id
        fprintf('us %d\n', usubs(j)); drawnow;
        subjid = ['naturalsound_us' num2str(usubs(j))];
        
        % surface directory
        surf_directory = [params('rootdir') 'freesurfer/fsaverage/linear_analysis_surf/' subjid '/'];
        if ~exist(surf_directory,'dir')
            mkdir(surf_directory);
        end
        
        % runs for that subject
        runs = read_runs('naturalsound',usubs(j),'main_v3_combined');
        
        % hemispheres
        sv = zeros(165, nsurfpts*2);
        for i = 1:2
            
            % auditory cortex mask
            mask = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
            
            % voxel responses
            for q = 1:length(runs)
                x = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r' num2str(runs(q)) '/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp3-6.mgz']);
                dat = squeeze(x.vol)';
                sv(:,mask.vnums+1 + (i-1)*nsurfpts) = sv(:,mask.vnums+1 + (i-1)*nsurfpts) + dat(:,mask.vnums+1)/length(runs);
            end
        end
        
        % sound response
        [~,p] = ttest(sv,0,0.05,'right');
        % [~,p] = ttest(sv,0,'tail','right');
        sound_responseP(:,j) = -log10(p);
        
        % select voxels with a significant sound response
        sv_selected = sv(:,sound_responseP(:,j) > soundthresh);
        norm_coeff = sum(sv_selected,1);
        sv_selected_norm = sv_selected ./ (ones(size(sv_selected,1),1)*norm_coeff);
        sv_selected_norm_demean = sv_selected_norm - grandmean*ones(1,size(sv_selected,2));
        
        % Regression analysis
        Uica_norm = Une ./ (ones(165,1)*std(Une));
        Uica_norm = [Une ./ (ones(165,1)*std(Une)), grandmean/std(grandmean)];
        [b, ~, logP, ~] = regress_SNH([Uica_norm], sv_selected);
        b_without_norm = b;
        %         b_without_norm = b .* (ones(nPCs,1)*norm_coeff);
        %         [b, ~, logP, ~] = regress_SNH([Uica_norm, grandmean/rms(grandmean)], sv_selected);
        
        % assign weights and compute standard deviation across voxels (used for plotting)
        ica_weights_soundthresh(:,sound_responseP(:,j) > soundthresh,j) = b_without_norm;
        weight_std = std(b_without_norm');
        
        % assign p-values
        ica_logP(:,sound_responseP(:,j) > soundthresh,j) = logP;
        
        normalization_coefficients(sound_responseP(:,j) > soundthresh,j) = norm_coeff;
        
        % write files
        for k = 1:nPCs
            
            if k <= nPCs
                name = component_names{k};
            else
                name = 'grandmean';
            end
            
            for i = 1:2
                
                fname = [surf_directory idstring 'beta-soundthresh' num2str(soundthresh) '-' name '-' hemis{i} '-r' sprintf('%d',runs) '.mgz'];
                x = ica_weights_soundthresh(k,(1:nsurfpts) + (i-1)*nsurfpts,j);
                xi = x~=0;
                xn = x;
                bounds = [-2 2];
                xn(x<bounds(1)*weight_std(k)) = 1;
                xn(x>bounds(2)*weight_std(k)) = 2;
                xn(x>bounds(1)*weight_std(k) & x<bounds(2)*weight_std(k) & xi) = interp1(weight_std(k)*bounds,[1 2], x(x>bounds(1)*weight_std(k) & x<bounds(2)*weight_std(k) & xi));
                MRIwrite_SNH(xn, fname, hemis{i});
                
                % save screenshot
                surface_image = [surface_projections_directory idstring 'beta-soundthresh'  num2str(soundthresh) '-' name '-us' num2str(usubs(j)) '-' hemis{i} '.png'];
                if false && ~exist(surface_image, 'file')
                    %                 freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.00001 0.00002 weight_std(k)*3],'piecewise');
                    %                 %                 freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[0.00001 0.00002 weight_std(k)*3],'piecewise');
                    %                 freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[0.25 1.25 2]);
                    freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.25 1.25 2]);
                end
                
                % write beta weights thresholded by sound responsivity
                fname = [surf_directory idstring 'logP-' name '-' hemis{i} '-r' sprintf('%d',runs) '.mgz'];
                MRIwrite_SNH(ica_logP(k,(1:nsurfpts) + (i-1)*nsurfpts,j), fname, hemis{i});
                
                % save screenshot
                surface_image = [surface_projections_directory idstring 'logP-' name '-us' num2str(usubs(j)) '-' hemis{i}  '.png'];
                if false && ~exist(surface_image, 'file')
                    freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[3 4.5 6]);
                    %                 freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[3 4.5 6]);
                end
            end
        end
    end
    save(whole_brain_projection_file, 'sound_responseP', 'ica_weights_soundthresh', 'ica_logP', 'normalization_coefficients')
else
    load(whole_brain_projection_file);
end

% surface directory
surf_directory = [params('rootdir') 'freesurfer/fsaverage/linear_analysis_surf/naturalsound_us' sprintf('-%d',usubs) '/'];
if ~exist(surf_directory,'dir')
    mkdir(surf_directory);
end

% average across subjects
x1 = ica_weights_soundthresh;
x1(x1==0) = NaN;
x2 = nanmean(x1,3);
x2(isnan(x2)) = 0;
ica_weights_group = x2;

% average across subjects
x1 = normalization_coefficients;
x1(x1==0) = NaN;
x2 = nanmean(x1,2);
x2(isnan(x2)) = 0;
normalization_coefficients_group = x2;

% normalized weights
ica_weights_group_norm = ica_weights_group./(ones(nPCs+1,1)*normalization_coefficients_group');

% x = ica_weights_soundthresh(i,:,:);
% x = x(sound_responseP > soundthresh);
% x = x(:)/std(x(:));
% hist(x,100);

% threshold by sound responsivity
group_voxelthresh = sum(sound_responseP > soundthresh,2) > min_subjects - 1e-3;
ica_weights_group_norm(:,~group_voxelthresh) = 0;
weight_std = std(ica_weights_group_norm(:,group_voxelthresh)');

% histogram
figure;
[y,bins] = hist(ica_weights_group_norm(1:nPCs,group_voxelthresh)',100);
plot(bins,y/sum(y(:,1)),'LineWidth',2);
legend(component_names{:});
export_fig([surface_projections_directory 'histogram.pdf'],'-pdf','-nocrop','-transparent');

% group results
for k = 2
    
    if k <= nPCs
        name = component_names{k};
    else
        name = 'grandmean';
    end
    
    for i = 1:1
        % group beta figure
        fname = [surf_directory idstring 'beta-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-' name '-' hemis{i} '.mgz'];
        delete(fname);
        x = ica_weights_group_norm(k,(1:nsurfpts) + (i-1)*nsurfpts);
        xi = x~=0;
        xn = x;
        bounds = [-2 2];
        xn(x<bounds(1)*weight_std(k)) = 1;
        xn(x>bounds(2)*weight_std(k)) = 2;
        xn(x>bounds(1)*weight_std(k) & x<bounds(2)*weight_std(k) & xi) = interp1(weight_std(k)*bounds,[1 2], x(x>bounds(1)*weight_std(k) & x<bounds(2)*weight_std(k) & xi));
        MRIwrite_SNH(xn, fname, hemis{i});
        
        %         % for plotting
        %         t = 3;
        %         x = ica_weights_group(k,:,i)/weight_std(k);
        %         x(x<-t) = -t;
        %         x(x>t) = t;
        %         ica_normalized_forplotting = interp1([-2,-0.4,1.2,2]*(t/2),[1 2.66, 4.33, 6], x);
        %         ica_normalized_forplotting(~group_voxelthresh) = 0;
        %
        %         % group beta figure
        %         fname = [surf_directory idstring 'beta-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-component' num2str(k) '-' hemis{i} '-forplotting.mgz'];
        %         delete(fname);
        %         MRIwrite_SNH(ica_normalized_forplotting, fname, hemis{i});
        
        % save screenshot
        surface_image = [surface_projections_directory idstring 'beta-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-' name '-group-' hemis{i} '.png'];
        if true || ~exist(surface_image, 'file')
            %             freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.00001 0.00002 weight_std(k)*2],'piecewise');
            %             freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[0.00001 0.00002 weight_std(k)*2],'piecewise');
            %             freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.25 1.25 2]);
            freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[1 1.5 2]);
        end
    end
end

addpath([pwd '/fs']);
clear cat;
ica_weights_augmented = cat(1, 0*ones(size(ica_weights_group(1,:))), ica_weights_group(1,:), -ica_weights_group(1,:), ica_weights_group([5 2 4],:));
[~,x] = max(ica_weights_augmented);
component_assignments = x-1;
component_assignments(~group_voxelthresh) = 0;

for i = 1:2
    
    % write labels
    lb = read_label(['~/freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
    labelstring = [];
    color_order = 1:5;
    for k = color_order
        
        % write label
        component_label = [surf_directory idstring 'label-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-component' num2str(k) '-' hemis{i} '.label'];
        
        %         hemis{i} '.' exp '_' idstring '_component' num2str(k-1) '_allsubjects_mask_threshold' num2str(round(mask_threshold*100)) '_component_threshold' num2str(round(component_threshold*100)) '.label'];
        if true || ~exist(component_label, 'file') || optInputs(varargin,'overwrite')
            vnums = find(component_assignments==k)-1 - (i-1)*nsurfpts;
            [~,inds] = intersect(lb.vnums, vnums);
            if isempty(inds)
                inds = 1;
            end
            lb_new = lb;
            lb_new.vnums = lb.vnums(inds);
            lb_new.vras = lb.vras(inds,:);
            
            delete(component_label);
            write_label(lb_new,component_label,varargin{:});
        end
        
        % save screenshot of label
        %         surface_image = [figure_directory 'cluster' num2str(j) '_k' num2str(cluster_order(j))  '_label_allsubjects_minoverlap' num2str(round(min_overlap*100)) '_' hemis{i} '.png'];
        %         if ~exist(surface_image, 'file') || optInputs(varargin, 'overwrite')
        %             freeview3('fsaverage', hemis{i}, 'label', cluster_label,'screenshot',surface_image);
        %         end
        labelstring = [labelstring ' --l ' component_label]; %#ok<AGROW>
        
    end
    
    ctab = [pwd '/DimensionsColorLUT2.txt'];
    annot_file = [surf_directory idstring 'bestcomponent-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-' hemis{i} '-r' sprintf('%d',runs) '.annot'];
    %     annot_file = [fsaverage_surfdir hemis{i} '.' exp '_' idstring '_allcomponents_allsubjects_mask_threshold' num2str(round(mask_threshold*100)) '_component_threshold' num2str(round(component_threshold*100)) '.annot'];
    if true || ~exist(annot_file, 'file') || optInputs(varargin,'overwrite')
        delete(annot_file);
        unix_freesurfer_version(freesurfer_version, ['mris_label2annot --s fsaverage --h ' hemis{i} ' --ctab ' ctab ' --annot-path ' annot_file ' ' labelstring]);
    end
    
    %     if ~exist(annot_file_format,'file') || optInputs(varargin,'overwrite')
    annot_file_format = strrep(annot_file, '.annot', '_format.annot');
    delete(annot_file_format);
    [vertices, label, ct] = read_annotation(annot_file);
    ct.table(1,1:3) = 1*ones(1,3);
    write_annotation(annot_file_format, vertices, label, ct);
    
    surface_image = [surface_projections_directory 'bestcomponent-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-' hemis{i} '.png'];
    %     freeview3('fsaverage', hemis{i}, 'annot', annot_file_format,'screenshot',surface_image);
    freeview3('fsaverage', hemis{i}, 'annot', annot_file_format,'screenshot',surface_image);
    freeview3('fsaverage', hemis{i}, 'annot', annot_file_format)
end

%% Follow-up experiments


% hemispheres
hemis = {'rh','lh'};
usubs = [1 11 106 121 122 123 143 160 161 164];
soundthresh = 3;
    
pitch_tonotopy_response_allsubjects = nan(length(usubs), 8, 6);
music_speech_allsubjects = nan(length(usubs), 8, 6);

% Une_pitch_invert = Une;
% Une_pitch_invert(:,ismember(component_names,'pitch')) = -Une(:,ismember(component_names,'pitch'));
for j = 1:length(usubs)
    
    %% Natural sound localizer
    
    % subject id
    fprintf('us %d\n', usubs(j)); drawnow;
    subjid = ['naturalsound_us' num2str(usubs(j))];
    
    % surface directory
    surf_directory = [params('rootdir') 'freesurfer/fsaverage/linear_analysis_surf/' subjid '/'];
    if ~exist(surf_directory,'dir')
        mkdir(surf_directory);
    end
    
    % runs for that subject
    runs = read_runs('naturalsound',usubs(j),'main_v3_combined');
    
    % hemispheres
    sv = zeros(165, nsurfpts*2);
    for i = 1:2
        
        % auditory cortex mask
        mask = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
        
        % voxel responses
        for q = 1:length(runs)
            x = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r' num2str(runs(q)) '/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp3-6.mgz']);
            dat = squeeze(x.vol)';
            sv(:,mask.vnums+1 + (i-1)*nsurfpts) = sv(:,mask.vnums+1 + (i-1)*nsurfpts) + dat(:,mask.vnums+1)/length(runs);
        end
    end
        
    % sound response
    %     [~,p] = ttest(sv,0,'tail','right');
    [~,p] = ttest(sv,0,0.05,'right');
    sound_responseP = -log10(p);
    
    % select voxels with a significant sound response
    sv_selected = sv(:,sound_responseP > soundthresh);
    
    % Regression analysis
    Uica_norm = Une ./ (ones(165,1)*std(Une));
    xhat = pinv([Uica_norm, grandmean/rms(grandmean)])*sv_selected;
    Uest = sv_selected*pinv(xhat);
    
    min_response = min(Uest);
    max_response = max(Uest);
    response_range = max_response-min_response;
    
    %% Texture
    test_runs = read_runs('naturalsound',usubs(j), 'texture_combined');
    if ~isempty(test_runs)
        test_conds = read_conditions('naturalsound',usubs(j),'texture_combined',varargin{:});
        test_response = nan(length(test_runs), length(test_conds), nPCs+1);
        for q = 1:length(test_runs)
            test_data = nan(length(test_conds), nsurfpts*2);
            for i = 1:2
                % auditory cortex mask
                mask = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
                
                % voxel responses
                fname = [params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/texture_combined_r' num2str(test_runs(q)) '/' hemis{i} '.sigav_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp3-6.mgz'];
                x = MRIread(fname);
                vol = squeeze(x.vol)';
                test_data(:,mask.vnums+1 + (i-1)*nsurfpts) = vol(:,mask.vnums+1);
            end
            test_data_selected = test_data(:,sound_responseP > soundthresh);
            xi = all(~isnan(test_data_selected));
            x = test_data_selected(:,xi)*pinv(xhat(:,xi));
            test_response(q,:,:) = (x - ones(size(x,1),1)*min_response) ./ (ones(size(x,1),1)*response_range);
        end
        
        test_response_mean = squeeze(mean(test_response,1));
        stim_by_condition_by_component = reshape(test_response_mean,[36, 5, nPCs+1]);        
        
        % select subset of stimuli
        [~,C] = naturalsound_acoustic_category_regressors;
        load('selected_stimuli_36stims.mat');
        xi = find(ismember(C.ids,ids));
        UneSubset = [Une(xi,:), grandmean(xi)]; 
        
        % correlation
        condition_correlation_by_component = nan(5,nPCs+1);
        condition_stderr_by_component = nan(2,5,nPCs+1);
        for q = 1:nPCs+1
            for z = 1:5
                x = stim_by_condition_by_component(:,z,q);
                [condition_correlation_by_component(z,q),~,~,~,condition_stderr_by_component(:,z,q)] = bootstrap_correlation(UneSubset(:,q),x);
            end
        end
        
        condition_names = {'orig-1','orig-2','marg','modpower','fullmodel'};
        % plot correlation
        for q = 1:nPCs+1
            lerr = condition_stderr_by_component(1,:,q)';
            uerr = condition_stderr_by_component(2,:,q)';
            mybar(condition_correlation_by_component(:,q).^2, condition_names, condition_names([3,4,5,1,2]), 'errorbar2', lerr.^2, uerr.^2);
            if q <= nPCs
                name = component_names{q};
            else
                name = 'grandmean';
            end
            title(name);
            box off;
            if ~exist([figure_directory 'texture/us' num2str(usubs(j))],'dir');
                mkdir([figure_directory 'texture/us' num2str(usubs(j))]);
            end
            ylabel('Correlation (R^2)');
            fname = [figure_directory 'texture/us' num2str(usubs(j)) '/correlation_' name '_us' num2str(usubs(j)) '.pdf'];
            if exist(fname,'file');
                delete(fname);
            end
            export_fig(fname,'-pdf','-nocrop','-transparent');
        end
        
        % mean response amplitude
        for q = 1:nPCs+1
            if q <= nPCs
                name = component_names{q};
            else
                name = 'grandmean';
            end
            x = stim_by_condition_by_component(:,:,q);
            condition_names = {'orig-1','orig-2','marg','modpower','fullmodel'};
            mybar(mean(x), condition_names, condition_names([3,4,5,1,2]), 'errorbar', stderr_withsub_corrected(x));
            title(name);
            fname = [figure_directory 'texture/us' num2str(usubs(j)) '/response_amplitude_' name '_us' num2str(usubs(j)) '.pdf'];
            if exist(fname,'file');
                delete(fname);
            end
            export_fig(fname,'-pdf','-nocrop','-transparent');

        end
        
        % response pattern
        for q = 1:nPCs+1
            if q <= nPCs
                name = component_names{q};
            else
                name = 'grandmean';
            end
            condition_names = {'orig-1','orig-2','marg','modpower','fullmodel'};
            for z = 1:5
                y = stim_by_condition_by_component(:,z,q);
                response_profile_figure(y, Csubset, stim_names,[-0.2 1.2]);
                export_fig([figure_directory 'texture/us' num2str(usubs(j)) '/response_profile_' name '_' condition_names{z} '.eps'],'-eps','-nocrop','-transparent');
            end
        end
    end
    close all;
    
    %% scrambling experiment
    if any(usubs(j) == [1 11 106 121 123 143 160 164])
        if usubs(j) == 143
            runtype = 'scrambling_russian';
        else
            runtype = 'scrambling';
        end
        if usubs(j) == 164
            test_runs = 2:3;
        else
            test_runs = read_runs('naturalsound',usubs(j), runtype);
        end
        test_conds = read_conditions('naturalsound',usubs(j),runtype,varargin{:});
        test_response = nan(length(test_runs), length(test_conds), nPCs+1);
        
        for q = 1:length(test_runs)
            test_data = nan(length(test_conds), nsurfpts*2);
            for i = 1:2
                % auditory cortex mask
                mask = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
                
                % voxel responses
                fname = [params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/' runtype '_r' num2str(test_runs(q)) '/' hemis{i} '.sigav_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp3-6.mgz'];
                x = MRIread(fname);
                vol = squeeze(x.vol)';
                test_data(:,mask.vnums+1 + (i-1)*nsurfpts) = vol(:,mask.vnums+1);
            end
            test_data_selected = test_data(:,sound_responseP > soundthresh);
            xi = all(~isnan(test_data_selected));
            x = test_data_selected(:,xi)*pinv(xhat(:,xi));
            test_response(q,:,:) = (x - ones(size(x,1),1)*min_response) ./ (ones(size(x,1),1)*response_range);
        end
        
        test_response_mean = squeeze(mean(test_response));
        test_response_stderr = stderr_withsub_corrected(test_response);
        for q = 1:nPCs+1
            if q <= nPCs
                name = component_names{q};
            else
                name = 'grandmean';
            end
            if ~exist([figure_directory 'scrambling/'],'dir');
                mkdir([figure_directory 'scrambling/']);
            end
            mybar(test_response_mean(:,q), test_conds, test_conds, 'errorbar', test_response_stderr(:,q),'ylim',[-0.2 1.2]);
            title(name);
            box off;
            export_fig([figure_directory 'scrambling/' name '_us' num2str(usubs(j)) '.pdf'],'-pdf','-nocrop','-transparent');
        end
        music_speech_allsubjects(j,:,:) = test_response_mean;
    end
    close all;
    
    %% pitch and tonotopy localizer
    if usubs(j) == 11
        test_runs = read_runs('pitch_resthr_v4',usubs(j),'localizer','first_Nruns',6);
        test_conds = read_conditions('pitch_resthr_v4',usubs(j),'localizer',varargin{:});
    else
        test_runs = read_runs('naturalsound',usubs(j),'localizer','first_Nruns',6);
        test_conds = read_conditions('naturalsound',usubs(j),'localizer',varargin{:});
    end
    test_response = nan(length(test_runs), length(test_conds), nPCs+1);
    for q = 1:length(test_runs)
        test_data = nan(length(test_conds), nsurfpts*2);
        for i = 1:2
            % auditory cortex mask
            mask = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
            
            % voxel responses
            if usubs(j) == 11
                fname = [params('rootdir') 'freesurfer/fsaverage/fla/pitch_resthr_v4_us' num2str(usubs(j)) '/localizer_r' num2str(test_runs(q)) '/' hemis{i} '.sigav_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp3-6.mgz'];
            else
                fname = [params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/localizer_r' num2str(test_runs(q)) '/' hemis{i} '.sigav_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp3-6.mgz'];
            end
            x = MRIread(fname);
            vol = squeeze(x.vol)';
            test_data(:,mask.vnums+1 + (i-1)*nsurfpts) = vol(:,mask.vnums+1);
        end
        test_data_selected = test_data(:,sound_responseP > soundthresh);
        xi = all(~isnan(test_data_selected));
        x = test_data_selected(:,xi)*pinv(xhat(:,xi));
        test_response(q,:,:) = (x - ones(size(x,1),1)*min_response) ./ (ones(size(x,1),1)*response_range);
    end
    
    if ~exist([figure_directory 'pitch_tonotopy_localizer/'],'dir');
        mkdir([figure_directory 'pitch_tonotopy_localizer/']);
    end
    pitch_tonotopy_response_mean = squeeze(mean(test_response));
    pitch_tonotopy_response_stderr = stderr_withsub_corrected(test_response);
    for q = 1:size(xhat,1)
        mybar(pitch_tonotopy_response_mean(:,q), test_conds, test_conds, 'errorbar', pitch_tonotopy_response_stderr(:,q),'ylim',[-0.2 1.2]);
        title(component_names{q});
        box off;
        
        export_fig([figure_directory 'pitch_tonotopy_localizer/' component_names{q} '_us' num2str(usubs(j)) '.pdf'],'-pdf','-nocrop','-transparent');
    end
    
    pitch_tonotopy_response_allsubjects(j,:,:) = pitch_tonotopy_response_mean;
    close all;

    %% localizer run just in me awhile ago
    %     if usubs(j) == 1
    %         test_runs = 1:6;
    %         conds = read_conditions('naturalsound',usubs(j),'natsoundloc',varargin{:});
    %         test_response = nan(length(test_runs), length(conds), size(xhat,1));
    %         for q = 1:length(test_runs)
    %             dat_hemi = cell(1,2);
    %             for i = 1:length(hemis)
    %                 fname = [params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/natsoundloc_r' num2str(test_runs(q)) '/' hemis{i} '.sigav_' num2str(100*fwhm, '%.0f') 'mm_hpfilt-Infsec-order2_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.mgz'];
    %                 x = MRIread(fname);
    %                 dat_hemi{i} = squeeze(x.vol)';
    %             end
    %             dat = cat(2,dat_hemi{1},dat_hemi{2});
    %             dat_selected = dat(:,sound_responseP > soundthresh);
    %             dat_selected_demeaned = demean_subjects(dat_selected, ones(1,size(dat_selected,2)));
    %             x = dat_selected_demeaned*pinv(xhat);
    %             test_response(q,:,:) = (x - ones(size(x,1),1)*min_response) ./ (ones(size(x,1),1)*response_range);
    %         end
    %
    %         test_response_mean = squeeze(mean(test_response));
    %         test_response_stderr = stderr_withsub_corrected(test_response);
    %         for q = 1:size(xhat,1)
    %             %             if strcmp(component_names{q},'pitch');
    %             %                 mybar(-test_response_mean(:,q), conds, conds, 'errorbar', test_response_stderr(:,q));
    %             %             else
    %             mybar(test_response_mean(:,q), conds, conds, 'errorbar', test_response_stderr(:,q),'ylim',[-0.2 1.2]);
    %             %             end
    %             title(component_names{q});
    %             box off;
    %             export_fig([figure_directory 'music_speech_localizer/' component_names{q} '_us' num2str(usubs(j)) '.pdf'],'-pdf','-nocrop','-transparent');
    %         end
    %     end
    %         [b, ~, logP, ~] = regress_SNH(Uica_norm, sv_norm_demeaned);
    %         [b, ~, logP, ~] = regress_SNH([Uica_norm, grandmean/std(grandmean)], sv_selected);
end

group_pitch_tonotopy_mean = squeeze(mean(pitch_tonotopy_response_allsubjects));
group_pitch_tonotopy_stderr = stderr_withsub_corrected(pitch_tonotopy_response_allsubjects);
for q = 1:size(xhat,1)
    %     if strcmp(component_names{q},'pitch');
    %         mybar(-group_pitch_tonotopy_mean(:,q), conds, conds, 'errorbar', group_pitch_tonotopy_stderr(:,q));
    %     else
    mybar(group_pitch_tonotopy_mean(:,q), conds, conds, 'errorbar', group_pitch_tonotopy_stderr(:,q),'ylim',[-0.2 1.2]);
    %     end
    title(component_names{q});
    box off;
    export_fig([figure_directory 'pitch_tonotopy_localizer/group_' component_names{q} '.pdf'],'-pdf','-nocrop','-transparent');
end

conds = read_conditions('naturalsound',usubs(j),'scrambling',varargin{:});
group_scrambling_mean = squeeze(mean(music_speech_allsubjects(ismember(usubs, usubs_three_scans),:,:)));
group_scrambling_stderr = stderr_withsub_corrected(music_speech_allsubjects(ismember(usubs, usubs_three_scans),:,:));
for q = 1:size(xhat,1)
    %     if strcmp(component_names{q},'pitch');
    %         mybar(-group_pitch_tonotopy_mean(:,q), conds, conds, 'errorbar', group_pitch_tonotopy_stderr(:,q));
    %     else
    mybar(group_scrambling_mean(:,q), conds, conds, 'errorbar', group_scrambling_stderr(:,q),'ylim',[-0.2 1.2]);
    %     end
    title(component_names{q});
    box off;
    export_fig([figure_directory 'scrambling/group_' component_names{q} '.pdf'],'-pdf','-nocrop','-transparent');
end

% component_name = ['linear_analysis_surf_naturalsound_speech'];
% roi_pscAll(1, 'naturalsound', 'main_v3_combined',  'sigav_plateau',5,'withsub_sem');

%% ICA analyses measuring reliability in left-out data

nPCs = 6;
nreps = 1000;
seed = 1;
NegEntOpt_matfile = [figure_directory 'NegEntOpt-Scan12-' num2str(nPCs) '-' num2str(nreps) '_seed' num2str(seed) '.mat'];
if ~exist(NegEntOpt_matfile,'file')
    [U,S,V] = svd(v12(:,ismember(si, usubs_three_scans)),'econ');
    VrScan12 = V(:,1:nPCs)'*sqrt(size(V,1)-1);
    UrScan12 = U(:,1:nPCs)/sqrt(size(VrScan12,1)-1);
    SrScan12 = S(1:nPCs,1:nPCs);
    [Vopt_best, h_pairs_best, negEnt_best, Vopt_all, h_pairs_all, negEnt_all] = ICA_NegEnt(VrScan12, nreps, seed, 1);
    save(NegEntOpt_matfile, 'Vopt_best', 'h_pairs_best', 'negEnt_best', 'Vopt_all', 'h_pairs_all', 'negEnt_all');
else
    [U,S,V] = svd(v12(:,ismember(si, usubs_three_scans)),'econ');
    VrScan12 = V(:,1:nPCs)'*sqrt(size(V,1)-1);
    UrScan12 = U(:,1:nPCs)/sqrt(size(V,1)-1);
    SrScan12 = S(1:nPCs,1:nPCs);
    load(NegEntOpt_matfile);
end

% scan12_component_names = {'pitch','speech','music','tonotopy','envsounds'};

% reformat matrices
VneScan12 = Vopt_best;
UneScan12 = UrScan12*SrScan12*VrScan12*pinv(VneScan12);
xi = sign(skewness(VneScan12'))==-1;
UneScan12(:,xi) = -UneScan12(:,xi);
VneScan12(xi,:) = -VneScan12(xi,:);

% 
% UneScan12 = UneScan12(:,[1 2 3 4 6]);
% VneScan12 = VneScan12([1 2 3 4 6],:);
% UneScan12(:,1) = -UneScan12(:,1);
% VneScan12(1,:) = -VneScan12(1,:);
% nPCs = 5;

Vica_leftout = pinv(UneScan12)*v3(:,ismember(si, usubs_three_scans));
for i = 1:size(Vica_leftout,1)
    Vica_leftout(i,:) = Vica_leftout(i,:)/std(Vica_leftout(i,:));
end

Vpca_leftout = pinv(UrScan12)*v3(:,ismember(si, usubs_three_scans));
for i = 1:size(Vpca_leftout,1)
    Vpca_leftout(i,:) = Vpca_leftout(i,:)/std(Vpca_leftout(i,:));
end

plot_rotation_negentropy_v2(VneScan12,[0 0.4]);
export_fig([figure_directory 'NegEnt-NotLeftOut-ICA-' num2str(nPCs) '.pdf'],'-nocrop','-transparent','-pdf');

plot_rotation_negentropy_v2(Vica_leftout,[0 0.4]);
export_fig([figure_directory 'NegEnt-LeftOut-ICA-' num2str(nPCs) '.pdf'],'-nocrop','-transparent','-pdf');

plot_rotation_negentropy_v2(Vpca_leftout(1:5,:),[0 0.4]);
export_fig([figure_directory 'NegEnt-LeftOut-PCA-' num2str(nPCs) '.pdf'],'-nocrop','-transparent','-pdf');

naturalsound_acoustic_category_figures(UneScan12, figure_directory, ['ICA-' num2str(nPCs) '-Scan12-'],'corr_bounds',[-1 1]);

%% Reliability of components

Vne_threescans = VneScan12;

% correlation of random projections
x = rand(size(v12_threescans,2),10000);
randprojection1 = v12_threescans*x;
randprojection2 = v3_threescans*x;
r_randprojection = fastcorr3(randprojection1,randprojection2);

% reliability of ICA components
UneScan1 = v1_threescans*pinv(Vne_threescans);
UneScan2 = v2_threescans*pinv(Vne_threescans);
UneScan3 = v3_threescans*pinv(Vne_threescans);
r_ica_decomposition = fastcorr3(UneScan12,UneScan3);
% names = {'speech', 'pitch', 'tonotopy', 'envsounds', 'music'};

% split-half correlation
for i = 1:size(UneScan12,2)
    [~, C] = naturalsound_acoustic_category_regressors;
    figure;
    x = [zscore(UneScan12(:,i)), zscore(UneScan3(:,i))];
    x = reshape(norm01(x),[165 2]);
    colorhandles = nan(1,length(C.category_labels));
    for j = 1:165
        ct = C.category_assignments(j);
        h = plot(x(j,1),x(j,2),'o','Color',C.colors(ct,:),'LineWidth',3);
        if j == 1;
            hold on;
        end
        colorhandles(ct) = h;
    end
    plot(0:0.1:1, 0:0.1:1, 'k', 'LineWidth',2);
    xlim([-0.05 1.05]); ylim([-0.05 1.05])
    title(sprintf('r = %.2f', corr(x(:,1),x(:,2))));
    set(gca,'FontSize',14,'FontWeight','Demi');
    box off;
    export_fig([figure_directory 'ICA-splithalf-scatter-' scan12_component_names{i} '.pdf'],'-pdf','-transparent');
end

% split-half correlation compared with random projections
figure;
[y,x] = hist(r_randprojection(:).^2,100);
plot(x,y/sum(y),'k','LineWidth',2);
hold on;
plot(median(r_randprojection(:).^2)*[1 1], [0 max(y/sum(y))],'k--','LineWidth',2);
cols = colormap('lines(5)');
hold on;
for j = 1:nPCs
    plot(r_ica_decomposition(j).^2*[1 1], [0 max(y/sum(y))],'Color',cols(j,:),'LineWidth',2);
end
title('Test-Retest Reliability');
xlabel('R^2');
xlim([0 1]);
legend([{'random projections','median projection'},scan12_component_names],'Location','NorthWest');
box off;
export_fig([figure_directory 'test_retest_reliability.pdf'],'-pdf','-nocrop','-transparent');

% raw voxel reliability
figure;
r_allvoxels = fastcorr([v1_threescans,v2_threescans],[v3_threescans,v3_threescans]);
[y,x] = hist(r_allvoxels(:),100);
plot(x,y/sum(y),'k','LineWidth',2);
hold on;
plot(median(r_allvoxels(:))*[1 1], [0 max(y/sum(y))],'k--','LineWidth',2);
xlabel('R');
title('Voxel Reliability for Single Scan');
box off;
export_fig([figure_directory 'voxel_reliability_single_scan.pdf'],'-pdf','-nocrop','-transparent');

% 

%% Example voxels

ResetRandStream(1)
xi = Shuffle(1:size(v,2));
v_raw = v_unfilt_mean(:,voxel_selection);
for i = 1:10
    figure;
    response = v_raw(:,xi(i));
    bar(1:165,100*response, 'FaceColor',0.25*[1 1 1]);
    ylim([-0.5 3.5]); xlim([0 166]);
    naturalsound_acoustic_category_figures(response, figure_directory, ['example_voxel' num2str(xi(1)) '-']);
    export_fig([figure_directory 'example_voxel' num2str(xi(i)) '.pdf'],'-pdf','-transparent');
end

%% Example stimulus correlations

% stim_inds = [find(strcmp('stim465_song_from_a_musical',conds)) find(strcmp('stim505_chinese',conds))];
% stim_inds = [find(strcmp('stim268_piano',conds)) find(strcmp('stim505_chinese',conds))];
conds = read_conditions('naturalsound',1,'main_v3_combined');
stim_inds = [find(strcmp('stim232_man_speaking',conds)) find(strcmp('stim505_chinese',conds))];
xi = Shuffle(1:size(v,2));
scatter(100*v(stim_inds(1),xi), 100*v(stim_inds(2),xi),1);
xlabel('Man Speaking'); ylabel('Chinese Speech');
box off;
export_fig([figure_directory 'example_correlation.pdf'],'-pdf','-transparent');

%% Clustering

% reformat data for clustering analysis
v_wmean = v+grandmean*ones(1,size(v,2));
vox_response_matrix_bysubject = cell(length(usubs),1);
for i = 1:length(usubs)
    vox_response_matrix_bysubject{i} = v_wmean(:,usubs(i) == si)';
end

% clustering
k = 5;
n_repetitions = 20;
random_seed = 1;
ResetRandStream(random_seed);
[grpRes, indRes, cluster_idstring] = response_cluster_discover('naturalsound', usubs, vox_response_matrix_bysubject, v_wmean', si, k, n_repetitions, random_seed, figure_directory, 1);
[~,cluster_order] = sort(grpRes.p,'descend');
close all;
naturalsound_acoustic_category_figures(grpRes.m(cluster_order,:)', figure_directory, ['clustering-' num2str(k) '-']);

% cluster assignments
[~,cluster_order] = sort(grpRes.p,'descend');
cluster_assignments = nan(size(v,2),1);
for i =1:length(cluster_order)
    xi = grpRes.clusters == cluster_order(i);
    cluster_assignments(xi) = i;
end

% cluster colors
fid = fopen([pwd '/ClusterColorLUT2.txt'],'r');
x = textscan(fid, '%d%s%d%d%d%d', 'HeaderLines',4); fclose(fid);
colors = colormap('lines');
colors_clusters = colors(1:k,:);

% plot in PCA space
xi = Shuffle(1:size(v,2));
[U,S,V] = svd(v,'econ');
for i = 1:2:5
    figure;
    colormap(colors_clusters);
    scatter(V(xi,i),V(xi,i+1),1,cluster_assignments(xi));
    h = colorbar; set(h,'YTick',1:k);
    xlabel(['PCA' num2str(i)]);
    ylabel(['PCA' num2str(i+1)]);
    box off;
    export_fig([figure_directory 'pca-components' num2str(i) '-' num2str(i+1) '_clusters-k' num2str(k) '.pdf'], '-pdf');
end

% response matrix
clusResponseMatrix = nan(165,k);
for i = 1:k
    clusResponseMatrix(:,i) = mean(v(:,cluster_assignments==i),2);
end

%% Behavioral data

conds = read_conditions('naturalsound',1,'main_v3_combined');
hits_allsubjects = nan(length(conds),length(usubs));
for i = 1:length(usubs)
    [~,~,allscans] = read_runs('naturalsound',usubs(i),'main_v3_combined');
    hits = nan(length(conds), length(allscans));
    for j = 1:length(allscans)
        allruns = read_runs('naturalsound',usubs(i),'main_v3','scans',allscans(j));
        for k = 1:length(allruns)
            b = read_behav('naturalsound',usubs(i),'main_v3',allruns(k));
            xi = find(b.targ==1);
            for q = 1:length(xi)
                x = b.conds{xi(q)};
                if strcmp(x,'stim283_contemporary_r&b');
                    x = 'stim283_contemporary_rb';
                end
                condition_index = strcmp(x, conds);
                if sum(condition_index)~=1
                    error('Bad condition index');
                end
                hits(condition_index, j) = ~isnan(b.rkey(xi(q)));
            end
        end
    end
    hits_allsubjects(:,i) = mean(hits,2);
end
mean_hitrate = mean(hits_allsubjects,2);
r = corr(Une, mean_hitrate);
mybar(r', component_names, component_names);
ylabel('Hit rate correlation (r)');
box off;
export_fig([figure_directory 'ICA_hit-rate-correlation.pdf'],'-pdf','-transparent');

%% Cochleogram plots

conds = read_conditions('naturalsound',1,'main_v3_combined');
for i = 1:size(Une,2)
    [~,xi] = sort(Une(:,i), 'descend');
    ordered_stimuli = conds(xi);
    pdf_list = cell(size(ordered_stimuli));
    for j = 1:length(ordered_stimuli);
        pdf_list{j} = [params('rootdir') 'naturalsound/figures/cochleograms/' ordered_stimuli{j} '.pdf'];
        if j <= 5 || j >= 161
            copyfile(pdf_list{j},[figure_directory 'cochleograms_' component_names{i} '_stim' num2str(j) '.pdf'])
        end
    end
    append_pdfs([figure_directory 'cochleograms_' component_names{i} '.pdf'],pdf_list{:})
end

%% Ordered stimuli

conds = read_conditions('naturalsound',1,'main_v3_combined');
for i = 1:size(Une,2)
    [~,xi] = sort(Une(:,i), 'descend');
    ordered_stimuli = conds(xi);
    sr = 44100;
    stim = zeros(round(165*(2+0.2)*sr),1);
    for j = 1:length(ordered_stimuli);
        wav = wavread([params('rootdir') 'naturalsound/stims/final_stims165/' ordered_stimuli{j} '.wav']);
        stim((1:length(wav)) + round((j-1)*(2+0.2)*sr)) = wav;
    end
    wavwrite(stim, sr, 16, [figure_directory 'ordered_stimuli_' component_names{i} '.wav'])
end

%% PCA for subsets of subjects

reliability = nan(1,length(usubs));
for i = 1:length(usubs)
    %     reliability(i) = mean(corrP(voxel_selection & si_unfilt' == usubs(i)));
    reliability(i) = mean(corrP(~isnan(corrP) & si_unfilt' == usubs(i)));
end
[~,xi] = sort(reliability,'descend');

nvoxels_persubject = nan(1,length(usubs));
for i = 1:length(usubs)
    nvoxels_persubject(i) = sum(si == usubs(i));
end
[~,xi] = sort(nvoxels_persubject,'descend');
usubs_sg1 = usubs(xi(2:2:end));
usubs_sg2 = usubs(xi(1:2:end));

si_sg1 = si(ismember(si,usubs_sg1));
v_sg1 = v(:,ismember(si,usubs_sg1));
v1_sg1 = v1(:,ismember(si,usubs_sg1));
v2_sg1 = v2(:,ismember(si,usubs_sg1));

sg1_directory = [figure_directory 'sg1/'];
if ~exist(sg1_directory, 'dir');
    mkdir(sg1_directory);
end
pca_reliability(v_sg1, v1_sg1, v2_sg1, usubs_sg1, si_sg1, sg1_directory);

si_sg2 = si(ismember(si,usubs_sg2));
v_sg2 = v(:,ismember(si,usubs_sg2));
v1_sg2 = v1(:,ismember(si,usubs_sg2));
v2_sg2 = v2(:,ismember(si,usubs_sg2));

sg2_directory = [figure_directory 'sg2/'];
if ~exist(sg2_directory, 'dir');
    mkdir(sg2_directory);
end
pca_reliability(v_sg2, v1_sg2, v2_sg2, usubs_sg2, si_sg2, sg2_directory);

usubs_best = [1 164 11 143 121];
si_sg_best = si(ismember(si,usubs_best));
v_sg_best = v(:,ismember(si,usubs_best));
v1_sg_best = v1(:,ismember(si,usubs_best));
v2_sg_best = v2(:,ismember(si,usubs_best));
v3_sg_best = v3(:,ismember(si,usubs_best));
v12_sg_best = v12(:,ismember(si,usubs_best));

sg_best_directory = [figure_directory 'sg_best5/'];
if ~exist(sg_best_directory, 'dir');
    mkdir(sg_best_directory);
end
pca_reliability(v_sg_best, v1_sg_best, v2_sg_best, usubs_best, si_sg_best, sg_best_directory);

pca_reliability_mediancorr(v_sg_best, v12_sg_best, v3_sg_best, usubs_best, si_sg_best, [sg_best_directory 'scans123-']);
pca_reliability_mediancorr(v12_sg_best, v12_sg_best, v3_sg_best, usubs_best, si_sg_best, [sg_best_directory 'scans12-']);

pca_reliability_mediancorr(v, v12_sg_best, v3_sg_best, usubs_best, si_sg_best, [figure_directory '5b-scan123-' ]);
pca_reliability_mediancorr(v12_sg_best, v12_sg_best, v3_sg_best, usubs_best, si_sg_best, [figure_directory '5b-scan12-' ]);

component_subject_specificity(v_sg_best, si_sg_best, usubs_best, [figure_directory '5b-scan123-']);

%% ICA for subsets of subjects

nPCs = 5;
[U1,S1,V1] = svd(v12_sg_best,'econ');
V1r = V1(:,1:nPCs)'*sqrt(size(V1,1)-1);
U1r = U1(:,1:nPCs)/sqrt(size(V1,1)-1);
S1r = S1(1:nPCs,1:nPCs);
[Vne1] = ICA_NegEnt(V1r, 10, 1, 1);

Une1 = U1r*S1r*V1r*pinv(Vne1);
xi = sign(skewness(Vne1'))==-1;
Une1(:,xi) = -Une1(:,xi);
Vne1(xi,:) = -Vne1(xi,:);
UneNorm1 = nan(size(Une1));
for i = 1:nPCs
    UneNorm1(:,i) = (Une1(:,i) - min(Une1(:,i)))/(max(Une1(:,i))-min(Une1(:,i)));
end

% plot figures
naturalsound_acoustic_category_figures(UneNorm1, sg1_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-Norm-'],'corr_bounds',[-1 1]);

% plot negentropy figures
plot_rotation_negentropy(Vne1,[0 0.4]);
export_fig([sg1_directory 'ICA-' num2str(nPCs) '-NegEntOpt-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');

% plot histogram by subject
plothist_bysubject(Vne1, si_sg1, [sg1_directory 'ICA-' num2str(nPCs) '-NegEntOpt-hist-bysubject']);

nPCs = 5;
[U2,S2,V2] = svd(v_sg2,'econ');
V2r = V2(:,1:nPCs)'*sqrt(size(V2,1)-1);
U2r = U2(:,1:nPCs)/sqrt(size(V2,1)-1);
S2r = S2(1:nPCs,1:nPCs);
[Vne2] = ICA_NegEnt(V2r, 100, 1, 1);

Une2 = U2r*S2r*V2r*pinv(Vne2);
xi = sign(skewness(Vne2'))==-1;
Une2(:,xi) = -Une2(:,xi);
Vne2(xi,:) = -Vne2(xi,:);
UneNorm2 = nan(size(Une2));
for i = 1:nPCs
    UneNorm2(:,i) = (Une2(:,i) - min(Une2(:,i)))/(max(Une2(:,i))-min(Une2(:,i)));
end

% plot figures
naturalsound_acoustic_category_figures(UneNorm2, sg2_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-Norm-'],'corr_bounds',[-1 1]);

% plot negentropy figures
plot_rotation_negentropy(Vne2,[0 0.4]);
export_fig([sg2_directory 'ICA-' num2str(nPCs) '-NegEntOpt-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');

% best subjects in terms of test-retest reliability
nPCs = 5;
[Ubest,Sbest,Vbest] = svd(v12_sg_best,'econ');
VbestNorm = Vbest'*sqrt(size(Vbest,1)-1);
Vbestr = Vbest(:,1:nPCs)'*sqrt(size(Vbest,1)-1);
Ubestr = Ubest(:,1:nPCs)/sqrt(size(Vbest,1)-1);
Sbestr = Sbest(1:nPCs,1:nPCs);
[Vnebest] = ICA_NegEnt(Vbestr, 10, 1, 1);

% PCA histogram
plothist_bysubject(VbestNorm(1:10,:), si_sg_best, [sg_best_directory 'PCA-10-hist-bysubject']);

% PCA statistics
% plothist_bysubject(Vnebest, si_sg1, [sg1_directory 'ICA-' num2str(nPCs) '-NegEntOpt-hist-bysubject']);

Unebest = Ubestr*Sbestr*Vbestr*pinv(Vnebest);
xi = sign(skewness(Vnebest'))==-1;
Unebest(:,xi) = -Unebest(:,xi);
Vnebest(xi,:) = -Vnebest(xi,:);
UneNormbest = nan(size(Unebest));
for i = 1:nPCs
    UneNormbest(:,i) = (Unebest(:,i) - min(Unebest(:,i)))/(max(Unebest(:,i))-min(Unebest(:,i)));
end

naturalsound_acoustic_category_figures(UneNormbest, sg_best_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-Norm-'],'corr_bounds',[-1 1]);
UneNorm1
naturalsound_acoustic_category_figures(UneNormbest, sg_best_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-Norm-'],'corr_bounds',[-1 1]);

% plot figures
plothist_bysubject(Vnebest, si_sg_best, [sg_best_directory 'ICA-' num2str(nPCs) '-NegEntOpt-hist-bysubject']);

% plot negentropy figures
plot_rotation_negentropy(Vnebest,[0 0.4]);
export_fig([sg_best_directory 'ICA-' num2str(nPCs) '-NegEntOpt-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');

%% Individual PCA analysis

for i = 1:length(usubs)
    us = usubs(i);
    individual_subject_directory = [figure_directory 'us' num2str(us) '/'];
    if ~exist(individual_subject_directory, 'dir');
        mkdir(individual_subject_directory);
    end
    v1_us = v12(:,si==us);
    v2_us = v3(:,si==us);
    
    r = nan(1,165);
    scan1 = v1_us;
    scan2 = v2_us;
    [U_loo,~,~] = svd(scan1,'econ');
    xhat = U_loo'*scan1;
    for j = 1:165
        prediction = U_loo(:,1:j)*xhat(1:j,:);
        r(j) = fastcorr3(prediction(:), scan2(:));
    end
    
    % individual subject plot
    r2 = [r.^2];
    figure;
    plot(r2([1:20,165])','-o');
    xlim([0 22]);
    ylim([0 1]);
    set(gca, 'XTick',[1:21],'XTickLabel',[1:20,165])
    xlabel('Number of PCs'); ylabel('Explained Variance in Left-Out Data');
    box off;
    export_fig([individual_subject_directory 'PCA_loo.pdf'],'-pdf','-nocrop','-transparent');
end

%% Individual Subject ICA Analysis

nPCs = 5;
for i = 1:length(usubs)
    us = usubs(i);
    individual_subject_directory = [figure_directory 'us' num2str(us) '/'];
    if ~exist(individual_subject_directory, 'dir');
        mkdir(individual_subject_directory);
    end
    v_us = v(:,si==us);
    
    [Uus,Sus,Vus] = svd(v_us,'econ');
    Vusr = Vus(:,1:nPCs)'*sqrt(size(Vus,1)-1);
    Uusr = Uus(:,1:nPCs)/sqrt(size(Vus,1)-1);
    Susr = Sus(1:nPCs,1:nPCs);
    [Vne_us] = ICA_NegEnt(Vusr, 10, 1, 1);
    
    Une_us = Uusr*Susr*Vusr*pinv(Vne_us);
    xi = sign(skewness(Vne_us'))==-1;
    Une_us(:,xi) = -Une_us(:,xi);
    Vne_us(xi,:) = -Vne_us(xi,:);
    Une_usNorm = nan(size(Une_us));
    for j = 1:nPCs
        Une_usNorm(:,j) = (Une_us(:,j) - min(Une_us(:,j)))/(max(Une_us(:,j))-min(Une_us(:,j)));
    end
    naturalsound_acoustic_category_figures(Une_usNorm, individual_subject_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-Norm-'],'corr_bounds',[-1 1]);
    
end

%% Acoustic scatter plots
[R2, C, conds] = naturalsound_acoustic_category_regressors;
acoustic_scatter(UicaNewNorm(:,1), R2, C, '400', [0 1]);
export_fig([ica_directory 'LowFreqCorr.pdf'],'-pdf','-transparent');
acoustic_scatter(UicaNewNorm(:,2)+5, R2, C, 'Res', [0 1]+5);
export_fig([ica_directory 'ResCorr.pdf'],'-pdf','-transparent');

%% 3D ICA and PCA

figure;
[U,S,V] = svd(v,'econ');
nPCs = 3;
Vr = V(:,1:nPCs)'*sqrt(size(V,1)-1);
Ur = U(:,1:nPCs)/sqrt(size(V,1)-1);
Sr = S(1:nPCs,1:nPCs);
naturalsound_acoustic_category_figures(Ur, figure_directory, ['pca-' num2str(nPCs) '-RandSeed' num2str(i) '-']);

% Kurtosis
figure;
subplot(1,2,1);
bar(diag(Sr).^2,'FaceColor',[1 1 1]);
xlabel('PCA Components'); ylabel('Explained Variance');
kurtICA = abs(kurtosis(Vr')-3);
subplot(1,2,2);
bar(kurtICA,'FaceColor',[1 1 1]);
xlabel('PCA Components'); ylabel('Kurtosis');
export_fig([figure_directory 'PCA-' num2str(nPCs) '_exvar_kurtosis.pdf'],'-pdf','-nocrop');

% Histogram
plothist_bysubject(Vr, si, [figure_directory 'pca-' num2str(nPCs) '-hist-bysubject']);

% ICA
for i = 1:1
    ResetRandStream(i);
    [icasig, A, W] = fastica(Vr);
    rotMatICA = W;
    Uica = Ur*Sr*A;
    Vica = icasig;
    subplot(1,2,1);
    bar(var(Uica),'FaceColor',[1 1 1]);
    xlabel('ICA Components'); ylabel('Explained Variance');
    kurtICA = abs(kurtosis(Vica')-3);
    subplot(1,2,2);
    bar(kurtICA,'FaceColor',[1 1 1]);
    xlabel('ICA Components'); ylabel('Kurtosis');
    export_fig([figure_directory 'ICA-' num2str(nPCs) '-RandSeed' num2str(i) '_exvar_kurtosis.pdf'],'-pdf','-nocrop');
    % x = Ur*Sr*A;
    % response = (x-ones(size(x,1),1)*min(x)) ./ (ones(size(x,1),1)*(max(x)-min(x)));
    naturalsound_acoustic_category_figures(Uica, figure_directory, ['ica-' num2str(nPCs) '-RandSeed' num2str(i) '-']);
end
close all;

% histogram by subject for pca
plothist_bysubject(Vica, si, [figure_directory 'ica-' num2str(nPCs) '-hist-bysubject']);

% 2D Kurtosis
n = 40;
th = linspace(-pi/2,pi/2,n);
phi = linspace(-pi/2,pi/2,n);
kurt = nan(n,n);
for i = 1:n
    for j = 1:n
        [x,y,z] = sph2cart(th(i),phi(j),1);
        kurt(i,j) = abs(kurtosis( ([x,y,z]*Vr)' )-3);
    end
end
kurtnorm = norm01(kurt);
figure;
hold on;
c = colormap('jet(100)');
for i = 1:size(kurt,1)
    for j = 1:size(kurt,2)
        [x,y,z] = sph2cart(th(i),phi(j),1);
        plot(y,z,'o','Color', c(round(kurtnorm(i,j)*99)+1,:), 'LineWidth',5);
    end
end
h = colorbar; set(h,'YTick',10:10:100,'YTickLabel',interp1(1:100, linspace(min(kurt(:)),max(kurt(:)),100), 10:10:100));
for i = 1:3
    plot(rotMatICA(i,2)*sign(rotMatICA(i,1)),rotMatICA(i,3)*sign(rotMatICA(i,1)),'ko','LineWidth',5);
end
export_fig([figure_directory 'ICA-' num2str(nPCs) '_kurtosis2D.pdf'],'-pdf','-nocrop');

% 3D Kurtosis plot
n = 40;
th = linspace(-pi/2,pi/2,n)*2;
phi = linspace(-pi/2,pi/2,n)*2;
kurt = nan(n,n);
for i = 1:n
    for j = 1:n
        [x,y,z] = sph2cart(th(i),phi(j),1);
        kurt(i,j) = abs(kurtosis( ([x,y,z]*Vr)' )-3);
    end
end
kurtnorm = norm01(kurt);
figure;
hold on;
c = colormap('jet(100)');
for i = 1:size(kurt,1)
    for j = 1:size(kurt,2)
        [x,y,z] = sph2cart(th(i),phi(j),1);
        plot3(x,y,z,'o','Color', c(round(kurtnorm(i,j)*99)+1,:), 'LineWidth',10);
    end
end
h = colorbar; set(h,'YTick',10:10:100,'YTickLabel',interp1(1:100, linspace(min(kurt(:)),max(kurt(:)),100), 10:10:100));
%
% localmin = [0.8597, 0.4964, -0.1205];
% naturalsound_acoustic_category_figures(Ur*Sr*localmin', figure_directory, ['localmin-']);

%% 2D ICA and PCA

figure;
[U,S,V] = svd(v,'econ');
nPCs = 2;
Vr = V(:,1:nPCs)'*sqrt(size(V,1)-1);
Ur = U(:,1:nPCs)/sqrt(size(V,1)-1);
Sr = S(1:nPCs,1:nPCs);

% Kurtosis
figure;
subplot(1,2,1);
bar(diag(Sr).^2,'FaceColor',[1 1 1]);
xlabel('PCA Components'); ylabel('Explained Variance');
kurtICA = abs(kurtosis(Vr')-3);
subplot(1,2,2);
bar(kurtICA,'FaceColor',[1 1 1]);
xlabel('PCA Components'); ylabel('Kurtosis');
export_fig([figure_directory 'PCA-' num2str(nPCs) '_exvar_kurtosis.pdf'],'-pdf','-nocrop');

% Histogram
plothist_bysubject(Vr, si, [figure_directory 'pca-' num2str(nPCs) '-hist-bysubject']);

% ICA
for i = 1:1
    ResetRandStream(i);
    [icasig, A, W] = fastica(Vr);
    icasig(1,:) = -icasig(1,:);
    A(:,1) = -A(:,1);
    rotMatICA = pinv(A);
    Uica = Ur*Sr*A;
    Vica = icasig;
    subplot(1,2,1);
    bar(var(Uica),'FaceColor',[1 1 1]);
    xlabel('ICA Components'); ylabel('Explained Variance');
    kurtICA = abs(kurtosis(Vica')-3);
    subplot(1,2,2);
    bar(kurtICA,'FaceColor',[1 1 1]);
    xlabel('ICA Components'); ylabel('Kurtosis');
    export_fig([figure_directory 'ICA-' num2str(nPCs) '-RandSeed' num2str(i) '_exvar_kurtosis.pdf'],'-pdf','-nocrop');
    % x = Ur*Sr*A;
    % response = (x-ones(size(x,1),1)*min(x)) ./ (ones(size(x,1),1)*(max(x)-min(x)));
    naturalsound_acoustic_category_figures(Uica, figure_directory, ['ica-' num2str(nPCs) '-RandSeed' num2str(i) '-']);
end
close all;

% histogram by subject for pca
plothist_bysubject(Vica, si, [figure_directory 'ica-' num2str(nPCs) '-hist-bysubject']);

% 2D Kurtosis
n = 100;
th = linspace(0,2*pi,n);
kurt = nan(n,1);
for i = 1:n
    [x,y] = pol2cart(th(i),1);
    kurt(i) = abs(kurtosis( ([x,y]*Vr)' )-3);
end
kurtnorm = norm01(kurt);

figure
[x,y] = pol2cart(th',kurt);
plot(x,y,'k-o','LineWidth',2);
hold on;
quiver(zeros(2,1),zeros(2,1),rotMatICA(:,1),rotMatICA(:,2),'r','LineWidth',2);
xlim(max(kurt(:))*[-1 1]);ylim(max(kurt(:))*[-1 1]);
export_fig([figure_directory 'ICA-' num2str(nPCs) '_kurtosis2D.pdf'],'-pdf','-nocrop');

% localmin = [0.3327, 0.2138];
% localmin = localmin/norm(localmin);
% naturalsound_acoustic_category_figures(Ur*Sr*localmin', figure_directory, ['localmin-2D-']);

%% Mask

surface_projections_directory = [figure_directory 'surface_projections/'];
if ~exist(surface_projections_directory,'file')
    mkdir(surface_projections_directory);
end

% initialize variables
sound_responseP = zeros(nsurfpts, length(usubs), 2);
sound_response = zeros(nsurfpts, length(usubs), 2);
soundthresh = 3;
min_subjects = 3;

usubs = [1 11 106 121 122 123 143 160 161 164];
% grandmean = mean(v,2);
% hemispheres
hemis = {'rh','lh'};
runs = 1:2;
for i = 1:2
    
    mask = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.cortex.label']);
    for j = 1:length(usubs)
        
        % subject id
        fprintf('hemi %s, usub %d\n', hemis{i}, usubs(j)); drawnow;
        subjid = ['naturalsound_us' num2str(usubs(j))];
        
        % surface directory
        surf_directory = [params('rootdir') 'freesurfer/fsaverage/linear_analysis_surf/' subjid '/'];
        if ~exist(surf_directory,'dir')
            mkdir(surf_directory);
        end
        
        % voxel responses
        for q = 1:length(runs)
            x = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r' num2str(runs(q)) '/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.mgz']);
            dat = squeeze(x.vol)';
            if q == 1
                sv = zeros(size(dat));
            end
            sv(:,mask.vnums+1) = sv(:,mask.vnums+1) + dat(:,mask.vnums+1)/length(runs);
        end
        
        % sound response
        [~,p] = ttest(sv,0,'tail','right');
        sound_responseP(:,j,i) = -log10(p);
        sound_response(:,j,i) = mean(sv);
        
    end
    sound_responseP_ovmap = squeeze(sum(sound_responseP(:,:,i) > soundthresh,2));
    fname = [surface_projections_directory hemis{i} '.overlapmap.mgz'];
    MRIwrite_SNH(sound_responseP_ovmap, fname, hemis{i});
    %     freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[1, length(usubs)/4, length(usubs)/2],'label',[params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.hand-audctx.label']);
    freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',strrep(fname,'.mgz','.png'),'overlay_threshold',[1, length(usubs)/4, length(usubs)/2],'label',[params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.hand-audctx.label']);
    
    sound_response_psc = 100*squeeze(mean(sound_response(:,:,i),2));
    fname = [surface_projections_directory hemis{i} '.sound_response_psc.mgz'];
    MRIwrite_SNH(sound_response_psc, fname, hemis{i});
    %     freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[0.01, 0.15, 0.30],'piecewise','label',[params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.hand-audctx.label']);
    freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',strrep(fname,'.mgz','.png'),'overlay_threshold',[0.01, 0.15, 0.30],'piecewise','label',[params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.hand-audctx.label']);
end

for i = 1:2
    freeview3('fsaverage',hemis{i},'label',[params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.hand-audctx.label']);
end



%% Select stimuli with maximum variance v2

conds = read_conditions('naturalsound',1,'main_v3_combined');
n_stimuli = 40;
UneZscore = zscore(Une);
remaining_stimuli = 1:165;

for i = 1:165-n_stimuli
    
    n = length(remaining_stimuli);
    mean_variance = nan(1,n);
    for j = 1:n
        x = UneZscore(remaining_stimuli(setdiff(1:n,j)), :);
        mean_variance(j) = mean(var(x));
    end
    
    [~,xi] = max(mean_variance);
    remaining_stimuli = setdiff(remaining_stimuli, remaining_stimuli(xi));
end

selected_stimuli = remaining_stimuli;
component_names = {'inv-pitch','speech','lowfreq','music','envsounds'}; %{'pitch','speech','envsounds','music','lowfreq'};
figure;
set(gcf,'Position', [0 0 1400 800]);
subplot(1,3,1);
[y,bins] = hist(UneZscore, 10);
plot(bins,y/165,'LineWidth',2);
xlabel('Z values');
title('Histogram of Stimulus Values');
legend(component_names,'Location','Best');
set(gca, 'FontSize', 14);
box off;

subplot(1,3,2);
[y,bins] = hist(UneZscore(selected_stimuli,:), 10);
plot(bins,y/sum(y(:,1)),'LineWidth',2);
xlabel('Z values');
title('Histogram of Stimulus Values for Selected Subset')
set(gca, 'FontSize', 14);
box off;

subplot(1,3,3);
x = UneZscore(selected_stimuli,:);
[y,bins] = hist(x(:), 10);
plot(bins,y/sum(y(:,1)),'k','LineWidth',2);
xlabel('Z values');
title('All Stimulus Values for Selected Set')
set(gca, 'FontSize', 14);
box off;
export_fig([figure_directory 'histogram_stimulus_weights_' num2str(n_stimuli) 'stims.pdf'],'-pdf','-nocrop','-transparent');

% category figures
figure;
set(gcf,'Position', [0 0 1400 800]);
subplot(1,2,1);
[R2, C, conds] = naturalsound_acoustic_category_regressors;
[y,bins] = hist(C.category_assignments,11);
hold on;
for i = 1:11
    cat = C.plotting_order(i);
    bar(i,y(cat)/165,'FaceColor',C.colors(cat,:));
end
xlim([0 12]);
title('All Stimuli');
legend(C.category_labels(C.plotting_order),'Location','Best');
set(gca, 'FontSize', 14);
box off;

subplot(1,2,2);
[y,bins] = hist(C.category_assignments(selected_stimuli),11);
hold on;
for i = 1:11
    cat = C.plotting_order(i);
    bar(i,y(cat)/n_stimuli,'FaceColor',C.colors(cat,:));
end
xlim([0 12]);
title('Selected Stimuli');
set(gca, 'FontSize', 14);
box off;
export_fig([figure_directory 'category_histogram_' num2str(n_stimuli) 'stims.pdf'],'-pdf','-nocrop','-transparent');

stim_names = conds(selected_stimuli);
ids = C.ids(selected_stimuli);
save([figure_directory 'selected_stimuli_' num2str(n_stimuli) 'stims.mat'],'stim_names','ids','Csubset','R2subset');


stim_suffix = cell(1,60);
for i = 1:60
    x = regexp(stim_names{i},'_');
    stim_suffix{i} = stim_names{i}(x+1:end);
end

%%
[R2, C, conds] = naturalsound_acoustic_category_regressors('gammatone_Spec_SpecTempMod');


ids_texture = [28 29 33 41 72 78 80 83 86 92 102 150 211 224 232 254 268 283 387 394 401 414 418 461 501 503 504 506 508 512];
% ids_texture = [11 28 29 33 41 72 78 80 83 86 92 97 102 150 211 224 232 254 268 283 315 320 337 387 394 401 414 418 461 501 502 503 504 506 508 512];
selected_stimuli = find(ismember(C.ids, ids_texture));
stims = conds(selected_stimuli);

% acoustic and category figures for subset of stimuli
Csubset = C;
Csubset.category_assignments = C.category_assignments(selected_stimuli);
Csubset.ids = C.ids(selected_stimuli);
Csubset.category_regressors = C.category_regressors(selected_stimuli,:);

R2subset = R2;
R2subset.ids = R2.ids(selected_stimuli);
R2subset.stats = R2.stats(selected_stimuli,:);

conds_subset = conds(selected_stimuli);

for i = 1:size(UneNorm,2)
    y = UneNorm(selected_stimuli,i);
    response_profile_figure(y, Csubset, conds_subset);
    export_fig([figure_directory 'ICA-5-handpicked-subset-' num2str(length(selected_stimuli)) 'stims-' component_names{i} '_response_profile.eps'],'-eps','-nocrop','-transparent');
    %     acoustic_category_figures(y, Csubset, R2subset, [figure_directory 'ICA-5-stimulus-subset-' num2str(length(selected_stimuli)) 'stims-' component_names{i}''], varargin{:});
end

%% Grid smooth
% 
% 
% 
% % gaussian kernel for smoothing
% if fwhm_surf > 0
%     x = round(3*fwhm_surf/2);
%     x = x-mod(x,2)+1; % make odd
%     gaussian_kernel = fspecial('gaussian', [x,x], fwhm2std(fwhm_surf/2));
%     gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
% end
% 
% % smoothed data
% vgrid_smooth = cell(size(vgrid));
% for i = 1:length(hemis)
%     for j = 1:length(usubs)
%         fprintf('%d %s\n', usubs(j), hemis{i}); drawnow;
%         for q = 1:165
%             runs = read_runs('naturalsound', usubs(j), 'main_v3_combined');
%             for z = 1:length(runs)
%                 vgrid_smooth{i}(:,:,j,q,z) = conv2_setNaNs_tozero(vgrid{i}(:,:,j,q,z), gaussian_kernel); 
%             end
%         end
%     end
% end
% 
% % roll-out the grid of voxels for each sound/subject/scan to a vector
% vgrid_smooth_reshape_rh = reshape(vgrid_smooth{1}, [n_voxels_rh, length(usubs), 165, maxruns]);
% vgrid_smooth_reshape_lh = reshape(vgrid_smooth{2}, [n_voxels_lh, length(usubs), 165, maxruns]);
% 
% % combine right and left hemisphere voxels
% n_voxels = n_voxels_lh + n_voxels_rh;
% v_smooth_hemi_combined = cat(1, vgrid_smooth_reshape_rh, vgrid_smooth_reshape_lh);
% 
% % combine across subjects
% % each row: [rh_vox_s1, lh_vox_s1, rh_vox_s2, lh_vox_s2, ...]
% v_smooth_unfilt_allscans = nan([165, n_voxels*length(usubs), maxruns]);
% for i = 1:maxruns
%     v_smooth_unfilt_allscans(:,:,i) = reshape(v_smooth_hemi_combined(:,:,:,i), [n_voxels*length(usubs), 165])';
% end
% v_smooth_unfilt_mean = nanmean(v_smooth_unfilt_allscans,3);
% 
% % p-value sound threshold
% [~,x] = ttest(v_smooth_unfilt_mean, 0, 0.05, 'right');
% soundP_smooth = -log10(x);
% soundP_smooth_individual_scans = nan(length(soundP_smooth), maxruns);
% for i = 1:maxruns
%     [~,x] = ttest(v_smooth_unfilt_allscans(:,:,i), 0, 0.05, 'right');
%     soundP_smooth_individual_scans(:,i) = -log10(x);
% end
% 
% % cross-scan regression
% regressExVar_smooth = nan(size(soundP_smooth));
% for i = 1:size(v_smooth_unfilt_allscans,2)
%     %     if mod(i,1000)==0
%     %         fprintf('%d\n',i);
%     %     end
%     x = v_smooth_unfilt_allscans(:,i,1);
%     y = v_smooth_unfilt_allscans(:,i,2);
%     [~,~,~,regressE] = regress_SNH(x,y);
%     regressExVar_smooth(i) = 100*(sum(y.^2) - sum(regressE.^2)) / sum(y.^2);
% end
%% PCA eigen value spectrum
% [U,S,V] = svd(v,'econ');
% eigs = diag(S).^2;
% exvar = eigs/sum(eigs);
% figure;
% bar(1:20,100*exvar(1:20),'FaceColor',[1 1 1]);
% xlabel('PCA Components'); ylabel('% Explained Variance');
% export_fig([figure_directory 'PCA_eigspec.pdf'],'-pdf','-nocrop');
% naturalsound_acoustic_category_figures(U(:,1:5), figure_directory, 'pca-');
% close all;

%% Explained variance of with-subject PCAs in left-out Data
% r = nan(length(usubs),165);
% for i = 1:length(usubs)
%     fprintf('subject %d\n',usubs(i)); drawnow;
%     %     [U_loo,~,~] = svd(v(:,si ~= usubs(i)),'econ');
%     scan1 = v1(:,si == usubs(i));
%     scan2 = v2(:,si == usubs(i));
%     [U_loo,~,~] = svd(scan1,'econ');
%     xhat = U_loo'*scan1;
%     for j = 1:165
%         prediction = U_loo(:,1:j)*xhat(1:j,:);
%         r(i,j) = fastcorr3(prediction(:), scan2(:));
%     end
% end
%
% r2 = [r.^2];
% figure;
% plot(r2','-o');
% ylim([0 max(r2(:))*1.2]); xlim([0 165]);
% xlabel('Number of PCs'); ylabel('Explained Variance in Left-Out Data');
% export_fig([figure_directory 'PCA_withinsubject_exvar_individual_subjects.pdf'],'-pdf','-nocrop');
%
% figure;
% plot(mean(r2),'-');
% errorbar(1:165,mean(r2),stderr_withsub(r2),'k.-');
% ylim([0 max(mean(r2))*1.2]); xlim([0 165]);
% xlabel('Number of PCs'); ylabel('Explained Variance in Left-Out Data');
% export_fig([figure_directory 'PCA_withinsubject_loso_exvar_mean.pdf'],'-pdf','-nocrop');
%
% figure;
% bar(mean(r2(:,1:20)),'FaceColor',[1 1 1]);
% hold on;
% errorbar(1:20,mean(r2(:,1:20)),stderr_withsub(r2(:,1:20)),'k.');
% ylim([0 max(mean(r2))*1.2]); xlim([0 21]);
% xlabel('Number of PCs'); ylabel('Explained Variance in Left-Out Data');
% export_fig([figure_directory 'PCA_withinsubject_loso_exvar_mean_20PCs.pdf'],'-pdf','-nocrop');
%
% figure;
% r2diff = diff([zeros(length(usubs),1),r2],[],2);
% bar(mean(r2diff(:,1:20)),'FaceColor',[1 1 1]);
% hold on;
% errorbar(1:20,mean(r2diff(:,1:20)),stderr_withsub(r2diff(:,1:20)),'k.');
% xlabel('PCs');ylabel('Change in Explained Variance in Data from Leftout Subject');
% export_fig([figure_directory 'PCA_withinsubject_loso_change_exvar_20PCs.pdf'],'-pdf','-nocrop');
% close all;

%% Explained variance of pca components left-out data

% [U_loo,~,~] = svd(v1,'econ');
% xhat = U_loo'*v1;
% r = nan(1,165);
% for j = 1:165
%     prediction = U_loo(:,1:j)*xhat(1:j,:);
%     r(j) = fastcorr3(prediction(:), v2(:));
% end
%
% r2 = [r.^2];
% figure;
% plot(r2','-o');
% ylim([0 max(r2(:))*1.2]); xlim([0 165]);
% xlabel('Number of PCs'); ylabel('Explained Variance in Left-Out Data');
% export_fig([figure_directory 'PCA_exvar_leftout_scan.pdf'],'-pdf','-nocrop');
%
% figure;
% bar(r2(1:20)','FaceColor',[1 1 1]);
% ylim([0 max(r2(:))*1.2]); xlim([0 21]);
% xlabel('Number of PCs'); ylabel('Explained Variance in Left-Out Data');
% export_fig([figure_directory 'PCA_exvar_leftout_scan_20PCs.pdf'],'-pdf','-nocrop');
%
% figure;
% r2diff = diff([0,r2],[],2);
% bar(r2diff(:,1:20),'FaceColor',[1 1 1]);
% hold on;
% xlabel('PCs');ylabel('Change in Explained Variance in Data from Leftout Subject');
% export_fig([figure_directory 'PCA_change_exvar_leftout_scan_20PCs.pdf'],'-pdf','-nocrop');
% close all;
%% Regress IC components in whole brain OLD
%
% % block timing, runs
% blocktp = 3:6;
% runs = 1:2;
% idstring = 'ica-5-norm-demean-';
%
% surface_projections_directory = [figure_directory 'surface_projections/'];
% if ~exist(surface_projections_directory,'file')
%     mkdir(surface_projections_directory);
% end
%
% % initialize variables
% sound_responseP = nan(nsurfpts, length(usubs), 2);
% % ica_weights = nan(5, nsurfpts, length(usubs), 2);
% ica_weights_soundthresh = nan(5, nsurfpts, length(usubs), 2);
% ica_logP = nan(5, nsurfpts, length(usubs), 2);
% ica_weights_group = nan(5, nsurfpts, 2);
% % ica_logP_group = nan(5, nsurfpts, 2);
% soundthresh = 3;
% min_subjects = 3;
%
% % hemispheres
% hemis = {'rh','lh'};
% for i = 1:2
%
%     mask = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.audctx.label']);
%     for j = 1:length(usubs)
%
%         % surface directory
%         surf_directory = [params('rootdir') 'freesurfer/fsaverage/fla/linear_analysis_surf/usub' num2str(usubs(j)) '/'];
%         if ~exist(surf_directory,'dir')
%             mkdir(surf_directory);
%         end
%
%         % subject id
%         fprintf('hemi %s, usub %d\n', hemis{i}, usubs(j)); drawnow;
%         subjid = ['naturalsound_us' num2str(usubs(j))];
%
%         % voxel responses
%         x1 = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r1/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.mgz']);
%         x2 = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r2/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.mgz']);
%         x3 = squeeze((x1.vol + x2.vol)/2)';
%         sv = zeros(size(x3));
%         sv(:,mask.vnums+1) = x3(:,mask.vnums+1);
%
%         % sound response
%         [~,p] = ttest(sv,0,'tail','right');
%         sound_responseP(:,j,i) = -log10(p);
%
%         % select voxels with a significant sound response
%         sv_selected = sv(:,sound_responseP(:,j,i) > soundthresh);
%
%         % normed and demeaned
%         sv_norm_demeaned = sumnorm_demean_subjects(sv_selected, ones(1,size(sv_selected,2)));
%
%         % Regression analysis
%         Uica_norm = Uica ./ (ones(165,1)*std(Uica));
%         [b, ~, logP, ~] = regress_SNH(Uica_norm, zscore(sv));
%
%         % assign betas
%         b(isnan(b)) = 0;
%         b(:,sound_responseP(:,j,i) < soundthresh) = 0;
%         ica_weights_soundthresh(:,:,j,i) = b;
%         weight_std = std(b(:,sound_responseP(:,j,i) > soundthresh)');
%
%         % assign p-values
%         logP(isnan(logP)) = 0;
%         logP(:,sound_responseP(:,j,i) < soundthresh) = 0;
%         ica_logP(:,:,j,i) = logP;
%         %         b_sound_thresh = b;
%         %         b_sound_thresh(:,sound_responseP(:,j,i) < soundthresh) = 0;
%         %         ica_weights(:,:,j,i) = b;
%
%         % write files
%         for k = 1:size(Uica,2)
%
%             % write beta weights
%             %             fname = [surf_directory idstring 'beta-component-' num2str(k) '-' hemis{i} '.mgz'];
%             %             MRIwrite_SNH(ica_weights(:,:,j,i), fname, hemis{i});
%
%             % write beta weights thresholded by sound responsivity
%             fname = [surf_directory idstring 'beta-soundthresh' num2str(soundthresh) '-component' num2str(k) '-' hemis{i} '.mgz'];
%             MRIwrite_SNH(ica_weights_soundthresh(k,:,j,i), fname, hemis{i});
%
%             % save screenshot
%             surface_image = [surface_projections_directory idstring 'beta-soundthresh'  num2str(soundthresh) '-component' num2str(k) '-us' num2str(usubs(j)) '-' hemis{i} '.png'];
%             if ~exist(surface_image, 'file')
%                 freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.00001 0.00002 weight_std(k)*3],'piecewise');
%                 %                 freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[0.00001 0.00002 0.5],'piecewise');
%             end
%
%             % write beta weights thresholded by sound responsivity
%             fname = [surf_directory idstring 'logP-component' num2str(k) '-' hemis{i} '.mgz'];
%             MRIwrite_SNH(ica_logP(k,:,j,i), fname, hemis{i});
%
%             % save screenshot
%             surface_image = [surface_projections_directory idstring 'logP-component' num2str(k) '-us' num2str(usubs(j)) '-' hemis{i}  '.png'];
%             if ~exist(surface_image, 'file')
%                 freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[3 4.5 6]);
%                 %                 freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[3 4.5 6]);
%             end
%
%             %             ica_weight_surf = x1;
%             %             dims = size(x1.vol);
%             %             ica_weight_surf.vol = nan([dims(1:3),1]);
%             %             ica_weight_surf.vol(1,:,1,1) = ica_weights(k,:,j,i)';
%             %             ica_weight_surf.nframes = size(ica_weight_surf.vol,4);
%             %             ica_weight_surf.fspec = surface_file;
%             %             MRIwrite(ica_weight_surf, surface_file);
%             %             surface_file = [surf_directory idstring 'logP-component-' num2str(k) '.mgz'];
%             %             ica_logP_surf = x1;
%             %             dims = size(x1.vol);
%             %             ica_logP_surf.vol = nan(dims(1:3));
%             %             ica_logP_surf.vol(1,:,1) = ica_logP(k,:,j,i)';
%             %             ica_logP_surf.nframes = size(ica_logP_surf.vol,4);
%             %             ica_logP_surf.fspec = surface_file;
%             %             MRIwrite(ica_logP_surf, surface_file);
%
%         end
%     end
%
%     % surface directory
%     surf_directory = [params('rootdir') 'freesurfer/fsaverage/fla/linear_analysis_surf/us' sprintf('-%d',usubs) '/'];
%     if ~exist(surf_directory,'dir')
%         mkdir(surf_directory);
%     end
%
%     x1 = ica_weights_soundthresh(:,:,:,i);
%     x1(x1==0) = NaN;
%     x2 = nanmean(x1,3);
%     x2(isnan(x2)) = 0;
%     ica_weights_group(:,:,i) = x2;
%     xi = sum(sound_responseP(:,:,i) > soundthresh,2) < min_subjects - 1e-3;
%     ica_weights_group(:,xi,i) = 0;
%     weight_std = std(ica_weights_group(:,~xi,i)');
%
%     % group results
%     for k = 1:size(Uica,2)
%
%         % group beta figure
%         fname = [surf_directory idstring 'beta-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-component' num2str(k) '-' hemis{i} '.mgz'];
%         MRIwrite_SNH(ica_weights_group(k,:,i), fname, hemis{i});
%
%         % save screenshot
%         surface_image = [surface_projections_directory idstring 'beta-soundthresh' num2str(soundthresh) '-minsubjects' num2str(min_subjects) '-component' num2str(k) '-group-' hemis{i} '.png'];
%         if true || ~exist(surface_image, 'file')
%             freeview3('fsaverage',hemis{i},'overlay',fname,'screenshot',surface_image,'overlay_threshold',[0.00001 0.00002 weight_std(k)*3],'piecewise');
%             %             freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[0.00001 0.00002 0.2],'piecewise');
%         end
%
%         %         [~,p] = ttest(ica_weights(:,:,:,i),0,'dim',3);
%         %         p(isnan(p)) = 0;
%         %         ica_logP_group(:,:,i) = -log10(p);
%         %
%         %         fname = [surf_directory idstring 'logP-component-' num2str(k) '-' hemis{i} '.mgz'];
%         %         MRIwrite_SNH(ica_logP_group(:,:,i), fname, hemis{i});
%         %
%         %         % save screenshot
%         %         surface_image = [surface_projections_directory idstring 'logP-component-' num2str(k) '-group-' hemis{i} '.png'];
%         %         if false || ~exist(surface_image, 'file')
%         %             %                 freeview3('fsaverage',hemis{i},'overlay',surface_file,'screenshot',surface_image,'overlay_threshold',[0.001 1.5 3],'piecewise');
%         %             freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[3 4.5 6]);
%         %         end
%
%         %         surface_file = ;
%         %         ica_weight_surf = x1;
%         %         dims = size(x1.vol);
%         %         ica_weight_surf.vol = nan(dims(1:3));
%         %         ica_weight_surf.vol(1,:,1) = ica_weights_group(:,:,i)';
%         %         ica_weight_surf.nframes = size(ica_weight_surf.vol,4);
%         %         ica_weight_surf.fspec = surface_file;
%         %         MRIwrite(ica_weight_surf, surface_file);
%
%         % save screenshot
%         %         surface_image = [surface_projections_directory idstring 'logP-component-' num2str(k) '-us' num2str(usubs) '.png'];
%         %         if false || ~exist(surface_image, 'file')
%         %             %                 freeview3('fsaverage',hemis{i},'overlay',surface_file,'screenshot',surface_image,'overlay_threshold',[0.001 1.5 3],'piecewise');
%         %             freeview3('fsaverage',hemis{i},'overlay',surface_file,'overlay_threshold',[3 4.5 6]);
%         %         end
%         %
%         %         surface_file = [surf_directory idstring 'logP-component-' num2str(k) '.mgz'];
%         %         ica_logP_surf = x1;
%         %         dims = size(x1.vol);
%         %         ica_logP_surf.vol = nan(dims(1:3));
%         %         ica_logP_surf.vol(1,:,1) = ica_logP(k,:,j,i)';
%         %         ica_logP_surf.nframes = size(ica_logP_surf.vol,4);
%         %         ica_logP_surf.fspec = surface_file;
%         %         MRIwrite(ica_logP_surf, surface_file);
%         %
%         %         % save screenshot
%         %         surface_image = [surface_projections_directory idstring 'logP-component-' num2str(k) '-us' num2str(usubs) '.png'];
%         %         if false || ~exist(surface_image, 'file')
%         %             %                 freeview3('fsaverage',hemis{i},'overlay',surface_file,'screenshot',surface_image,'overlay_threshold',[0.001 1.5 3],'piecewise');
%         %             freeview3('fsaverage',hemis{i},'overlay',surface_file,'overlay_threshold',[3 4.5 6]);
%         %         end
%         %
%     end
% end



%% Regress pitch and music and low-frequencies

[R2, C, conds, R3] = naturalsound_acoustic_category_regressors;
regressors = zscore([R2.stats(:,ismember(R2.stat_names,{'Res'})), any(C.category_regressors(:,ismember(C.category_labels,{'EngSpeech','ForSpeech'})),2), any(C.category_regressors(:,ismember(C.category_labels,{'Song','Music'})),2)]);
regressor_names = {'Res','Speech','Music'};

[pitchmusB,se,pitchmusP,e] = regress_SNH(regressors,v);

% convert to grid representation
pitchmusB_grid = cell(size(pitchmusB,1), 2);
pitchmusB_grid_group = cell(size(pitchmusB,1), 2);
for i = 1:size(pitchmusB,1);
    % n_voxels x subject matrix
    x1 = nan(n_voxels*length(usubs),1);
    x1(voxel_selection) = pitchmusB(i,:)';
    x2 = reshape(x1, [n_voxels, length(usubs)]);
    
    % split left and right hemispheres, and convert back to a grid
    pitchmusB_grid{i,1} = reshape(x2(1:n_voxels_rh, :), [size(grid_x{1}), length(usubs)]);
    pitchmusB_grid{i,2} = reshape(x2((1:n_voxels_lh) + n_voxels_rh, :), [size(grid_x{2}), length(usubs)]);
    pitchmusB_grid_group{i,1} = nanmean(pitchmusB_grid{i,1},3);
    pitchmusB_grid_group{i,2} = nanmean(pitchmusB_grid{i,2},3);
    pitchmusB_grid_group{i,1}(sum(~isnan(pitchmusB_grid{i,1}),3) < 2.001) = NaN;
    pitchmusB_grid_group{i,2}(sum(~isnan(pitchmusB_grid{i,2}),3) < 2.001) = NaN;
end

% grid directory
grid_directory = [figure_directory 'Res-Speech-Music/'];
if ~exist(grid_directory, 'dir')
    mkdir(grid_directory);
end

% plot ic grids
for i = 1:size(pitchmusB,1);
    close all;
    figure(1);clf(1);
    set(gcf,'Position',[0 0 1000, 500])
    subplot(1,2,1);
    %     imagesc(fliplr(ic_group{i,2}),[-2 2]);
    surf(grid_x{2}, grid_y{2}, pitchmusB_grid_group{i,2});
    view(0,90);
    title('Left Hemisphere');
    subplot(1,2,2);
    surf(grid_x{1}, grid_y{1}, pitchmusB_grid_group{i,1});
    view(0,90);
    title('Right Hemisphere');
    drawnow;
    export_fig([grid_directory regressor_names{i} '-group-vras-coordinates.png'],'-png');
    %     print(gcf, '-dpng', [grid_directory 'ic-' num2str(i) '-group-vras-coordinates.png'], '-r100');
    
    close all;
    figure(1);clf(1);
    set(gcf,'Position',[0 0 1000, 500]);
    subplot(1,2,1);
    imagesc(rot90(pitchmusB_grid_group{i,2},3),[-1.5 1.5]*1e-3);
    %     surf(grid_x{2}, grid_y{2}, ic_group{i,2});
    title('Left Hemisphere');
    subplot(1,2,2);
    imagesc(rot90(pitchmusB_grid_group{i,1},-1),[-1.5 1.5]*1e-3);
    title('Right Hemisphere');
    drawnow;
    export_fig([grid_directory regressor_names{i} '-group.png'],'-png');
    %     print(gcf, '-dpng', [grid_directory 'ic-' num2str(i) '-group.png'], '-r100');
    
    for j = 1:length(usubs)
        close all;
        figure(1);clf(1);
        set(gcf,'Position',[0 0 1000, 500]);
        subplot(1,2,1);
        imagesc(rot90(pitchmusB_grid{i,2}(:,:,j),3),[-1.5 1.5]*1e-3);
        %     surf(grid_x{2}, grid_y{2}, ic_group{i,2});
        title('Left Hemisphere');
        subplot(1,2,2);
        imagesc(rot90(pitchmusB_grid{i,1}(:,:,j),-1),[-1.5 1.5]*1e-3);
        title('Right Hemisphere');
        drawnow;
        %         print(gcf, '-dpng', [grid_directory 'ic-' num2str(i) '-us' num2str(usubs(j)) '.png'], '-r100');
        export_fig([grid_directory regressor_names{i} '-us' num2str(usubs(j)) '.png'],'-png');
    end
end

%
% % plot ic grids
% for i = 1:size(Vica,1);
%     figure(1);clf(1);
%     subplot(1,2,1);
%     imagesc(fliplr(pitchmusB_grid_group{i,2}),[-2 2]);
%     title('Right Hemisphere');
%     subplot(1,2,2);
%     imagesc(flipud(pitchmusB_grid_group{i,1}),[-2 2]);
%     title('Left Hemisphere');
%     export_fig([grid_directory 'ic-' num2str(i) '-group.png'],'-png','-r100');
%
%     for j = 1:length(usubs)
%         figure(1);clf(1);
%         subplot(1,2,1);
%         imagesc(fliplr(pitchmusB_grid{i,2}(:,:,j)),[-2 2]);
%         title('Right Hemisphere');
%         subplot(1,2,2);
%         imagesc(flipud(pitchmusB_grid{i,1}(:,:,j)),[-2 2]);
%         title('Left Hemisphere');
%         export_fig([grid_directory 'ic-' num2str(i) '-us' num2str(usubs(j)) '.png'],'-png','-r100');
%     end
% end

%%

% plot ic grids
x = reshape(voxel_selection, [n_voxels, length(usubs)]);
sound_response_mask{1} = sum(reshape(x(1:n_voxels_rh, :), [size(grid_x{1}), length(usubs)]),2);
sound_response_mask{2} = sum(reshape(x((1:n_voxels_lh) + n_voxels_rh, :), [size(grid_x{2}), length(usubs)]),2);
ic_project_axis = nan(10,length(usubs),size(Vica,1));
pitchmus_project_axis = nan(10,length(usubs),size(pitchmusB,1));

for k = 1:length(hemis)
    
    if strcmp(hemi,'rh')
        %     bverts = [104804 143214];
        %     bverts = [52628 66233];
        %     bverts = [88767 137574];
        %     bpts = [patch.vras(bi(1),[1 2])' patch.vras(bi(2),[1 2])'];
        bpts = [[20,20]', [40 40]'];
    elseif strcmp(hemi,'lh')
        %     bverts = [38431 26238];
        %     bverts = [11400 60981];
        %     bpts = [patch.vras(bi(1),[1 2])' patch.vras(bi(2),[1 2])'];
        bpts = [[20,20]', [40 40]'];
    end
    
    % project voxels onto vector specified by boundary points
    projection_vector = bpts(:,2) - bpts(:,1);
    projection_vector = projection_vector/norm(projection_vector);
    vras_grid_coordinates = [grid_x{i}(:)'; grid_y{i}(:)'] - bpts(:,1)*ones(1,numel(grid_x{i}));
    voxel_axis_weights = projection_vector'*vras_grid_coordinates;
    
    % sort voxels
    [~,voxel_axis_order] = sort(voxel_axis_weights,'ascend');
    sound_response_CDF = cumsum(sound_response_mask{k}(voxel_axis_order));
    sound_response_CDF = sound_response_CDF/sum(sound_response_CDF);
    
    %     ngroups = 10;
    %     for i = 1:ngroups
    %         xi = sound_response_CDF > ((i-1)*totalarea/ngroups) & sound_response_CDF <= (i*totalarea/ngroups);
    %         lb_new.vnums = sort( verts_sort( inds ) );
    %         [~,x] = intersect(lb.vnums, lb_new.vnums);
    %         lb_new.vras = lb.vras(x,:);
    %         %     write_label(lb_new,['~/freesurfer/fsaverage/label/' hemi '.stp-r' num2str(i) 'of' num2str(nlabels) '.label'],'overwrite');
    %         write_label(lb_new,['~/freesurfer/fsaverage/label/' hemi '.stp-aud-pm2al-r' num2str(i) 'of' num2str(nlabels) '.label'],'overwrite');
    %         annot_labels(lb_new.vnums+1) = colortable.table(i+1,5);
    %     end
    
    for i = 1:size(Vica,1);
        for j = 1:length(usubs)
            x = ic{i,k}(:,:,j);
            ic_axis_sort = x(voxel_axis_order);
            x = pitchmusB_grid{i,k}(:,:,j);
            pitchmus_axis_sort = x(voxel_axis_order);
            for q = 1:ngroups
                xi = sound_response_CDF > ((q-1)/ngroups) & sound_response_CDF <= (q/ngroups);
                ic_project_axis(q,j,i) = sum(ic_axis_sort(xi));
                pitchmus_project_axis(q,j,i) = sum(pitchmus_axis_sort(xi));
            end
        end
    end
end



%
% ica_dimension_bysubject = reshape(ica_dimension_grid, [n_voxels, length(usubs)]);
%
% ica_dimension_grid = cell(1,2);
% ica_dimension_grid{1} =

% figure;
% vnorm_recon = interp2(xi,yi,vgrid,coords(:,1),coords(:,2));
% vnorm_recon(isnan(vnorm_recon)) = 0;
% scatter(coords(:,1),coords(:,2),3,vnorm_recon);
%
%

%%

% addpath(genpath([pwd '/InfoTheory']))
% addpath(genpath([pwd '/ITE-0.56_code']));

%
% [U,S,V] = svd(v,'econ');
% nPCs = 3;
% Vica = V(:,1:nPCs)*sqrt(size(V,1)-1);
% Uica = U(:,1:nPCs)/sqrt(size(V,1)-1);
% Sica = S(1:nPCs,1:nPCs);
%
% ResetRandStream(0);
% close all;
% icasig = Vica';
% [~,~,si_unique] = unique(si);
% for i = 1:nPCs
%   binsize = 100;
%   bins = linspace(-10,10,binsize);
%   subject_hist = nan(binsize,length(usubs));
%   for j = 1:length(usubs)
%     x = hist(icasig(i,si_unique == j),bins);
%     subject_hist(:,j) = x/sum(x);
%   end
%   colormap(['jet(' num2str(length(usubs)) ')']);
%   plot(bins,subject_hist);
%   print(gcf, '-dpdf', [figure_directory 'pca-' num2str(nPCs) '-hist-bysubject-' num2str(i) '.pdf'], '-r100');
% end
% close all;
%
% n = 40;
% th = linspace(-pi/2,pi/2,n);
% phi = linspace(-pi/2,pi/2,n);
% meankurt = nan(n,n);
% minkurt = nan(n,n);
% rotMat = nan(3,3,n,n);
% for i = 1:n
%     for j = 1:n
%         rotMat(:,:,i,j) = rotationmat3D(th(i),[0 1 0])*rotationmat3D(phi(j),[0 0 1]);
%         meankurt(i,j) = mean(abs(kurtosis(Vica*rotMat(:,:,i,j)')-3));
%         minkurt(i,j) = min(abs(kurtosis(Vica*rotMat(:,:,i,j)')-3));
%     end
% end
%
% meankurt_norm = norm01(meankurt);
% minkurt_norm = norm01(minkurt);
%
% figure;
% hold on;
% c = colormap('jet(100)');
% for i = 1:n
%     for j = 1:n
%       vec = rotMat(:,:,i,j)*[1 0 0]';
%       plot(vec(2),vec(3),'o','Color', c(round(meankurt_norm(i,j)*99)+1,:), 'LineWidth',5);
%
% %       plot3(vec(1),vec(2),vec(3),'o','Color', c(round(meankurt_norm(i,j)*99)+1,:), 'LineWidth',5);
% %             %             surf
% % %       if meankurt_norm(i,j) > 0
% %       vec = rotMat(:,:,i,j)*[1 0 0]'*meankurt_norm(i,j);
% %         quiver3(0,0,0,vec(1),vec(2),vec(2),'Color', c(round(meankurt_norm(i,j)*99)+1,:), 'LineWidth',3)
%     end
% end
% export_fig([figure_directory 'ICA-' num2str(nPCs) '_mean_kurtosis_radian_plot.pdf'],'-pdf','-nocrop','-transparent');
%
%
% figure;
% hold on;
% c = colormap('jet(100)');
% for i = 1:n
%     for j = 1:n
%       vec = rotMat(:,:,i,j)*[1 0 0]';
%       plot(vec(2),vec(3),'o','Color', c(round(minkurt_norm(i,j)*99)+1,:), 'LineWidth',5);
%
% %       plot3(vec(1),vec(2),vec(3),'o','Color', c(round(meankurt_norm(i,j)*99)+1,:), 'LineWidth',5);
% %             %             surf
% % %       if meankurt_norm(i,j) > 0
% %       vec = rotMat(:,:,i,j)*[1 0 0]'*meankurt_norm(i,j);
% %         quiver3(0,0,0,vec(1),vec(2),vec(2),'Color', c(round(meankurt_norm(i,j)*99)+1,:), 'LineWidth',3)
%     end
% end
% export_fig([figure_directory 'ICA-' num2str(nPCs) '_min_kurtosis_radian_plot.pdf'],'-pdf','-nocrop','-transparent');
%
% figure;
% imagesc(minkurt_norm);
% export_fig([figure_directory 'ICA-' num2str(nPCs) '_min_kurtosis_cartesian_plot.png'],'-png','-nocrop','-transparent','-r100');
%
% figure;
% imagesc(meankurt_norm);
% export_fig([figure_directory 'ICA-' num2str(nPCs) '_mean_kurtosis_cartesian_plot.png'],'-png','-nocrop','-transparent','-r100');
%
% figure;
% quiver3(zeros(1,3),zeros(1,3),zeros(1,3),rotMat(1,:,10,9),rotMat(2,:,10,9),rotMat(3,:,10,9),'r');
% hold on;
% quiver3(zeros(1,3),zeros(1,3),zeros(1,3),-rotMat(1,:,10,9),-rotMat(2,:,10,9),-rotMat(3,:,10,9),'r');
% quiver3(zeros(1,3),zeros(1,3),zeros(1,3),rotMat(1,:,29,9),rotMat(2,:,29,9),rotMat(3,:,29,9),'b');
% quiver3(zeros(1,3),zeros(1,3),zeros(1,3),-rotMat(1,:,29,9),-rotMat(2,:,29,9),-rotMat(3,:,29,9),'b');
% export_fig([figure_directory 'ICA-' num2str(nPCs) '_best_rotations.pdf'],'-pdf','-nocrop','-transparent');
%
% naturalsound_acoustic_category_figures(Uica*Sica*rotMat(:,:,10,9), figure_directory, ['ica-rot-10-9-' num2str(nPCs) '-']);
% naturalsound_acoustic_category_figures(Uica*Sica*rotMat(:,:,29,9), figure_directory, ['ica-rot-29-9-' num2str(nPCs) '-']);
% naturalsound_acoustic_category_figures(Uica*Sica*rotMat(:,:,8,17), figure_directory, ['ica-rot-8-17-' num2str(nPCs) '-']);
% naturalsound_acoustic_category_figures(Uica*Sica*rotMat(:,:,27,17), figure_directory, ['ica-rot-27-17-' num2str(nPCs) '-']);
%
% % scatter3(Vica*rotMat(1,:,10,9)', Vica*rotMat(2,:,10,9)', Vica*rotMat(3,:,10,9)',1);
%
% best_indices = [10,9;29,9;8,17;27,17];
% for q = 1:size(best_indices)
%   icasig = rotMat(:,:,best_indices(q,1),best_indices(q,2))*Vica';
%   ResetRandStream(0);
%   close all;
%   [~,~,si_unique] = unique(si);
%   for i = 1:nPCs
%     binsize = 100;
%     bins = linspace(-10,10,binsize);
%     subject_hist = nan(binsize,length(usubs));
%     for j = 1:length(usubs)
%       x = hist(icasig(i,si_unique == j),bins);
%       subject_hist(:,j) = x/sum(x);
%     end
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     plot(bins,subject_hist);
%     print(gcf, '-dpdf', [figure_directory 'ica-rot-' num2str(best_indices(q,1)) '-' num2str(best_indices(q,2)) '-' num2str(nPCs) '-hist-bysubject-' num2str(i) '.pdf'], '-r100');
%   end
%   close all;
% end

%% ICA on first N components
%
% [U,S,V] = svd(v,'econ');
% nPCs = 3;
% Vica = V(:,1:nPCs);
% Uica = U(:,1:nPCs);
% Sica = S(1:nPCs,1:nPCs);
% for i = 1:5
%   ResetRandStream(i);
%   [icasig, A, W] = fastica(Vica');
%   [~,xi] = sort(var(Uica*Sica*A),'descend');
%   A = A(:,xi); icasig = icasig(xi,:);
%   % A(:,2:3) = -A(:,2:3);
%   % icasig(2:3,:) = -icasig(2:3,:);
%   figure;
%   exvar = var(Uica*Sica*A);
%   bar(100*exvar,'FaceColor',[1 1 1]);
%   xlabel('ICA Components'); ylabel('Explained Variance');
%   export_fig([figure_directory 'ICA-' num2str(nPCs) '-RandSeed' num2str(i) '_exvar.pdf'],'-pdf','-nocrop');
%   % x = Uica*Sica*A;
%   % response = (x-ones(size(x,1),1)*min(x)) ./ (ones(size(x,1),1)*(max(x)-min(x)));
%   naturalsound_acoustic_category_figures(Uica*Sica*A, figure_directory, ['ica-' num2str(nPCs) '-RandSeed' num2str(i) '-']);
% end

% %% Residual of grandmean projected onto first two pca components
% resid = grandmean - U(:,2)*pinv(U(:,2))*grandmean;
% naturalsound_acoustic_category_figures(resid, figure_directory, 'grandmean-resid-pc2-');


% %% ICA histograms
% icasig = VicaRot';
% ResetRandStream(0);
% close all;
% [~,~,si_unique] = unique(si);
% for i = 1:nPCs
%     binsize = 100;
%     bins = linspace(-10,10,binsize);
%     subject_hist = nan(binsize,length(usubs));
%     for j = 1:length(usubs)
%         x = hist(icasig(i,si_unique == j),bins);
%         subject_hist(:,j) = x/sum(x);
%     end
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     plot(bins,subject_hist);
%     print(gcf, '-dpdf', [figure_directory 'ica-' num2str(nPCs) '-hist-12-9-bysubject-' num2str(i) '.pdf'], '-r100');
% end
% close all;

%%
% mask_threshold = 0.1;
% idstring = ['ica-3-' DataHash(icasig)];
% surfproj_dimensions('naturalsound', usubs, icasig, ai, si, idstring, mask_threshold, component_threshold, figure_directory, varargin{:})
%
% mask_threshold = 0.1;
% component_threshold = 0.1;
% x = Uica*Sica*A;
% icasig_norm = icasig .* (std(x)' * ones(1,size(icasig,2)));
% icasig_norm = icasig_norm / std(icasig_norm(:));
% idstring = ['ica-3-norm-' DataHash(icasig_norm)];
% surfproj_dimensions('naturalsound', usubs, icasig_norm, ai, si, idstring, mask_threshold, component_threshold, figure_directory, varargin{:})

% %% ICA on subject-specific components
% [U,S,V] = svd(v,'econ');
% xi = [1 2 4 5];
% Vica = V(:,xi);
% Uica = U(:,xi);
% Sica = S(xi,xi);
% [icasig, A, W] = fastica(Vica');
% [~,xi] = sort(var(Uica*Sica*A),'descend');
% A = A(:,xi); icasig = icasig(xi,:);
% naturalsound_acoustic_category_figures(Uica*Sica*A, figure_directory, ['ica-1-2-4-5-']);
%
% ResetRandStream(0);
% close all;
% [~,~,si_unique] = unique(si);
% for i = 1:size(icasig,1)
%     binsize = 100;
%     bins = linspace(-10,10,binsize);
%     subject_hist = nan(binsize,length(usubs));
%     for j = 1:length(usubs)
%         x = hist(icasig(i,si_unique == j),bins);
%         subject_hist(:,j) = x/sum(x);
%     end
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     plot(bins,subject_hist);
%     print(gcf, '-dpdf', [figure_directory 'ica-1-2-4-5-hist-bysubject' num2str(i) '.pdf'], '-r100');
% end
% close all;


%%


%
% %% Plot PCA and ICA scatter plot
% [U,S,V] = svd(v,'econ');
% figure;
% set(gcf, 'Position', [0 0 1000 400]);
% xi = Shuffle(1:size(V,1));
% subplot(1,2,1);
% colormap('jet(10)');
% [~,~,si_unique] = unique(si);
% scatter(V(xi,1)*sqrt(size(V,1)),V(xi,2)*sqrt(size(V,1)),1,si_unique(xi));
% xlabel(['PCA' num2str(1)]);
% ylabel(['PCA' num2str(2)]);
% title('PCA');
% xlim([-6 6]); ylim([-6 6]);
%
% [icasig, A, W] = fastica(V(:,1:2)');
% subplot(1,2,2);
% colormap('jet(10)');
% scatter(icasig(1,xi),icasig(2,xi),1,si_unique(xi));
% xlabel(['ICA' num2str(1)]);
% ylabel(['ICA' num2str(2)]);
% title('ICA');
% xlim([-6 6]); ylim([-6 6]);
% export_fig([figure_directory 'pca-ica-scatter-1-2.pdf'], '-pdf','-nocrop');
%
% x = U(:,1:2)*S(1:2,1:2)*A;
% response = (x-ones(size(x,1),1)*min(x)) ./ (ones(size(x,1),1)*(max(x)-min(x)));
% naturalsound_acoustic_category_figures(response, figure_directory, 'ica-2-');
%
% %%
% % [icasig, A, W] = fastica(v,'lasteig');
%
% %%
%
%
% [U,S,V] = svd(v(:,si == 123),'econ');
% % nPCs = 20;
% nPCs = 3;
% Vica = V(:,1:nPCs);
% Uica = U(:,1:nPCs);
% Sica = S(1:nPCs,1:nPCs);
% [icasig, A, W] = fastica(Vica');
% [~,xi] = sort(var(Uica*Sica*A),'descend');
% A = A(:,xi); icasig = icasig(xi,:);
% naturalsound_acoustic_category_figures(Uica*Sica*A, figure_directory, ['ica-' num2str(nPCs) '-us123-']);
%
% % export_fig([figure_directory 'pca-ica-scatter-1-2.pdf'], '-pdf','-nocrop');
%
% %%
% [icasig, A, W] = fastica(v);
% [~,xi] = sort(var(A), 'descend');
% A = A(:,xi); icasig = icasig(xi,:);
% naturalsound_acoustic_category_figures(A(:,1:10), pwd, 'tmp-');
%
%
%
% % naturalsound_acoustic_category_figures(Uica*Sica*A, figure_directory, ['ica-3-4-']);
%
%
% %
% %
% xi = Shuffle(1:size(V,1));
% % subplot(1,2,2);
% colormap('jet(10)');
% scatter(icasig(1,xi),icasig(2,xi),1,si_unique(xi));
% % scatter3(icasig(1,xi),icasig(2,xi),icasig(3,xi),1,si_unique(xi));
% xlabel('X');ylabel('Y');zlabel('Z');
% % xlabel(['ICA' num2str(1)]);
% % ylabel(['ICA' num2str(2)]);
% title('ICA');
% xlim([-6 6]); ylim([-6 6]);
% % export_fig([figure_directory 'pca-ica-scatter-1-2.pdf'], '-pdf','-nocrop');
%
% % naturalsound_acoustic_category_figures(Uica*A, figure_directory, ['ica-' num2str(nPCs) '-v2-']);
%
% %% ICA on subject-specific components
% [U,S,V] = svd(v,'econ');
% xi = 3:12;
% Vica = V(:,xi);
% Uica = U(:,xi);
% Sica = S(xi,xi);
% [icasig, A, W] = fastica(Vica');
%
% ResetRandStream(0);
% close all;
% [~,~,si_unique] = unique(si);
% for i = 1:size(icasig,1)
%     binsize = 100;
%     bins = linspace(-10,10,binsize);
%     subject_hist = nan(binsize,length(usubs));
%     for j = 1:length(usubs)
%         x = hist(icasig(i,si_unique == j),bins);
%         subject_hist(:,j) = x/sum(x);
%     end
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     plot(bins,subject_hist);
%     print(gcf, '-dpdf', [figure_directory 'ica-3-12-hist-bysubject' num2str(i) '.pdf'], '-r100');
% end
% close all;
%
% xi = Shuffle(1:size(icasig,2));
% for i = 1:2:9
%     figure;
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     scatter(icasig(i,xi),icasig(i+1,xi),ones(size(V,1),1)*3,si_unique(xi));
%     xlabel(['PCA' num2str(i)]);
%     ylabel(['PCA' num2str(i+1)]);
%     h = colorbar; set(h,'YTick',1:length(usubs));
%     print(gcf, '-dpdf', [figure_directory 'ica-scatter-3-12-subjectlabeled-' num2str(i) '-' num2str(i+1) '.pdf'], '-r100');
% end
% close all;
%
% %% ICA on general components
% [U,S,V] = svd(v,'econ');
% xi = [3:4];
% Vica = V(:,xi);
% Uica = U(:,xi);
% Sica = S(xi,xi);
% [icasig, A, W] = fastica(Vica');
% [~,xi] = sort(var(Uica*Sica*A),'descend');
% A = A(:,xi); icasig = icasig(xi,:);
% figure;
% exvar = var(Uica*Sica*A);
% plot(100*exvar,'o-','LineWidth',1);
% xlabel('ICA Components'); ylabel('% Explained Variance');
% export_fig([figure_directory 'ICA-3-4_eigspec.pdf'],'-pdf','-nocrop');
% naturalsound_acoustic_category_figures(Uica*Sica*A, figure_directory, ['ica-3-4-']);
%
% ResetRandStream(0);
% close all;
% [~,~,si_unique] = unique(si);
% for i = 1:size(icasig,1)
%     binsize = 100;
%     bins = linspace(-10,10,binsize);
%     subject_hist = nan(binsize,length(usubs));
%     for j = 1:length(usubs)
%         x = hist(icasig(i,si_unique == j),bins);
%         subject_hist(:,j) = x/sum(x);
%     end
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     plot(bins,subject_hist);
%     print(gcf, '-dpdf', [figure_directory 'ica-1-2-13-20-hist-bysubject' num2str(i) '.pdf'], '-r100');
% end
% close all;
%
% xi = Shuffle(1:size(icasig,2));
% for i = 1:2:9
%     figure;
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     scatter(icasig(i,xi),icasig(i+1,xi),ones(size(V,1),1)*3,si_unique(xi));
%     xlabel(['PCA' num2str(i)]);
%     ylabel(['PCA' num2str(i+1)]);
%     h = colorbar; set(h,'YTick',1:length(usubs));
%     print(gcf, '-dpdf', [figure_directory 'ica-scatter-1-2-13-20-subjectlabeled-' num2str(i) '-' num2str(i+1) '.pdf'], '-r100');
% end
% close all;
%
% %% PCA eigen value spectrum for two splits of data
% % 121   123    11   160   143   122     1   164   106   161
% % 76.3319   73.0574   72.8342   72.6270   71.3884   70.6876   70.1991   69.2998   68.1907   63.5919
% % usubs_splithalf = {[121, 123, 11, 160, 143], [122, 1, 164, 106, 161]};
% usubs_splithalf = {[121, 11, 143, 1, 106], [123, 160, 122, 164, 161]};
% for i = 1:2
%     [U,S,V] = svd(v(:,ismember(si, usubs_splithalf{i})),'econ');
%     eigs = diag(S).^2;
%     exvar = eigs/sum(eigs);
%     plot(100*exvar,'o-','LineWidth',1);
%     xlabel('PCA Components'); ylabel('% Explained Variance');
%     export_fig([figure_directory 'PCA_eigspec_us' sprintf('%d',usubs_splithalf{i}) '.pdf'],'-pdf','-nocrop');
%     naturalsound_acoustic_category_figures(U(:,1:5), figure_directory, ['pca-us' sprintf('%d',usubs_splithalf{i}) '-']);
% end
%
% %% Reliability of PC's in left out data
% [U,S,V] = svd(v(:,ismember(si, usubs_splithalf{1})),'econ');
%
% data1 = v1(:,ismember(si, usubs_splithalf{2}));
% data2 = v2(:,ismember(si, usubs_splithalf{2}));
% V1hat = U'*data1;
% V2hat = U'*data2;
% corr = nan(165,1);
% for j = 1:165
%     p1 = U(:,j)*V1hat(j,:);
%     p2 = U(:,j)*V2hat(j,:);
%     x1 = fastcorr3(p1(:),data2(:));
%     x2 = fastcorr3(p2(:),data1(:));
%     corr(j) = (x1 + x2)/2;
% end
%
% [explained_variance_sort,xi] = sort(sign(corr(1:50)).*corr(1:50).^2,'descend');
% Usort = U(:,xi);
% bar(explained_variance_sort,'FaceColor',[1 1 1]);
% xlabel('PCs');ylabel('Explained Variance in Data Left Out Subjects');
% export_fig([figure_directory 'PCA_splithalf_explained_variance_sort'  '.pdf'],'-pdf','-nocrop');
% naturalsound_acoustic_category_figures(Usort(:,1:5), figure_directory, ['pca-us' sprintf('%d',usubs_splithalf{1}) '-reliability-sorted']);


%%
% for i = 1:2
%     [U,S,V] = svd(v,'econ');
%     figure;
%     set(gcf, 'Position', [0 0 1000 400]);
%     xi = Shuffle(1:size(V,1));
%     subplot(1,2,1);
%     colormap('jet(10)');
%     [~,~,si_unique] = unique(si);
%     scatter(V(xi,1)*sqrt(size(V,1)),V(xi,2)*sqrt(size(V,1)),1,si_unique(xi));
%     xlabel(['PCA' num2str(1)]);
%     ylabel(['PCA' num2str(2)]);
%     title('PCA');
%     xlim([-6 6]); ylim([-6 6]);
%
%     [icasig, A, W] = fastica(V(:,1:2)');
%     subplot(1,2,2);
%     colormap('jet(10)');
%     scatter(icasig(1,xi),icasig(2,xi),1,si_unique(xi));
%     xlabel(['ICA' num2str(1)]);
%     ylabel(['ICA' num2str(2)]);
%     title('ICA');
%     xlim([-6 6]); ylim([-6 6]);
%     export_fig([figure_directory 'pca-ica-scatter-1-2.pdf'], '-pdf','-nocrop');
% end


%% Subject-specific clustering index for PCA components
% nPCs = 20;
% matfile = [figure_directory 'subject_specificity' num2str(nPCs) '.mat'];
% if ~exist(matfile, 'file')
%     [U,S,V] = svd(v,'econ');
%     subject_specificity = nan(1,nPCs);
%     for i = 1:nPCs
%         i
%         drawnow;
%         tic;
%         d = pdist(V(:,i));
%         subject_nonmatch = (pdist(si) > 0);
%         subject_specificity(i) = mean(d(~subject_nonmatch))/mean(d(subject_nonmatch));
%         toc;
%     end
%     save(matfile,'subject_specificity');
% else
%     load(matfile);
% end
% plot(1./subject_specificity,'k-o');
% export_fig([figure_directory 'subject_specificity' num2str(nPCs) '.pdf'],'-pdf','-nocrop');

%% 2D PCA plots showing voxels colored by subjects
% [U,S,V] = svd(v,'econ');
% xi = Shuffle(1:size(V,1));
% [~,~,unique_si] = unique(si);
% nPCs = 20;
% for i = 1:nPCs
%     figure;
%     colormap(['jet(' num2str(length(usubs)) ')']);
%     scatter(S(i,i)*V(xi,i),S(i+1,i+1)*V(xi,i+1),ones(size(V,1),1)*3,unique_si(xi));
%     xlabel(['PCA' num2str(i)]);
%     ylabel(['PCA' num2str(i+1)]);
%     h = colorbar; set(h,'YTick',1:length(usubs));
%     print(gcf, '-dpdf', [figure_directory 'pca-scatter-subjectlabeled-' num2str(i) '-' num2str(i+1) '.pdf'], '-r100');
% end
% close all;


%% Subject-specificity of PCs based on half the data


% %%
%
% % subject_specificity = nan(1,165);
% % for i = [47:165]
% %     i
% %     drawnow;
% %     tic;
% %     d = pdist(V(:,i));
% %     subject_nonmatch = (pdist(si) > 0);
% %     subject_specificity(i) = mean(d(~subject_nonmatch))/mean(d(subject_nonmatch));
% %     toc;
% % end
% % save([figure_directory 'subject_specificity.mat'],'subject_specificity');
%
% % load([figure_directory 'subject_specificity.mat']);
% % plot(1./subject_specificity(1:50),'k-o');
% % export_fig([figure_directory 'subject_specificity.pdf'],'-pdf','-nocrop');
% %
% % subject_specific_variance = nan(length(usubs),165);
% % for i = 1:length(usubs)
% %     subject_specific_variance(i,:) = mean(V(si == usubs(i),:).^2,1);
% % end
% % skew = sum(subject_specific_variance.^3);
%
% %% PCA
% [U,S,V] = svd(v,'econ');
% eigs = diag(S).^2;
% exvar = eigs/sum(eigs);
% [Ushuf,Sshuf,Vshuf] = svd(v(xi),'econ');
% eigshuf = diag(Sshuf).^2;
% exvar_shuf = eigshuf/sum(eigshuf);
% figure;
% plot(100*exvar(1:50),'ko-','LineWidth',2);
% hold on;
% plot(100*exvar_shuf(1:50),'ro-','LineWidth',1);
% % errorbar(1:25,100*mean(exvar_shuf(1:25,:)'),100*std(exvar_shuf(1:25,:)'),'LineWidth',1);
% xlabel('PCA Components'); ylabel('% Explained Variance');
% export_fig([figure_directory 'PCA_eigspec.pdf'],'-pdf','-nocrop');
%
% %% ICA
%
% % k = 11;
% icvar = var(A);
% [~,xi] = sort(icvar,'descend');
% A = A(:,xi); icasig = icasig(xi,:); icvar = icvar(xi);
% figure;
% plot(100*icvar,'o-','LineWidth',1);
% xlabel('ICA Components'); ylabel('% Explained Variance');
% export_fig([figure_directory 'ICA_eigspec.pdf'],'-pdf','-nocrop');
% naturalsound_acoustic_category_figures(A, figure_directory, 'ica-11-');
%
%
%
%
%
% %%
% [U,S,V] = svd(v1,'econ');
% projection1 = U'*v1;
% projection2 = U'*v2;
% cross_scan_corr = fastcorr3(projection1',projection2');
% plot(cross_scan_corr);
% %%
%
%
%
% eigs = diag(S).^2;
% exvar = eigs/sum(eigs);
% plot(100*exvar,'o-','LineWidth',1);
% xlabel('PCA Components'); ylabel('% Explained Variance');
% export_fig([figure_directory 'PCA_eigspec.pdf'],'-pdf','-nocrop');
%
%
% %%
%
% [icasig, A, W] = fastica(v,'approach','symm');
% eigs = var(A); eigs = eigs/sum(eigs);
% [~,xi] = sort(eigs,'descend');
% A = A(:,xi); icasig = icasig(xi,:); eigs = eigs(xi);
% figure;
% plot(100*eigs,'o-','LineWidth',1);
% xlabel('ICA Components'); ylabel('% Explained Variance');
% export_fig([figure_directory 'ICA_eigspec.pdf'],'-pdf','-nocrop');
% naturalsound_acoustic_category_figures(A(:,1:20), figure_directory, 'ica-165-');
%
% for k = 20
%     ResetRandStream(0);
%     close all;
%     [~,~,si_unique] = unique(si);
%     for i = 1:k
%         binsize = 100;
%         bins = linspace(-10,10,binsize);
%         subject_hist = nan(binsize,length(usubs));
%         for j = 1:length(usubs)
%             x = hist(icasig(i,si_unique == j),bins);
%             subject_hist(:,j) = x/sum(x);
%         end
%         colormap(['jet(' num2str(length(usubs)) ')']);
%         plot(bins,subject_hist);
%         print(gcf, '-dpdf', [figure_directory 'ica-165-hist-bysubject-c' num2str(i) '.pdf'], '-r100');
%     end
% end
%

% exvar_shuf = nan(165,10);
% for i = 1:10
%     tic;
%     xi = Shuffle((1:165)'*ones(1,size(v,2)) + ones(165,1)*(0:size(v,2)-1)*165);
%     [Ushuf,Sshuf,Vshuf] = svd(v(xi),'econ');
%     eigshuf = diag(Sshuf).^2;
%     exvar_shuf(:,i) = eigshuf/sum(eigshuf);
%     toc;
% end

%      = fastcorr3(projection1',projection2').^2*(var(projection1') + var(projection2'))';

%     Ua(:,:,i) = U;

%     subplot(1,2,1);
%     matchcorr1 = corr(U_loo,U);
%     imagesc(matchcorr1(1:165,1:165))
% match eigenvectors
%     matchcorr1 = corr(U_loo,U);
%     [matching,cost] = Hungarian(1-matchcorr1);
%     x = (1:165)'*ones(1,165);
%     xi = x(logical(matching));



%     x = (1:165)'*ones(1,165);
%     xi = x(logical(matching));
%     matchcorr2 = corr(U,U_loo(:,xi));
%     subplot(1,2,2);
%     imagesc(matchcorr2(1:30,1:30));

%% Explained variance of ICA components left-out subjects
% addpath(genpath([pwd '/functional_systems/']))
% % [U,S,V] = svd(v,'econ');
% % U_all = nan(165,165,10);
% r = nan(165,10);
% for i = 1:length(usubs)
%     i
%     drawnow;
%     %     [U_loo,S_loo,~] = svd(v(:,si ~= usubs(i)),'econ');
%     %     xhat = U_loo'*v1(:,si == usubs(i));
%     for j = 1:20
%
%         prediction = U_loo(:,1:j)*S_loo(1:j,1:j)*xhat(1:j,:);
%         scan2 = v2(:,si == usubs(i));
%         r(j,i) = fastcorr3(prediction(:), scan2(:)).^2;
%     end
% end
% x = diff([zeros(1,10);r(1:165,:)]);
% bar(mean(x(1:20,:),2),'FaceColor',[1 1 1]);
% hold on;
% errorbar(1:20,mean(x(1:20,:)'),stderr_withsub(x(1:20,:)'),'k.');
% xlabel('PCs');ylabel('Explained Variance in Data from Leftout Subject');
% export_fig([figure_directory 'PCA_eigspec_loso.pdf'],'-pdf','-nocrop');

% %% Explained variance in individual subjects
%
% r = nan(165,10);
% for i = 1:length(usubs)
%     xhat = U'*v1(:,si == usubs(i));
%     for j = 1:165
%         prediction = U(:,j)*xhat(j,:);
%         scan2 = v2(:,si == usubs(i));
%         r(j,i) = fastcorr3(prediction(:), scan2(:));
%     end
% end
% % x = diff([zeros(1,10);r(1:165,:)]);
% explained_variance_mean = sign(mean(r')).*mean(r').^2;
% explained_variance_sem = sign(mean(r')).*mean(r');
% bar(mean(),'FaceColor',[1 1 1]);
% hold on;
% errorbar(1:20,mean(x(1:20,:)'),stderr_withsub(x(1:20,:)'),'k.');
% xlabel('PCs');ylabel('Explained Variance in Individual Subjects');
% export_fig([figure_directory 'PCA_eigspec_loso.pdf'],'-pdf','-nocrop');

%%
% subjid = ['naturalsound_us1'];
% blocktp = 3:6;
%
%
% % right
% i = 1;
% hemis = {'rh','lh'};
%
% % contraint roi
% mask{i} = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.audctx.label']);
%
% % patch coordinates and surface voxel indices
% patch{i} = read_patch(['~/freesurfer/fsaverage/surf/' hemis{i} '.cortex.patch.flat']);
% [~,ai,ib] = intersect(abs(patch{i}.vnums),mask{i}.vnums);
% vras{i} = patch{i}.vras(ai,1:2);
% vi{i} = mask{i}.vnums(ib);% voxel responses
%
% % voxel responses
% x = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r1/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.mgz']);
% vr{i} = squeeze(x.vol(:,vi{i}+1,:,:))';
%
% % bounds for regridding
% bounds{i} = [floor(min(vras{i})); ceil(max(vras{i}))];
% [grid_x{i},grid_y{i}] = meshgrid(bounds{i}(1,1):2:bounds{i}(2,1),(bounds{i}(1,2):2:bounds{i}(2,2)));
% vgrid{i} = nan(size(grid_x{i},1),size(grid_y{i},2), size(vr{i},1));
% for j = 1:165
%     vgrid{i}(:,:,j) = griddata(vras{i}(:,1),vras{i}(:,2),vr{i}(j,:),grid_x{i},grid_y{i});
% end
%
% for j = 1:2
%     figure;
%     scatter(vras{i}(:,1),vras{i}(:,2),3,vr{i}(j,:));
%     figure;
%     %     surf(grid_x{i},grid_y{i},vgrid{i}(:,:,j));
%     imagesc(flipud(vgrid{i}(:,:,j)));
% end

% % left
% mask_lh = read_label([params('rootdir') 'freesurfer/fsaverage/label/lh.audctx.label']);
% patch_lh = read_patch('~/freesurfer/fsaverage/surf/lh.cortex2.patch.flat');
% [~,lia,lib] = intersect(abs(patch_lh.vnums),mask_lh.vnums);
% vras_lh = patch_lh.vras(lia,1:2);
% vn_lh = mask_lh.vnums(lib);
% x = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r1/lh.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.mgz']);
% v_lh = squeeze(x.vol(:,vn_lh+1,:,:))';
%
% % combined
% % vras = [vras_rh; vras_lh];
% % vn = [vn_rh; vn_lh];
% % hemi = [2*ones(length(vn_rh),1); 1*ones(length(vn_lh),1)];
%
%
% vgrid = griddata(vras_lh(:,1),vras_lh(:,2),v_lh(1,:)',xi_lh,xi_lh);
% vgrid(isnan(vgrid)) = 0;
%
% bounds_rh = [floor(min(vras_rh)); ceil(max(vras_rh))];
%
%
% %
%
% figure;
% scatter(coords(:,1),coords(:,2),3,vnorm(1,:));
%
% bounds = [floor(min(coords)); ceil(max(coords))];
% % bounds = [[50 62]',[-10 2]'];
% [xi,yi] = meshgrid(bounds(1,1):2:bounds(2,1),(bounds(1,2):2:bounds(2,2)));
% vgrid = griddata(coords(:,1),coords(:,2),vnorm(1,:),xi,yi);
% vgrid(isnan(vgrid)) = 0;
% figure;
% surf(xi,yi,vgrid);
%
% figure
% acf = xcorr2(vgrid);
% imagesc(acf/max(acf(:)),[0 1]); colorbar;
%
% figure;
% vnorm_recon = interp2(xi,yi,vgrid,coords(:,1),coords(:,2));
% vnorm_recon(isnan(vnorm_recon)) = 0;
% scatter(coords(:,1),coords(:,2),3,vnorm_recon);
%
%
%
% patch_lh = read_patch(['~/freesurfer/fsaverage/surf/lh.cortex2.patch.flat']);
%
% % s1_rh = v_unfilt_s12(:,1:length(mask_rh.vnums));
% % s1_rh = s1_rh ./ (ones(size(s1_rh,1),1)*sqrt(mean(s1_rh.^2,1)));
%
% % xgrid = patch_rh.vras(xi,1)*ones(1,sum(xi));
% % ygrid = ones(sum(xi),1)*patch_rh.vras(xi,2)';
% [~,ia,ib] = intersect(patch_rh.vnums, mask_rh.vnums);
% x = s1_rh(29,ib);
% % x = mean_response(ib);
% x = (x-min(x))/(max(x)-min(x));
%
% grid_response = griddata(patch_rh.vras(ia,1), patch_rh.vras(ia,2),x,xi,yi);
% grid_response(isnan(grid_response)) = 0;
%
% acf = xcorr2(grid_response);
%
% meshgrid(-3:.1:3,-3:.1:3);

% patch_lh = read_patch(['~/freesurfer/fsaverage/surf/lh.cortex2.patch.flat']);

% surf(patch_rh.vras(xi,1),patch_rh.vras(xi,2),mean_response(1:sum(xi)));
%
% x = reshape(v_unfilt_s12',[nvertices_persubject,10,165]);
% mean_response = mean(mean(x,3),2);
% mean_response = mean_response(1:length(mask_rh.vnums));
%% right
% mask_rh = read_label([params('rootdir') 'freesurfer/fsaverage/label/rh.audctx.label']);
% patch_rh = read_patch('~/freesurfer/fsaverage/surf/rh.cortex.patch.flat');
% [~,ria,ib] = intersect(abs(patch_rh.vnums),mask_rh.vnums);
% vras_rh = patch_rh.vras(ria,1:2);
% vn_rh = mask_rh.vnums(ib);
%
% % left
% mask_lh = read_label([params('rootdir') 'freesurfer/fsaverage/label/lh.audctx.label']);
% patch_lh = read_patch('~/freesurfer/fsaverage/surf/lh.cortex2.patch.flat');
% [~,lia,lib] = intersect(abs(patch_lh.vnums),mask_lh.vnums);
% vras_lh = patch_lh.vras(lia,1:2);
% vn_lh = mask_lh.vnums(lib);
%
% % combined
% vras = [vras_rh; vras_lh];
% rl_vn = [vn_rh; vn_lh];
%
% v = squeeze(v_rh.vol(:,vi+1,:,:))';
% vnorm = v ./ (ones(165,1)*sqrt(sum(v.^2)));
% nvertices_persubject = sum(length(mask_rh.vnums)+length(mask_lh.vnums));
% v_unfilt = nan(nvertices_persubject, length(usubs), 165, length(runs));


%% OLD ICA Analyses
%
% % ICA
% % approaches = {'symm','defl'};
% approaches = {'symm'};
% % nonlinearities = {'pow3','tanh','gauss'};
% nonlinearities = {'pow3'};
% sumKurt = nan(length(approaches), length(nonlinearities),1);
% sumNegEntropy = nan(length(approaches), length(nonlinearities),1);
% for j = 1:length(approaches)
%     for k = 1:length(nonlinearities)
%         for i = 1:1
%             ica_directory = [figure_directory 'ica-' num2str(nPCs) '-' approaches{j} '-' nonlinearities{k} '/'];
%             if ~exist(ica_directory,'dir')
%                 mkdir(ica_directory);
%             end
%             ResetRandStream(i);
%             [icasig, A, W] = fastica(Vr,'approach',approaches{j},'g',nonlinearities{k},'epsilon',1e-10,'maxNumIterations',1e6);
%             %             nICs = 4;
%             nICs = nPCs;
%             xi = 1:nICs;
%             xi = xi(sign(skewness(icasig'))==-1);
%             A(:,xi) = -A(:,xi);
%             icasig(xi,:) = -icasig(xi,:);
%             rotMatICA = pinv(A);
%             Uica = Ur*Sr*A(:,1:nICs);
%             Vica = icasig;
%             subplot(1,3,1);
%             bar(var(Uica),'FaceColor',[1 1 1]);
%             xlabel('ICA Components'); ylabel('Explained Variance');
%             xlim([0 6]);
%             kurtICA = abs(kurtosis(Vica')-3);
%             subplot(1,3,2);
%             bar(kurtICA,'FaceColor',[1 1 1]);
%             xlabel('ICA Components'); ylabel('Kurtosis');
%             sumKurt(j,k,i) = sum(kurtICA);
%             xlim([0 6]);
%             %             fprintf('%s, %s: 0.3f',approaches{j},nonlinearities{k},sum(kurtICA));
%             subplot(1,3,3)
%             h = log(sqrt(2*pi*exp(1))) - entropy_SNH(Vica',10,50);
%             sumNegEntropy(j,k,i) = sum(h);
%             bar(h,'FaceColor',[1 1 1]);
%             xlabel('ICA Components'); ylabel('Negentropy (Nats)');
%             xlim([0 6]);
%             export_fig([ica_directory 'RandSeed' num2str(i) '-component-' num2str(i) '_exvar_kurtosis_entropy.pdf'],'-pdf','-nocrop');
%
%             % x = Ur*Sr*A;
%             % response = (x-ones(size(x,1),1)*min(x)) ./ (ones(size(x,1),1)*(max(x)-min(x)));
%             %             close all;
%             %             naturalsound_acoustic_category_figures(Uica, ica_directory, ['RandSeed' num2str(i) '-component-' num2str(i) '-']);
%         end
%     end
% end
% close all;
%
% % Grid search
%
% % negentropy
% figure;
% set(gca,'FontSize',20);
% pairs = flipud(combnk(1:nPCs,2));
% th = linspace(0,2*pi,41*4);
% h_pairs = nan(length(th), size(pairs,1));
% pair_labels = cell(1,size(pairs,1));
% handles = nan(1,size(pairs,1));
% c = colormap(['jet(' num2str(size(pairs,1)) ')']);
% x = {'b','g','r','c','k','m'};
% hold on;
% for i = 1:size(pairs,1)
%     for j = 1:length(th)
%         rotMat = [cos(th(j)), -sin(th(j)); sin(th(j)), cos(th(j))];
%         Vrot = rotMat*Vr(pairs(i,:),:);
%         h_pairs(j,i) = sum(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vrot',14,70));
%     end
%     pair_labels{i} = [num2str(pairs(i,1)) '-' num2str(pairs(i,2))];
%     handles(i) = polar(th', h_pairs(:,i),'r');% plot(th/pi,h_pairs(:,i),'LineWidth',2,'Color',c(i,:));
% end
% % legend(handles,pair_labels,'Location','Best');
% % xlabel('Independence (Negentropy)');ylabel('Independence (Negentropy)');
% ylim([-0.25 0.25]); xlim([-0.25 0.25]);
% set(gca, 'XTick', [-0.2 0 0.2], 'YTick', [-0.2 0 0.2]);
% box off;
% export_fig([ica_directory 'pca-' num2str(nPCs) '-NegEntropyAllPairs.pdf'],'-pdf','-nocrop','-transparent');
%
% % negentropy polar
% figure;
% set(gca,'FontSize',20);
% pairs = flipud(combnk(1:nPCs,2));
% th = linspace(0,2*pi,41*4);
% h_pairs = nan(length(th), size(pairs,1));
% pair_labels = cell(1,size(pairs,1));
% handles = nan(1,size(pairs,1));
% c = colormap(['jet(' num2str(size(pairs,1)) ')']);
% x = {'b','g','r','c','k','m'};
% hold on;
% for i = 1:size(pairs,1)
%     for j = 1:length(th)
%         rotMat = [cos(th(j)), -sin(th(j)); sin(th(j)), cos(th(j))];
%         Vrot = rotMat*Vr(pairs(i,:),:);
%         h_pairs(j,i) = sum(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vrot',14,70));
%     end
%     pair_labels{i} = [num2str(pairs(i,1)) '-' num2str(pairs(i,2))];
%     handles(i) = polar(th', h_pairs(:,i),'r');% plot(th/pi,h_pairs(:,i),'LineWidth',2,'Color',c(i,:));
% end
% % legend(handles,pair_labels,'Location','Best');
% % xlabel('Independence (Negentropy)');ylabel('Independence (Negentropy)');
% ylim([-0.25 0.25]); xlim([-0.25 0.25]);
% set(gca, 'XTick', [-0.2 0 0.2], 'YTick', [-0.2 0 0.2]);
% box off;
% export_fig([ica_directory 'pca-' num2str(nPCs) '-NegEntropyAllPairs-Polar.pdf'],'-pdf','-nocrop','-transparent');
%
%
% % histogram by subject for pca
% plothist_bysubject(Vica, si, [ica_directory 'hist-bysubject']);
%
% % subjects
% [~,~,unique_subject_id] = unique(si);
% pairs = flipud(combnk(1:nPCs,2));
% xi = Shuffle(1:length(si));
% for i = 1:size(pairs,1)
%     colormap('jet(10)');
%     scatter(Vica(pairs(i,1),xi),Vica(pairs(i,2),xi),1,unique_subject_id(xi));
%     xlim([-8 8]); ylim([-8 8]);
%     export_fig([ica_directory 'RandSeed1-scatter-component-' num2str(pairs(i,1)) '-' num2str(pairs(i,2)) '.pdf'],'-pdf','-nocrop');
% end
%
% % negentropy
% figure;
% th = linspace(0,2*pi,41);
% h_pairs = nan(length(th), size(pairs,1));
% pair_labels = cell(1,size(pairs,1));
% handles = nan(1,size(pairs,1));
% c = colormap(['jet(' num2str(size(pairs,1)) ')']);
% hold on;
% for i = 1:size(pairs,1)
%     for j = 1:length(th)
%         rotMat = [cos(th(j)), -sin(th(j)); sin(th(j)), cos(th(j))];
%         Vrot = rotMat*Vica(pairs(i,:),:);
%         h_pairs(j,i) = sum(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vrot',10,50));
%     end
%     pair_labels{i} = [num2str(pairs(i,1)) '-' num2str(pairs(i,2))];
%     handles(i) = polar(th',h_pairs(:,i));%,'LineWidth',2,'Color',c(i,:));
% end
% legend(handles,pair_labels,'Location','NorthEastOutside');
% xlabel('Radians');ylabel('NegEntropy');
% % ylim([0 0.25])
% export_fig([ica_directory 'ica-' num2str(nPCs) '-NegEntropyAllPairs.pdf'],'-pdf','-nocrop');
%
% VicaNew = Vica;
% th = linspace(-pi/4,pi/4,41);
% n_iter = 10;
% negEnt = nan(1,n_iter+1);
% negEnt(1) = mean(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vica',10,50));
% for q = 1:n_iter
%     for i = 1:size(pairs,1)
%         for j = 1:length(th)
%             rotMat = [cos(th(j)), -sin(th(j)); sin(th(j)), cos(th(j))];
%             Vrot = rotMat*VicaNew(pairs(i,:),:);
%             h_pairs(j,i) = mean(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vrot',10,50));
%         end
%         [~,xi] = max(h_pairs(:,i));
%         rotMat = [cos(th(xi)), -sin(th(xi)); sin(th(xi)), cos(th(xi))];
%         VicaNew(pairs(i,:),:) = rotMat*VicaNew(pairs(i,:),:);
%     end
%     negEnt(q+1) = mean(log(sqrt(2*pi*exp(1))) - entropy_SNH(VicaNew',10,50));
% end
%
% figure;
% set(gca,'FontSize',20);
% th = linspace(0,pi*2,41*4);
% h_pairs = nan(length(th), size(pairs,1));
% pair_labels = cell(1,size(pairs,1));
% handles = nan(1,size(pairs,1));
% c = colormap(['jet(' num2str(size(pairs,1)) ')']);
% x = {'b','r','k','m'};
% for i = 1:size(pairs,1)
%     for j = 1:length(th)
%         rotMat = [cos(th(j)), -sin(th(j)); sin(th(j)), cos(th(j))];
%         Vrot = rotMat*VicaNew(pairs(i,:),:);
%         h_pairs(j,i) = sum(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vrot',10,50));
%     end
%     pair_labels{i} = [num2str(pairs(i,1)) '-' num2str(pairs(i,2))];
%     handles(i) = polar(th',h_pairs(:,i),x{mod(i-1,length(x))+1});%plot(th/pi,h_pairs(:,i),'LineWidth',2,'Color',c(i,:));
%     if i == 1;
%         hold on;
%     end
% end
%
% % legend(handles,pair_labels,'Location','NorthEastOutside');
% % xlabel('Radians');ylabel('NegEntropy');
% % set(gca,'XTick',[0 0.25 0.5],'YTick',[0.1 0.2]);
% ylim([-0.25 0.25]); xlim([-0.25 0.25])
% % box off;
% export_fig([ica_directory 'ica-' num2str(nPCs) '-NegEntropyAllPairs_Rotated.pdf'],'-pdf','-nocrop','-transparent');
%
% figure;
% nE_PCA = sum(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vr',10,50));
% nE_ICA = sum(log(sqrt(2*pi*exp(1))) - entropy_SNH(VicaNew',10,50));
%
% bar([nE_PCA, nE_ICA],'FaceColor',[1 1 1]);
%
% %
% Ri = Vica*pinv(VicaNew);
% UicaNew = Ur*Sr*A*Ri;
% subplot(1,3,1);
% bar(var(Uica),'FaceColor',[1 1 1]);
% xlabel('ICA Components'); ylabel('Explained Variance');
% xlim([0 6]);
% kurtICA = abs(kurtosis(VicaNew')-3);
% subplot(1,3,2);
% bar(kurtICA,'FaceColor',[1 1 1]);
% xlabel('ICA Components'); ylabel('Kurtosis');
% xlim([0 6]);
% subplot(1,3,3)
% h = log(sqrt(2*pi*exp(1))) - entropy_SNH(VicaNew',10,50);
% bar(h,'FaceColor',[1 1 1]);
% xlabel('ICA Components'); ylabel('Negentropy (Nats)');
% xlim([0 6]);
% export_fig([ica_directory 'RandSeed' num2str(i) '-Rotated-component-' num2str(i) '_exvar_kurtosis_entropy.eps'],'-eps','-nocrop','-transparent');
%
% UicaNewNorm = nan(size(UicaNew));
% for i = 1:nPCs
%     UicaNewNorm(:,i) = (UicaNew(:,i) - min(UicaNew(:,i)))/(max(UicaNew(:,i))-min(UicaNew(:,i)));
% end
%
%
% % export_fig([ica_directory 'RandSeed' num2str(i) '-Rotated-component-' num2str(i) '_exvar_kurtosis_entropy.pdf'],'-pdf','-nocrop');
% naturalsound_acoustic_category_figures(UicaNew, ica_directory, ['RandSeed' num2str(i) '-Rotated-component-' num2str(i) '-']);
%
% naturalsound_acoustic_category_figures(UicaNewNorm, ica_directory, ['RandSeed' num2str(i) '-MinMax-Rotated-component-' num2str(i) '-']);
%
% plothist_bysubject(VicaNew, si, [ica_directory 'RandSeed' num2str(i) '-Rotated-hist-bysubject']);
%
% % subjects
% [~,~,unique_subject_id] = unique(si);
% pairs = flipud(combnk(1:nPCs,2));
% xi = Shuffle(1:length(si));
% for i = 1:size(pairs,1)
%     colormap('jet(10)');
%     scatter(VicaNew(pairs(i,1),xi),VicaNew(pairs(i,2),xi),1,unique_subject_id(xi));
%     xlim([-8 8]); ylim([-8 8]);
%     export_fig([ica_directory 'RandSeed1-Rotated-scatter-component-' num2str(pairs(i,1)) '-' num2str(pairs(i,2)) '.pdf'],'-pdf','-nocrop');
% end
%
%
% %% Full Neg-entropy based optimization
%
% VicaNegE = Vr;
% th = linspace(-pi/4,pi/4,41);
% n_iter = 10;
% negEnt = nan(1,n_iter+1);
% negEnt(1) = mean(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vr',10,50));
% for q = 1:n_iter
%     h_pairs = nan(length(th),size(pairs,1));
%     for i = 1:size(pairs,1)
%         for j = 1:length(th)
%             rotMat = [cos(th(j)), -sin(th(j)); sin(th(j)), cos(th(j))];
%             Vrot = rotMat*VicaNegE(pairs(i,:),:);
%             h_pairs(j,i) = mean(log(sqrt(2*pi*exp(1))) - entropy_SNH(Vrot',10,50));
%         end
%         [~,xi] = max(h_pairs(:,i));
%         rotMat = [cos(th(xi)), -sin(th(xi)); sin(th(xi)), cos(th(xi))];
%         VicaNegE(pairs(i,:),:) = rotMat*VicaNegE(pairs(i,:),:);
%     end
%     negEnt(q+1) = mean(log(sqrt(2*pi*exp(1))) - entropy_SNH(VicaNegE',10,50));
% end
%
% xi = 1:nPCs;
% xi = xi(sign(skewness(VicaNegE'))==-1);
% VicaNegE(xi,:) = -VicaNegE(xi,:);
%
% UicaNegE = Ur*Sr*Vr*pinv(VicaNegE);
% UicaNegENorm = nan(size(UicaNegE));
% for i = 1:nPCs
%     UicaNegENorm(:,i) = (UicaNegE(:,i) - min(UicaNegE(:,i)))/(max(UicaNegE(:,i))-min(UicaNegE(:,i)));
% end
% naturalsound_acoustic_category_figures(UicaNegENorm, ica_directory, ['RandSeed' num2str(i) '-MinMax-NegEntOpt-component-' num2str(i) '-']);

%%
%
%
% naturalsound_acoustic_category_figures(UneNorm, figure_directory, ['ICA-' num2str(nPCs) '-NegEntOpt-NormProfile-'],'corr_bounds',[-1 1]);
%
% % plot negentropy figures
% plot_rotation_negentropy(Vne,[0 0.3]);
% export_fig([figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');
%
% % plot negentropy figures
% plot_rotation_negentropy(Vr,[0 0.3]);
% export_fig([figure_directory 'PCA-' num2str(nPCs) '-NegEntPairs.pdf'],'-nocrop','-transparent','-pdf');
%
% % rotation figure for factorized distribution
% ResetRandStream(1);
% Vfac = nan(size(Vne));
% for i = 1:size(Vne,1)
%     Vfac(i,:) = Shuffle(Vne(i,:));
% end
% plot_rotation_negentropy(Vfac,[0 0.3]);
% export_fig([figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt-NegEntPairs-Factorized.pdf'],'-nocrop','-transparent','-pdf');
%
% % scatter
% myscatter2(Vne',si,'allpairs',3,[figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt-scatter']);
%
% % plot histogram by subject
% plothist_bysubject(Vne, si, [figure_directory 'ICA-' num2str(nPCs) '-NegEntOpt-hist-bysubject']);
%
%
%
% [U2,S2,V2] = svd(v_sg2,'econ');



% Select first N PCs, rescale Vs to have unit variance
% nPCs = 5;

% V2r = V2(:,1:nPCs)'*sqrt(size(V2,1)-1);
% U2r = U2(:,1:nPCs)/sqrt(size(V2,1)-1);
% S2r = S2(1:nPCs,1:nPCs);



% %%
%
% subplot(1,3,1);
% bar(var(Uica),'FaceColor',[1 1 1]);
% xlabel('ICA Components'); ylabel('Explained Variance');
% xlim([0 6]);
% kurtICA = abs(kurtosis(VicaNew')-3);
% subplot(1,3,2);
% bar(kurtICA,'FaceColor',[1 1 1]);
% xlabel('ICA Components'); ylabel('Kurtosis');
% xlim([0 6]);
% subplot(1,3,3)
% h = log(sqrt(2*pi*exp(1))) - entropy_SNH(VicaNew',10,50);
% bar(h,'FaceColor',[1 1 1]);
% xlabel('ICA Components'); ylabel('Negentropy (Nats)');
% xlim([0 6]);
% export_fig([ica_directory 'RandSeed' num2str(i) '-Rotated-component-' num2str(i) '_exvar_kurtosis_entropy.eps'],'-eps','-nocrop','-transparent');
%
% UicaNewNorm = nan(size(UicaNew));
% for i = 1:nPCs
%     UicaNewNorm(:,i) = (UicaNew(:,i) - min(UicaNew(:,i)))/(max(UicaNew(:,i))-min(UicaNew(:,i)));
% end


%%
% scatter plot

%
% % 2D Kurtosis
% n = 40;
% th = linspace(-pi/2,pi/2,n);
% phi = linspace(-pi/2,pi/2,n);
% kurt = nan(n,n);
% for i = 1:n
%     for j = 1:n
%         [x,y,z] = sph2cart(th(i),phi(j),1);
%         kurt(i,j) = abs(kurtosis( ([x,y,z]*Vr)' )-3);
%     end
% end
% kurtnorm = norm01(kurt);
% figure;
% hold on;
% c = colormap('jet(100)');
% for i = 1:size(kurt,1)
%     for j = 1:size(kurt,2)
%         [x,y,z] = sph2cart(th(i),phi(j),1);
%         plot(y,z,'o','Color', c(round(kurtnorm(i,j)*99)+1,:), 'LineWidth',5);
%     end
% end
% h = colorbar; set(h,'YTick',10:10:100,'YTickLabel',interp1(1:100, linspace(min(kurt(:)),max(kurt(:)),100), 10:10:100));
% for i = 1:3
%     plot(rotMatICA(i,2)*sign(rotMatICA(i,1)),rotMatICA(i,3)*sign(rotMatICA(i,1)),'ko','LineWidth',5);
% end
% export_fig([figure_directory 'ICA-' num2str(nPCs) '_kurtosis2D.pdf'],'-pdf','-nocrop');
%
% % 3D Kurtosis plot
% n = 40;
% th = linspace(-pi/2,pi/2,n)*2;
% phi = linspace(-pi/2,pi/2,n)*2;
% kurt = nan(n,n);
% for i = 1:n
%     for j = 1:n
%         [x,y,z] = sph2cart(th(i),phi(j),1);
%         kurt(i,j) = abs(kurtosis( ([x,y,z]*Vr)' )-3);
%     end
% end
% kurtnorm = norm01(kurt);
% figure;
% hold on;
% c = colormap('jet(100)');
% for i = 1:size(kurt,1)
%     for j = 1:size(kurt,2)
%         [x,y,z] = sph2cart(th(i),phi(j),1);
%         plot3(x,y,z,'o','Color', c(round(kurtnorm(i,j)*99)+1,:), 'LineWidth',10);
%     end
% end
% h = colorbar; set(h,'YTick',10:10:100,'YTickLabel',interp1(1:100, linspace(min(kurt(:)),max(kurt(:)),100), 10:10:100));


%%

% % block timing, runs
% blocktp = 3:6;
% runs = 1:2;
%
% matfile = [figure_directory 'voxel_matrix_unfilt.mat'];
% if ~exist(matfile,'file')
%     hemis = {'rh','lh'};
%     patch = cell(1,2); vras = cell(1,2); vras_patch = cell(1,2); vi = cell(1,2); vi_patch = cell(1,2); mask = cell(1,2);
%     grid_x = cell(1,2); grid_y = cell(1,2); bounds = cell(1,2); vgrid = cell(1,2);
%     for i = 1:2
%
%
%
%         % mask
%         mask{i} = read_label([params('rootdir') 'freesurfer/fsaverage/label/' hemis{i} '.' roi '.label']);
%         %         vras_patch{i} = patch{i}.vras(:,1:2);
%         %         vi_patch{i} = abs(patch{i}.vnums);% voxel responses
%
%         % vertices
%         % sign of patch vertices indicates something about position
%         % relative to the faces, not relevant for selecting the vertices
%         [~,ai,ib] = intersect(abs(patch{i}.vnums),mask{i}.vnums);
%         vras{i} = patch{i}.vras(ai,1:2); % 2-D coordinates
%         vi{i} = mask{i}.vnums(ib); % indices within roi
%
%         % bounds for regridding
%         bounds{i} = [floor(min(vras{i})); ceil(max(vras{i}))]; % min and max of 2-d coordinates
%         [grid_x{i},grid_y{i}] = meshgrid(bounds{i}(1,1):2:bounds{i}(2,1),(bounds{i}(1,2):2:bounds{i}(2,2))); % regrid in units of mm
%
%         % voxel grid
%         vgrid{i} = nan(size(grid_x{i},1), size(grid_y{i},2), length(usubs), 165, length(runs));
%
%         for k = 1:length(runs)
%             for j = 1:length(usubs)
%
%                 % subject id
%                 fprintf('hemi %s, usub %d, scan %d\n', hemis{i}, usubs(j), runs(k)); drawnow;

%% Sparse PCA
%
% addpath(genpath([pwd '/spasm']));
% vnorm = v ./ (ones(size(v,1),1)*sqrt(sum(v.^2)));
% nSPCs = 5;
% [Vspca SD L D paths] = spca(vnorm,[],nSPCs,0,-(round(size(v,2)/(nSPCs/3))));
% Uspca = vnorm*pinv(Vspca');
% naturalsound_acoustic_category_figures(Uspca(:,1:5), figure_directory, ['SPCA-' num2str(nSPCs) '-']);



% %%
%
% nPCs = 5;
% resolution = 11;
% nd = nPCs-1;
% angles = cell(1, nd);
% angle_grid = cell(1,nd);
% for i = 1:nd
%     angles{i} = linspace(-pi/2,pi/2,resolution);
% end
%
% [angle_grid{:}] = ndgrid(angles{:});
% angle_matrix = nan(resolution^nd, nd);
% for i = 1:nd
%     angle_matrix(:,i) = angle_grid{i}(:);
% end
%
% for i = 1:resolution^nd
% %     A*eye() = []
% end
%
% %%
%
% angles = angle_matrix(1,:);
%
% R = nan(nd+1,nd+1);
% % rotation matrix
% for i = 1:nd+1
%     ang = angles;
%     ang(1:i-1) = ang(1:i-1)+pi/2;
%     R(1,i) = cos(ang(1));
%     R(nd+1,i) = prod(sin(ang(1:nd)));
%     for j = 2:nd
%         R(j,i) = cos(ang(j))*prod(sin(angles(1:j-1)));
%     end
%

%% K-SVD
% addpath(genpath([pwd '/ompbox10']));
% addpath(genpath([pwd '/ksvdbox13']));
%
% % for k = 2:5
% params.data = v;
% params.Tdata = k;
% params.dictsize = 10;
% [D,Gamma,err,gerr] = ksvd(params);
% naturalsound_acoustic_category_figures(D, figure_directory, ['KSVD-' num2str(5) '-']);
% end

%                 subjid = ['naturalsound_us' num2str(usubs(j))];
%
%                 % voxel responses
%                 x = MRIread([params('rootdir') 'freesurfer/fsaverage/fla/' subjid '/main_v3_combined_r' num2str(runs(k)) '/' hemis{i} '.sigav_demean_' num2str(100*fwhm, '%.0f') 'mm_tp' num2str(blocktp(1)) '-' num2str(blocktp(end)) '.mgz']);
%                 vr = squeeze(x.vol(:,vi{i}+1,:,:))'; % vertex indices are zero-indexed
%                 %                 vr = squeeze(x.vol(:,vi{i}+1,:,:))';
%                 %                 vr_patch = squeeze(x.vol(:,vi_patch{i}+1,:,:))';
%
%                 % add boundary points as NaN values, useful for
%                 % interpolation
%                 vras_boundary_pts = [vras{i}; bounds{i}(1,:); bounds{i}(2,:); bounds{i}(1,1), bounds{i}(2,2); bounds{i}(2,1), bounds{i}(1,2)];
%                 vr_boundary_pts = [vr,nan(165,4)];
%
%                 % grided nan mask
%                 mask = all(~isnan(vr_boundary_pts));
%                 grid_mask = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),double(mask),grid_x{i},grid_y{i},'natural');
%
%                 % regridding
%                 for q = 1:165
%                     fprintf('%d ', q); drawnow;
%                     grid_response = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),vr_boundary_pts(q,:),grid_x{i},grid_y{i},'natural');
%                     grid_response(grid_mask<0.99) = NaN;
%                     vgrid{i}(:,:,j,q,k) = grid_response;
%                     %                     figure;
%                     %                     surf(grid_x{i}, grid_y{i}, 100*grid_response,'EdgeColor','none');
%                     %                     caxis([0 3]);
%                     %                     view(0, 90);
%                 end
%
%                 %                 vgrid{i}(:,:,j,q,k) = griddata(vras{i}(:,1),vras{i}(:,2),vr(q,:),grid_x{i},grid_y{i});
%                 %                 vgrid{i}(:,:,j,q,k) = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),vr_boundary_pts(q,:),grid_x{i},grid_y{i});
%
%                 %                 g = griddata(vras_patch{i}(:,1),vras_patch{i}(:,2),vr_patch(1,:),grid_x{i},grid_y{i});
%                 %                 g = griddata(vras{i}(:,1),vras{i}(:,2),vr(1,:),grid_x{i},grid_y{i});
%
%                 %                 figure;
%                 %                 q = 1;
%                 %                 scatter(vras_boundary_pts(:,1),vras_boundary_pts(:,2),3,100*vr_boundary_pts(q,:));
%                 %                 caxis([0 3]);
%                 %
%                 %                 figure;
%                 %                 q = 1;
%                 %                 scatter(vras_boundary_pts(:,1),vras_boundary_pts(:,2),3,double(mask));
%                 %                 caxis([0 1]);
%                 %

%%

    %
    %     % sound response
    %     [~,p] = ttest(sv,0,'tail','right');
    %     sound_responseP = -log10(p);
    %
    %     % select voxels with a significant sound response
    %     sv_selected = sv(:,sound_responseP > soundthresh);
    %     %     sv_demeaned = demean_subjects(sv_selected, ones(1,size(sv_selected,2)));
    %     %     sv_selected_norm = sv_selected ./ (ones(165,1)*sum(sv_selected));
    %
    %     %             Uica_norm = Une ./ (ones(165,1)*std(Une));
    %     %         %         [b, ~, logP, ~] = regress_SNH(Uica_norm, sv_demeaned);
    %     %         [b, ~, logP, ~] = regress_SNH([Uica_norm, grandmean/std(grandmean)], sv_selected);
    %
    %     % Regression analysis
    %     UneUnitVar = Une ./ (ones(165,1)*std(Une));
    %     xhat = pinv([Uica_norm])*(sv_selected_norm-mean(sv_selected_norm,2)*ones(1,size(sv_selected_norm,2)));
    %     xhat = pinv([Uica_norm])*sv_norm_demeaned;
%%
   %         [~,p] = ttest(ica_weights(:,:,:,i),0,'dim',3);
    %         p(isnan(p)) = 0;
    %         ica_logP_group(:,:,i) = -log10(p);
    %
    %         fname = [surf_directory idstring 'logP-component-' num2str(k) '-' hemis{i} '.mgz'];
    %         MRIwrite_SNH(ica_logP_group(:,:,i), fname, hemis{i});
    %
    %         % save screenshot
    %         surface_image = [surface_projections_directory idstring 'logP-component-' num2str(k) '-group-' hemis{i} '.png'];
    %         if false || ~exist(surface_image, 'file')
    %             %                 freeview3('fsaverage',hemis{i},'overlay',surface_file,'screenshot',surface_image,'overlay_threshold',[0.001 1.5 3],'piecewise');
    %             freeview3('fsaverage',hemis{i},'overlay',fname,'overlay_threshold',[3 4.5 6]);
    %         end
    
    %         surface_file = ;
    %         ica_weight_surf = x1;
    %         dims = size(x1.vol);
    %         ica_weight_surf.vol = nan(dims(1:3));
    %         ica_weight_surf.vol(1,:,1) = ica_weights_group(:,:,i)';
    %         ica_weight_surf.nframes = size(ica_weight_surf.vol,4);
    %         ica_weight_surf.fspec = surface_file;
    %         MRIwrite(ica_weight_surf, surface_file);
    
    % save screenshot
    %         surface_image = [surface_projections_directory idstring 'logP-component-' num2str(k) '-us' num2str(usubs) '.png'];
    %         if false || ~exist(surface_image, 'file')
    %             %                 freeview3('fsaverage',hemis{i},'overlay',surface_file,'screenshot',surface_image,'overlay_threshold',[0.001 1.5 3],'piecewise');
    %             freeview3('fsaverage',hemis{i},'overlay',surface_file,'overlay_threshold',[3 4.5 6]);
    %         end
    %
    %         surface_file = [surf_directory idstring 'logP-component-' num2str(k) '.mgz'];
    %         ica_logP_surf = x1;
    %         dims = size(x1.vol);
    %         ica_logP_surf.vol = nan(dims(1:3));
    %         ica_logP_surf.vol(1,:,1) = ica_logP(k,:,j,i)';
    %         ica_logP_surf.nframes = size(ica_logP_surf.vol,4);
    %         ica_logP_surf.fspec = surface_file;
    %         MRIwrite(ica_logP_surf, surface_file);
    %
    %         % save screenshot
    %         surface_image = [surface_projections_directory idstring 'logP-component-' num2str(k) '-us' num2str(usubs) '.png'];
    %         if false || ~exist(surface_image, 'file')
    %             %                 freeview3('fsaverage',hemis{i},'overlay',surface_file,'screenshot',surface_image,'overlay_threshold',[0.001 1.5 3],'piecewise');
    %             freeview3('fsaverage',hemis{i},'overlay',surface_file,'overlay_threshold',[3 4.5 6]);
    %         end
    %

% %% Select stimuli with maximum variance
% 
% n_stimuli = 36;
% UneZscore = zscore(Une);
% available_stimuli = 1:165;
% selected_stimuli = nan(n_stimuli,1);
% max_variance = nan(n_stimuli/2,1);
% 
% for i = 1:n_stimuli/2
%     
%     n = length(available_stimuli);
%     xi = available_stimuli'*ones(1,n);
%     yi = ones(n,1)*available_stimuli;
%     
%     % npairs x nICs x nstimuli matrix
%     x = nan(n^2, size(UneZscore,2), i*2);
%     for j = 1:(i-1)*2
%         x(:,:,j) = UneZscore(selected_stimuli(j)*ones(n,n), :);
%     end
%     x(:,:,(i-1)*2+1) = UneZscore(xi, :);
%     x(:,:,(i-1)*2+2) = UneZscore(yi, :);
%     
%     mean_variance_allpairs = mean(var(x,[],3),2);
%     [max_variance(i),bestpair] = max(mean_variance_allpairs);
%     selected_stimuli((1:2) + (i-1)*2) = [xi(bestpair), yi(bestpair)];
%     available_stimuli = setdiff(1:165, selected_stimuli(1:i*2));
%     
% end
% 
% figure;
% set(gcf,'Position', [0 0 1400 800]);
% subplot(1,3,1);
% [y,bins] = hist(UneZscore, 10);
% plot(bins,y/165,'LineWidth',2);
% xlabel('Z values');
% title('Histogram of Stimulus Values');
% legend(component_names,'Location','Best');
% set(gca, 'FontSize', 14);
% box off;
% 
% subplot(1,3,2);
% [y,bins] = hist(UneZscore(selected_stimuli,:), 10);
% plot(bins,y/sum(y(:,1)),'LineWidth',2);
% xlabel('Z values');
% title('Histogram of Stimulus Values for Selected Subset')
% set(gca, 'FontSize', 14);
% box off;
% 
% subplot(1,3,3);
% x = UneZscore(selected_stimuli,:);
% [y,bins] = hist(x(:), 10);
% plot(bins,y/sum(y(:,1)),'k','LineWidth',2);
% xlabel('Z values');
% title('All Stimulus Values for Selected Set')
% set(gca, 'FontSize', 14);
% box off;
% export_fig([figure_directory 'histogram_stimulus_weights_' num2str(n_stimuli) 'stims.pdf'],'-pdf','-nocrop','-transparent');
% 
% % category figures
% figure;
% set(gcf,'Position', [0 0 1400 800]);
% subplot(1,2,1);
% [R2, C, conds] = naturalsound_acoustic_category_regressors;
% [y,bins] = hist(C.category_assignments,11);
% hold on;
% for i = 1:11
%     cat = C.plotting_order(i);
%     bar(i,y(cat)/165,'FaceColor',C.colors(cat,:));
% end
% xlim([0 12]);
% title('All Stimuli');
% legend(C.category_labels(C.plotting_order),'Location','Best');
% set(gca, 'FontSize', 14);
% box off;
% 
% subplot(1,2,2);
% [y,bins] = hist(C.category_assignments(selected_stimuli),11);
% hold on;
% for i = 1:11
%     cat = C.plotting_order(i);
%     bar(i,y(cat)/n_stimuli,'FaceColor',C.colors(cat,:));
% end
% xlim([0 12]);
% title('Selected Stimuli');
% set(gca, 'FontSize', 14);
% box off;
% export_fig([figure_directory 'category_histogram_' num2str(n_stimuli) 'stims.pdf'],'-pdf','-nocrop','-transparent');
% 
% % acoustic and category figures for subset of stimuli
% Csubset = C;
% Csubset.category_assignments = C.category_assignments(selected_stimuli);
% Csubset.ids = C.ids(selected_stimuli);
% Csubset.category_regressors = C.category_regressors(selected_stimuli,:);
% 
% R2subset = R2;
% R2subset.ids = R2.ids(selected_stimuli);
% R2subset.stats = R2.stats(selected_stimuli,:);
% 
% conds_subset = conds(selected_stimuli);
% 
% for i = 1:size(UneNorm,2)
%     y = UneNorm(selected_stimuli,i);
%     response_profile_figure(y, Csubset, conds_subset);
%     export_fig([figure_directory 'ICA-6-stimulus-subset-' num2str(n_stimuli) 'stims-' component_names{i} '_response_profile.eps'],'-eps','-nocrop','-transparent');
%     acoustic_category_figures(y, Csubset, R2subset, [figure_directory 'ICA-6-stimulus-subset-' num2str(n_stimuli) 'stims-' component_names{i}''], varargin{:});
% end
% 
% stim_names = conds(selected_stimuli);
% ids = C.ids(selected_stimuli);
% save([figure_directory 'selected_stimuli_' num2str(n_stimuli) 'stims.mat'],'stim_names','ids');
%                 %                 figure;
%                 %                 gx_mask = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),double(mask),grid_x{i},grid_y{i},'natural');
%                 %                 surf(grid_x{1}, grid_y{1}, gx_mask,'EdgeColor','none');
%                 %                 caxis([0 1]);
%                 %                 view(0, 90);
%                 %
%                 %                 figure;
%                 %                 gx = griddata(vras_boundary_pts(:,1),vras_boundary_pts(:,2),vr_boundary_pts(q,:),grid_x{i},grid_y{i},'natural');
%                 %                 surf(grid_x{1}, grid_y{1}, 100*gx.*double(gx_mask>0.99),'EdgeColor','none');
%                 %                 view(0, 90);
%                 %                 caxis([0 3]);
%
%                 %                 figure;
%                 %                 q = 1;
%                 %                 scatter(vras{i}(:,1),vras{i}(:,2),3,100*vr(q,:));
%                 %                 caxis([0 3]);
%                 %
%                 %                 figure;
%                 %                 q = 1;
%                 %                 scatter(vras{i}(:,1),vras{i}(:,2),3,double(~isnan(vr(q,:))));
%                 %                 caxis([0 3]);
%                 %
%                 %                 figure;
%                 %                 mask = all(~isnan(vr));
%                 %                 gx_mask = griddata(vras{i}(:,1),vras{i}(:,2),double(mask),grid_x{i},grid_y{i},'natural');
%                 %                 surf(grid_x{1}, grid_y{1}, gx_mask,'EdgeColor','none');
%                 %                 caxis([0 1])
%                 %                 view(0, 90);
%                 %                 figure;
%                 %                 gx = griddata(vras{i}(:,1),vras{i}(:,2),vr(q,:),grid_x{i},grid_y{i},'natural');
%                 %                 surf(grid_x{1}, grid_y{1}, 100*gx.*double(gx_mask>0.99),'EdgeColor','none');
%                 %                 view(0, 90);
%                 %                 caxis([0 3]);
%
%                 %                 figure;
%                 %                 scatter(vras{i}(:,1),vras{i}(:,2),3,vr(1,:));
%                 %                 figure;
%                 %                 scatter(vras_patch{i}(:,1),vras_patch{i}(:,2),3,vr_patch(1,:));
%                 %                 figure;
%                 %                 surf(grid_x{i},grid_y{i},g);
%                 %                 view(0, 90);
%                 %                 figure;
%                 %                 mesh(grid_x{i},grid_y{i},m);
%                 %                 view(0, 90);
%
%             end
%         end
%     end
%     save(matfile, 'vras','vi','mask','grid_x','grid_y','vgrid');
% else
%     load(matfile);
% end
