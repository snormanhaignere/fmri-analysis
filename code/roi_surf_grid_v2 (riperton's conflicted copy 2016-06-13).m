function [psc, conditions, n_voxels_per_run_and_threshold] = ...
    roi_surf_grid_v2(  us, grid_roi, grid_spacing_mm, ...
    localizer_info, test_info, fwhm, varargin )

fprintf('Running roi_surg_grid_v2\n\n');
drawnow;

% New version that is more coded more naturally such that the files that perform
% the relevant first and second level analyses are call directly by this
% function and return the appropriate files needed for the ROI analysis. This
% function also supports permtutation-based stats (with the flag 'perm_stats'),
% and returns the number of voxels in each ROI computed.
%
% fprintf('v2\n'); drawnow;

% % -- Example Arguments --
% % subject
% us = 158;
%
% % roi constraint region and spacing of the grid
% grid_roi = 'hand-audctx';
% grid_spacing_mm = 2.8571/2;
%
% % amount the data was smoothed
% fwhm = 2.8571;
%
% % information about the test data
% clear test_info;
% test_info.exp = 'pitch_localizer_monkey';
% test_info.runtype = 'pitchloc2_combined_split';
%
% % information about the localizer data/contrasts
% clear localizer_info;
% localizer_info = [test_info, test_info];
%
% localizer_info(1).contrast = 'all_pitchloc2';
% localizer_info(1).contrast_sign = 1;
% localizer_info(1).threshold_type = 'absolute';
% localizer_info(1).thresholds = 3;
% localizer_info(1).max_runs_to_use = [];
%
% localizer_info(2).contrast = 'f0100-200_vs_f0800-1600';
% localizer_info(2).contrast_sign = 1;
% localizer_info(2).threshold_type = 'relative';
% localizer_info(2).thresholds = 0.05:0.05:0.3;
% localizer_info(2).max_runs_to_use = [];

% default parameters
test_info = ...
    default_test_parameters(test_info, us, varargin{:});
localizer_info = ...
    default_localizer_parameters(localizer_info, us, test_info, varargin{:});

% number of thresholds for each localizer
n_localizers = length(localizer_info);
n_thresholds_per_localizer = nan(1, n_localizers);
for j = 1:n_localizers
    n_thresholds_per_localizer(j) = length(localizer_info(j).thresholds);
end

% initialize PSC matrix
% runs x conditions x thresholds
psc = nan([ length(test_info.runs), length(test_info.conditions), ...
    n_thresholds_per_localizer ]);

% number of voxels selected for different relative thresholds
n_voxels_per_run_and_threshold = ...
    nan([length(test_info.runs), n_thresholds_per_localizer]);

for k = 1:length(test_info.runs) % loop through runs
    
    if ~optInputs(varargin, 'suppress-output')
        % print information about this localizer
        fprintf('Test: %s, run %d\n', test_info.runtype, test_info.runs(k));
        drawnow;
    end
    
    psc_file = test_psc_file(test_info, us, test_info.runs(k), ...
        fwhm, grid_roi, grid_spacing_mm, varargin{:});
            
    % matrix of psc values for each voxels
    % condition x voxel matrix
    % voxels unwrapped from grid and concatenated across hemispheres
    % right hemisphere first, then left
    pscmap = load(psc_file, 'G');
    n_voxels = numel(pscmap.G.grid_x{1}) + numel(pscmap.G.grid_x{2});
    voxel_psc = nan(length(test_info.conditions), n_voxels);
    conditions_for_this_run = ...
        read_conditions(test_info.exp, us, test_info.runtype, ...
        'run', test_info.runs(k), varargin{:});
    
    for j = 1:length(test_info.conditions)
        xi = strcmp(test_info.conditions{j}, conditions_for_this_run);
        responses_rh = pscmap.G.grid_data{1}(:,:,xi);
        responses_lh = pscmap.G.grid_data{2}(:,:,xi);
        voxel_psc(xi, :) = [responses_rh(:)', responses_lh(:)'];
    end
    clear pscmap;
    
    % set values that are exactly zero to NaN
    % these values should not be included in the roi
    voxel_psc(voxel_psc==0) = NaN;
        
    % read in the localizer contrast matrix
    localizer_contrast_stat_matrix = nan(n_localizers, n_voxels);
    for j = 1:n_localizers
        
        % ensure non-independence
        if strcmp(test_info.exp, localizer_info(j).exp) && strcmp(test_info.runtype, localizer_info(j).runtype)
            localizer_runs_to_use = setdiff(localizer_info(j).runs, test_info.runs(k));
        else
            localizer_runs_to_use = localizer_info(j).runs;
        end
        
        % subsample the localizer runs
        % (e.g. in order to match response reliability)
        if ~isinf(localizer_info(j).max_runs_to_use)
            
            % distance of localizer runs to test run
            dist = (localizer_runs_to_use - test_info.runs(k)).^2;
            
            try
                % which runs to use based on distance
                switch localizer_info(j).use_nearest_or_farthest
                    case 'nearest'
                        [~,xi] = sort(dist, 'ascend');
                    case 'farthest'
                        fprintf('Using the farthest runs\n'); drawnow;
                        [~,xi] = sort(dist, 'descend');
                    otherwise
                        error('nearest_or_farthest_runs cannot be %s', ...
                            localizer_info(j).use_nearest_or_farthest);
                end
                localizer_runs_to_use = ...
                    localizer_runs_to_use(xi(1:localizer_info(j).max_runs_to_use));
                clear dist;
            catch
                keyboard
            end
            
        end
        
        % check there is at least one usable localizer run
        if isempty(localizer_runs_to_use)
            error('There needs to be at least one usable localizer run.\n');
        end
        
        % print the runs being used
        if ~optInputs(varargin, 'suppress-output')
            fprintf('Localizer: %s, runs %s\n', ...
                localizer_info(j).contrast, sprintf('%d',localizer_runs_to_use));
            drawnow;
        end
        
        fprintf('Finding pstat\n\n'); drawnow;
        
        % file with p-values for a given localizer
        pstat_file = localizer_pstat_file(...
            localizer_info(j), us, localizer_runs_to_use, ...
            fwhm, grid_spacing_mm, grid_roi, varargin{:});
        
        % fprintf('roi_surf_grid_v2.m\n'); drawnow;
        % keyboard;
        
        % load and store significance map
        statmap = load(pstat_file, 'G');
        localizer_contrast_stat_matrix(j,:) = ...
            [statmap.G.grid_data{1}(:)', statmap.G.grid_data{2}(:)'];
        clear statmap;
                
    end
    
    % set values that are exactly zero to NaN
    localizer_contrast_stat_matrix(localizer_contrast_stat_matrix==0) = NaN;
        
    % loop through all combinations of thresholds, selecting voxels, and
    % measuring mean PSC values
    for i = 1:prod(n_thresholds_per_localizer)
        
        % the index into the threshold vector for each localizer
        threshold_indices = cell(1,n_localizers);
        [threshold_indices{:}] = ind2sub(n_thresholds_per_localizer, i);
        
        % remove NaN voxels
        % logical all is applied to conditions for each voxel
        voxels_in_roi = find(...
            all(~isnan(voxel_psc),1) ...
            & all(~isnan(localizer_contrast_stat_matrix),1));
        
        % loop through the localizers
        for j = 1:n_localizers
            
            % relavent info for single localizer
            thresh = localizer_info(j).thresholds(threshold_indices{j});
            stat = localizer_contrast_stat_matrix(j,:) * sign(localizer_info(j).contrast_sign);
            
            % check all statistics are not NaN
            assert(all(~isnan(stat(voxels_in_roi(:)))));
            
            try
                % select voxels from those left
                switch localizer_info(j).threshold_type
                    case 'absolute'
                        xi = stat(voxels_in_roi) > thresh;
                        voxels_in_roi = voxels_in_roi(xi);
                    case 'relative'
                        n_voxels_to_select = round( thresh * length(voxels_in_roi) );
                        [~,xi] = sort(stat(voxels_in_roi), 'descend');
                        voxels_in_roi = voxels_in_roi(xi(1:n_voxels_to_select));
                    otherwise
                        error('Localizer "selection type" should be "absolute-threshold" or "relative-threshold" not %s', localizer_info(j).threshold_type);
                end
            catch
                keyboard
            end
            
        end
        
        n_voxels_after_selection = length(voxels_in_roi);
        n_voxels_per_run_and_threshold(k, threshold_indices{:}) = n_voxels_after_selection;
        
        if isempty(voxels_in_roi)
            warning('No Voxels in ROI');
            continue;
        end
        
        voxel_psc_in_roi = voxel_psc(:, voxels_in_roi);
        
        % check voxels do not have NaN values
        assert(all(~isnan(voxel_psc_in_roi(:))));
        
        % average PSC values within selected voxels
        psc(k,:,threshold_indices{:}) = mean(voxel_psc_in_roi, 2);
        
    end
end

% remove single dimensions for the thresholds
psc = reshape(psc, [length(test_info.runs), length(test_info.conditions), ...
    setdiff(n_thresholds_per_localizer,1)]);
n_voxels_per_run_and_threshold = reshape(n_voxels_per_run_and_threshold, ...
    [length(test_info.runs), setdiff(n_thresholds_per_localizer,1)]);

% return conditions used
conditions = test_info.conditions;

% find/compute desired psc file
function psc_file = test_psc_file(test_info, us, test_run, ...
    fwhm, grid_roi, grid_spacing_mm, varargin)

% input and output file to the first level analysis
fla_input_fname = [...
    'motcorr_smooth' num2str(100*fwhm, '%.0f') ...
    'mm_grid_' grid_roi '_' num2str(grid_spacing_mm) 'mm'];
fla_output_fname = [...
    'smooth' num2str(100*fwhm, '%.0f') ...
    'mm_grid_' grid_roi '_' num2str(grid_spacing_mm) 'mm_10whitematterPCs'];

% hrf type
if optInputs(varargin, 'monkey')
    hrf = 'MION_CUSTOM1';
else
    hrf = 'BOLD';
end

% first level analysis script
% returns a structure with the needed file
S = fla_matlab(...
    test_info.exp, us, test_info.runtype, ...
    test_run, hrf, 'downsampled_surface',...
    fla_input_fname, fla_output_fname, ...
    'n_whitematter_PCs', 10, 'motcorr_highres_2mm',...
    'contrasts', {}, ...
    'permute_condition_order', 'n_cpus', 1, 'n_perms', 100,...
    varargin{:});

% return the psc file
psc_file = S.psc_file;

% helper function that find the appropriate file with p-values
function pstat_file = localizer_pstat_file(...
    localizer_info, us, localizer_runs_to_use, ...
    fwhm, grid_spacing_mm, grid_roi, varargin)

% stat file from first or second level analysis
if length(localizer_runs_to_use) == 1 % first level analysis
    
    pstat_file = fla_pstat_file(...
        localizer_info, us, localizer_runs_to_use, ...
        fwhm, grid_spacing_mm, grid_roi, varargin{:});
    
else % second level analysis
    
    pstat_file = sla_pstat_file(...
        localizer_info, us, localizer_runs_to_use, ...
        fwhm, grid_spacing_mm, grid_roi, varargin{:});
    
end

% localizer file for a localizer with a single run
function pstat_file = fla_pstat_file(...
    localizer_info, us, localizer_runs_to_use, ...
    fwhm, grid_spacing_mm, grid_roi, varargin)

% input and output file to the first level analysis
fla_input_fname = [...
    'motcorr_smooth' num2str(100*fwhm, '%.0f') ...
    'mm_grid_' grid_roi '_' num2str(grid_spacing_mm) 'mm'];
fla_output_fname = [...
    'smooth' num2str(100*fwhm, '%.0f') ...
    'mm_grid_' grid_roi '_' num2str(grid_spacing_mm) 'mm_10whitematterPCs'];

% hrf type
if optInputs(varargin, 'monkey')
    hrf = 'MION_CUSTOM1';
else
    hrf = 'BOLD';
end

% first level analysis script
% returns a structure with the needed file
S = fla_matlab(...
    localizer_info.exp, us, localizer_info.runtype, ...
    localizer_runs_to_use, hrf, 'downsampled_surface',...
    fla_input_fname, fla_output_fname, ...
    'n_whitematter_PCs', 10, 'motcorr_highres_2mm',...
    'contrasts',{localizer_info.contrast},...
    'permute_condition_order', 'n_cpus', 1, 'n_perms', 10e3,...
    varargin{:});

clear fla_input_fname fla_output_fname;

% return parametric or permutation-based contrast map
switch localizer_info.stat_type
    case 'param'
        
        xi = ismember(S.p_file(:,1), localizer_info.contrast);
        pstat_file = S.p_file{xi,2};
        clear xi;
        
    case 'perm'
        
        xi = ismember(S.p_file_permuted_gaussfit(:,1), localizer_info.contrast);
        pstat_file = S.p_file_permuted_gaussfit{xi,2};
        clear xi;
        
    otherwise
        error('stat_type cannot be %s', localizer_info.stat_type)
        
end

% localizer file for a localizer with many runs
function pstat_file = sla_pstat_file(...
    localizer_info, us, localizer_runs_to_use, ...
    fwhm, grid_spacing_mm, grid_roi, varargin)

grid_file = ['smooth' num2str(fwhm*100, '%.0f') 'mm_grid_' grid_roi ...
    '_' num2str(grid_spacing_mm) 'mm_10whitematterPCs'];

% return parametric or permutation-based contrast map
switch localizer_info.stat_type
    case 'param'
        
        pstat_file = sla_matlab(...
            localizer_info.exp, us, ...
            localizer_info.runtype, localizer_info.contrast, ...
            localizer_runs_to_use,...
            'downsampled_surface',grid_file,grid_file,...
            'fixed_without_whitening',varargin{:});
        
    case 'perm'
        
        % create parameter structure
        clear P;
        for q = 1:length(localizer_runs_to_use);
            P(q).exp = localizer_info.exp; %#ok<*AGROW,*SAGROW>
            P(q).us = us;
            P(q).runtype = localizer_info.runtype;
            P(q).runs = localizer_runs_to_use(q);
            P(q).contrast = localizer_info.contrast;
            P(q).lower_level_directory_name = grid_file;
        end
        
        % string to identify the localizer_info runs to use
        % when there are more than 50 runs
        % the string is written in a more compressed form
        if length(localizer_runs_to_use) > 50
            localizer_run_string = [num2str(length(localizer_runs_to_use)) ...
                'r-' num2str(localizer_runs_to_use(1)) ...
                '-' num2str(localizer_runs_to_use(end)) ...
                '_' DataHash(localizer_runs_to_use)];
        else
            localizer_run_string = ['r' sprintf('%d',localizer_runs_to_use)];
        end
        
        % number of permutations per run
        if optInputs(varargin, 'monkey');
            n_perms = 100;
        else
            n_perms = 10e3;
        end
        
        n_smps = 10e3;
        
        tic;
        
        % output directory name
        tla_directory_name = [...
            localizer_info.contrast '_' ...
            localizer_info.exp '_' num2str(us) '_' ...
            localizer_info.runtype '_' localizer_run_string ...
            '_' num2str(n_perms) 'perms_' num2str(n_smps) 'smps'];
        
        fprintf('Running Permutation test\n\n');
        drawnow;
        
        % call tla_permutation
        [~, ~, ~, pstat_file] = ...
            tla_permutation_voxelwise_stats(...
            P,'downsampled_surface',tla_directory_name,...
            'no_stat_plots', 'n_smps', n_smps, 'n_perms', n_perms,...
            varargin{:});
        
        toc;
        
    otherwise
        error('stat_type cannot be %s', localizer_info.stat_type)
        
end

function localizer_info = ...
    default_localizer_parameters(localizer_info, us, test_info, varargin)

n_localizers = length(localizer_info);

for j = 1:n_localizers
    if ~isfield(localizer_info(j), 'exp') || isempty(localizer_info(j).exp)
        localizer_info(j).exp = test_info.exp;
    end
    
    if ~isfield(localizer_info(j), 'runtype') || ...
            isempty(localizer_info(j).runtype)
        localizer_info(j).runtype = test_info.runtype;
    end
    
    if ~isfield(localizer_info(j), 'runs') || isempty(localizer_info(j).runs)
        localizer_info(j).runs = read_runs(...
            localizer_info(j).exp, us, localizer_info(j).runtype, varargin{:});
    end
    
    if ~isfield(localizer_info(j), 'max_runs_to_use') || ...
            isempty(localizer_info(j).max_runs_to_use)
        localizer_info(j).max_runs_to_use = inf;
    end
    
    if ~isfield(localizer_info(j), 'contrast_sign') || ...
            isempty(localizer_info(j).contrast_sign)
        localizer_info(j).contrast_sign = 1;
    end
    
    if ~isfield(localizer_info(j), 'use_nearest_or_farthest') || ...
            isempty(localizer_info(j).use_nearest_or_farthest)
        localizer_info(j).use_nearest_or_farthest = 'nearest';
    end
end

function test_info = default_test_parameters(test_info, us, varargin)

% runs to use
if ~isfield(test_info, 'runs')
    test_info.runs = read_runs(...
        test_info.exp, us, test_info.runtype, varargin{:});
end

% conditions to measure responses to
if ~isfield(test_info, 'conditions')
    test_info.conditions = read_conditions(...
        test_info.exp, us, test_info.runtype, varargin{:});
end