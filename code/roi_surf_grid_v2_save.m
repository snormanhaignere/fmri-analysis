function [psc, conditions] = ...
    roi_surf_grid_v2(  us, grid_roi, grid_spacing_mm, ...
    localizer_info, test_info, fwhm, varargin )

% New version that is more coded a bit more naturally so that the link the
% contrast and psc files that are loaded are more directly tied to underlying
% analyses.
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

% -- Test data parameters --

% default parameters for test data, if not specified
if ~isfield(test_info, 'runs')
    test_info.runs = read_runs(...
        test_info.exp, us, test_info.runtype, varargin{:});
end
if ~isfield(test_info, 'conditions')
    test_info.conditions = read_conditions(...
        test_info.exp, us, test_info.runtype, varargin{:});
end

% directory with information from the first level analysis of the test data
if optInputs(varargin, 'monkey')
    test_fla_directory = [params('rootdir') 'freesurfer/'...
        test_info.exp '_us' num2str(us) '/fla_matlab'];
else
    test_fla_directory = [params('rootdir') 'freesurfer/' ...
        'fsaverage/fla_matlab/' test_info.exp '_us' num2str(us)];
end

% -- Localizer data parameters --

% default parameters for localizer data, if not specified
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
end

% number of thresholds for each localizer
n_thresholds_per_localizer = nan(1, n_localizers);
for j = 1:n_localizers
    n_thresholds_per_localizer(j) = length(localizer_info(j).thresholds);
end

% initialize PSC matrix
% runs x conditions x thresholds
psc = nan([ length(test_info.runs), length(test_info.conditions), n_thresholds_per_localizer ]);

for k = 1:length(test_info.runs) % loop through runs
    
    % mat file with psc values
    psc_mat_file = [test_fla_directory '/' test_info.runtype '_r' num2str(test_info.runs(k)) ...
        '/smooth' num2str(100*fwhm, '%.0f') 'mm_grid_' grid_roi '_' num2str(grid_spacing_mm) 'mm_10whitematterPCs_downsampled'...
        '/psc.mat'];
    
    % matrix of psc values for each voxels
    % condition x voxel matrix
    % voxels unwrapped from grid and concatenated across hemispheres
    % right hemisphere first, then left
    pscmap = load(psc_mat_file, 'G');
    n_voxels = numel(pscmap.G.grid_x{1}) + numel(pscmap.G.grid_x{2});
    voxel_psc = nan(length(test_info.conditions), n_voxels);
    conditions_for_this_run = read_conditions(test_info.exp, us, test_info.runtype, 'run', test_info.runs(k), varargin{:});
    for j = 1:length(test_info.conditions)
        xi = strcmp(test_info.conditions{j}, conditions_for_this_run);
        responses_rh = pscmap.G.grid_data{1}(:,:,xi);
        responses_lh = pscmap.G.grid_data{2}(:,:,xi);
        voxel_psc(xi, :) = [responses_rh(:)', responses_lh(:)'];
    end
    clear PSC;
    
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
            [~,xi] = sort((localizer_runs_to_use - test_info.runs(k)).^2, 'ascend');
            localizer_runs_to_use = localizer_runs_to_use(xi(1:localizer_info(j).max_runs_to_use));
        end
        
        % check there is at least one usable localizer run
        if isempty(localizer_runs_to_use)
            error('There needs to be at least one usable localizer run.\n');
        end
        
        % file with p-values for a given localizer
        pstat_file = find_pstat_file(...
            localizer_info(j), us, localizer_runs_to_use, ...
            fwhm, grid_spacing_mm, grid_roi, varargin{:});
        
        % load and store significance map
        statmap = load(pstat_file, 'G');
        localizer_contrast_stat_matrix(j,:) = [statmap.G.grid_data{1}(:)', statmap.G.grid_data{2}(:)'];
        clear statmap;
        
    end
    
    for i = 1:prod(n_thresholds_per_localizer) % loop through all combinations of thresholds
        
        % the index into the threshold vector for each localizer
        threshold_indices = cell(1,n_localizers);
        [threshold_indices{:}] = ind2sub(n_thresholds_per_localizer, i);
        
        % remove NaN voxels
        % logical all is applied to conditions for each voxel
        voxels_in_roi = find(all(~isnan(voxel_psc),1));
        n_starting_voxels = size(voxels_in_roi,2);
        
        % loop through the localizers
        for j = 1:n_localizers
            
            % relavent info for single localizer
            thresh = localizer_info(j).thresholds(threshold_indices{j});
            stat = localizer_contrast_stat_matrix(j,:) * sign(localizer_info(j).contrast_sign);
            
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
        end
        
        n_voxels_after_selection = length(voxels_in_roi);
        
        if ~optInputs(varargin, 'suppress-output')
            
            % print information about this localizer
            fprintf('Test: %s, run %d\n', test_info.runtype, test_info.runs(k));
            for j = 1:n_localizers
                fprintf('Localizer: %s, runs %s\n', localizer_info(j).contrast, sprintf('%d',localizer_runs_to_use));
                fprintf('Threshold: %.3f\n', localizer_info(j).thresholds(threshold_indices{j}));
            end
            fprintf('starting number of voxels: %d\n', n_starting_voxels);
            fprintf('number of voxels after selection: %d\n\n\n',n_voxels_after_selection);
            drawnow;
            
        end
        
        % average PSC values within selected voxels
        psc(k,:,threshold_indices{:}) = mean(voxel_psc(:, voxels_in_roi), 2);
    end
end

% remove single dimensions for the thresholds
psc = reshape(psc, [length(test_info.runs), length(test_info.conditions), setdiff(n_thresholds_per_localizer,1)]);

conditions = test_info.conditions;

% helper function that find the appropriate file with p-values
function pstat_file = find_pstat_file(...
    localizer_info, us, localizer_runs_to_use, ...
    fwhm, grid_spacing_mm, grid_roi, varargin)

% stat file from first or second level analysis
if length(localizer_runs_to_use) == 1 % first level analysis
    
    % input and output file to the first level analysis
    fla_input_fname = [...
        'motcorr_smooth' num2str(100*fwhm, '%.0f') ...
        'mm_grid_hand-audctx_' num2str(grid_spacing_mm) 'mm'];
    fla_output_fname = [...
        'smooth' num2str(100*fwhm, '%.0f') ...
        'mm_grid_hand-audctx_' num2str(grid_spacing_mm) 'mm_10whitematterPCs'];
    
    % first level analysis script
    % returns a structure with the needed file
    S = fla_matlab(...
        localizer_info.exp, us, localizer_info.runtype, ...
        localizer_runs_to_use, 'MION_CUSTOM1', 'downsampled_surface',...
        fla_input_fname, fla_output_fname, ...
        'n_whitematter_PCs', 10, 'motcorr_highres_2mm',...
        'contrasts',{localizer_info.contrast},...
        'permute_condition_order', 'n_cpus', 1, 'n_perms', 10000,...
        varargin{:});
    clear fla_input_fname fla_output_fname;
    
    % p-file based on ols or permutation test
    if optInputs(varargin, 'perm_stats')
        
        xi = ismember(S.p_file_permuted_gaussfit(:,1), localizer_info.contrast);
        pstat_file = S.p_file_permuted_gaussfit{xi,2};
        clear xi;
        
    else
        
        xi = ismember(S.p_file(:,1), localizer_info.contrast);
        pstat_file = S.p_file{xi,2};
        clear xi;
        
    end
    
else % second level analysis
    
    grid_file = ['smooth' num2str(fwhm*100, '%.0f') 'mm_grid_' grid_roi ...
        '_' num2str(grid_spacing_mm) 'mm_10whitematterPCs'];
    
    if optInputs(varargin, 'perm_stats') % permutation test
        
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
        % when there are more than 50 runs, the string is written in a more compressed form
        if length(localizer_runs_to_use) > 50
            localizer_run_string = [num2str(length(localizer_runs_to_use)) 'r-' num2str(localizer_runs_to_use(1)) '-' num2str(localizer_runs_to_use(end)) '_' DataHash(localizer_runs_to_use)];
        else
            localizer_run_string = ['r' sprintf('%d',localizer_runs_to_use)];
        end
        
        % output directory name
        tla_directory_name = [...
            localizer_info.contrast '_' ...
            localizer_info.exp '_' num2str(us) '_' ...
            localizer_info.runtype '_' localizer_run_string];
        
        % call tla_permutation
        [~, ~, ~, pstat_file] = ...
            tla_permutation_voxelwise_stats(...
            P,'downsampled_surface',tla_directory_name,...
            'no_stat_plots', varargin{:});
        
    else % ols solution
        
        pstat_file = sla_matlab(...
            localizer_info.exp, us, ...
            localizer_info.runtype, localizer_info.contrast, ...
            localizer_runs_to_use,...
            'downsampled_surface',grid_file,grid_file,...
            'fixed_without_whitening',varargin{:});
        
    end
end

