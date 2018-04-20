function activation_size_fdrcorr = tla_permutation_bootstrapped_activation_size_fdrcorr(P,volume_or_surface,tla_directory_name,varargin)

% function clustsizes_bootstrapped_fdr = tla_permutation_bootstrapped_activation_size_fdrcorr(P,volume_or_surface,tla_directory_name,varargin)
% 
% Activation extent using FDR correction with bootstrapping to calculate error bars
%
% -- Example from Amusia experiment --
% clear P;
% usubs_amusia  = [45,49,51,53,55,57,59,71,73,75,171]';
% tla_directory_name = 'amusia_group_harm_vs_noise_5mm';
% for i = 1:11;
%     P(i).exp = 'amusia'; %#ok<*SAGROW>
%     P(i).us = usubs_amusia(i);
%     P(i).runtype = 'localizer';
%     P(i).runs = 1;
%     P(i).contrast = 'harm_vs_noise';
%     P(i).lower_level_directory_name = 'smooth500mm_grid_hand-stp-stg_1.5mm_10whitematterPCs';
% end
% volume_or_surface = 'downsampled_surface';
% tla_permutation_voxelwise_stats(P,volume_or_surface,tla_directory_name,'n_smps',10e3);
% activation_size_fdrcorr = tla_permutation_bootstrapped_activation_size_fdrcorr(P,volume_or_surface,tla_directory_name,'n_smps',10e3,'sigmap','fdr_q',0.05);

% scripts directories
source_directory = strrep(which('fla_matlab.m'),'fla_matlab.m','');
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
addpath([source_directory 'export_fig']);

% default number of samples used to compute p-values
n_smps = 10e3;
if optInputs(varargin, 'n_smps')
    n_smps = varargin{optInputs(varargin, 'n_smps')+1};
end

% p-value threshold, which determines the cluster size
voxel_pthresh = 3;
if optInputs(varargin, 'voxel_pthresh')
    voxel_pthresh = varargin{optInputs(varargin, 'voxel_pthresh') + 1};
end

% cluster-corrected p-value
fdr_q = 1.3;
if optInputs(varargin, 'fdr_q')
    fdr_q = varargin{optInputs(varargin, 'fdr_q') + 1};
end

% number of surface points on the fsaverage template brain
nsurfpts = 163842;

% files
switch volume_or_surface
    case 'volume'
        
        error('Need to setup volume analysis');
        
    case 'surface'
        
        error('Need to setup surface analysis');
        
    case 'downsampled_surface'
        
        tla_directory = [params('rootdir') 'freesurfer/fsaverage/tla_matlab/' tla_directory_name '_downsampled_hash' DataHash(P) '/'];
        if ~exist(tla_directory,'dir');
            mkdir(tla_directory);
        end
        
        % p-values based on counting the number of times the null
        % distribution exceeds the measured value
        voxelwise_p_file =  [tla_directory  'pstat_permutation_' num2str(n_smps) '.mat'];
        mat_file_with_all_useful_statistics =  [tla_directory  'all_statistics_' num2str(n_smps) '.mat'];
        
    otherwise
        error('Error in tla_permutation: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"...');
        
end

% check if cluster sizes have already been computed
X = whos('-file', mat_file_with_all_useful_statistics);
activation_sizes_already_computed = false;
for i = 1:length(X)
    if strcmp(X(i).name, 'activation_size_fdrcorr');
        activation_sizes_already_computed = true;
    end
end

if ~activation_sizes_already_computed || optInputs(varargin, 'overwrite') % if not, calculate via bootstrap
    
    % stats from previous analyses
    load(mat_file_with_all_useful_statistics, 'lower_level_z', 'nullmean_voxels_with_no_NaNs', 'nullstd_voxels_with_no_NaNs', 'voxels_with_no_NaNs');
    
    % initialize the 2D image for right and left hemispheres
    load(voxelwise_p_file);
    rh_grid = nan(size(G.grid_data{1}));
    lh_grid = nan(size(G.grid_data{2}));
   
    tic;
    n_subjects = length(P);
    subject_smps = randi(n_subjects, [n_subjects, n_smps]); % subjects sampled with replacement
    activation_size_fdrcorr = zeros(n_smps,2);
    for i = 1:n_smps
        
        % bootstrapped p-value map based on gaussian fits to the null
        bootstrapped_zstat = mean(lower_level_z(subject_smps(:,i), voxels_with_no_NaNs));
        x = (bootstrapped_zstat - nullmean_voxels_with_no_NaNs) ./ nullstd_voxels_with_no_NaNs;
        bootstrapped_pmap_gaussfit = nan(1, length(voxels_with_no_NaNs));
        bootstrapped_pmap_gaussfit(voxels_with_no_NaNs) = -log10(2*normcdf(-abs(x), 0, 1)) .* sign(x);
        
        % grid z-map
        rh_grid(:) = bootstrapped_pmap_gaussfit(1:numel(rh_grid));
        lh_grid(:) = bootstrapped_pmap_gaussfit(numel(rh_grid)+1:end);
                
        % remove NaNs and remove sign
        p_values = [rh_grid(:); rh_grid(:)];
        p_values = abs(p_values(isfinite(p_values)));
        
        % sort
        p_values = sort(p_values,'descend');
        V = length(p_values);
        I = (1:V)';
        above_thresh = p_values >= -log10(I/V) + fdr_q;
        fdr_pvalue = p_values(find( above_thresh, 1, 'last'));
        
        activation_size_fdrcorr(i,:) = [sum(rh_grid(:) > fdr_pvalue), sum(lh_grid(:) > fdr_pvalue)];
    end
    toc;
    
    save(mat_file_with_all_useful_statistics, '-append', 'activation_size_fdrcorr');
    
else
    load(mat_file_with_all_useful_statistics, 'activation_size_fdrcorr');
end
