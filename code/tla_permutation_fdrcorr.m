function mat_file_with_all_useful_statistics = tla_permutation_fdrcorr(P,volume_or_surface,tla_directory_name,varargin)

% function mat_file_with_all_useful_statistics = tla_permutation_fdrcorr(P,volume_or_surface,tla_directory_name,varargin)
% 
% Corrects voxelwise p-map created from a permutation test, 
% using standard voxelwise FDR (e.g. Genovese, Lazar, Nichols, 2002).
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
%     % P(i).lower_level_directory_name = 'smooth300mm_grid_hand-stp-stg_1.5mm';
%     P(i).lower_level_directory_name = 'smooth500mm_grid_hand-stp-stg_1.5mm_10whitematterPCs';
% end
% volume_or_surface = 'downsampled_surface';
% tla_permutation_voxelwise_stats(P,volume_or_surface,tla_directory_name,'n_smps',10e3);
% tla_permutation_fdrcorr(P,volume_or_surface,tla_directory_name,'n_smps',10e3,'sigmap','fdr_q',1.3);

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
        
        % p-values based on gaussian fits to the null distribution
        voxelwise_p_gaussfit_file =  [tla_directory  'pstat_permutation_' num2str(n_smps) '_gaussfit.mat'];
        voxelwise_p_gaussfit_rh_inflated_file =  [tla_directory  'rh.pstat_permutation_' num2str(n_smps) '_gaussfit.mgz'];
        voxelwise_p_gaussfit_lh_inflated_file =  [tla_directory  'lh.pstat_permutation_' num2str(n_smps) '_gaussfit.mgz'];
        
        mat_file_with_all_useful_statistics =  [tla_directory  'all_statistics_' num2str(n_smps) '.mat'];
        
        p_fdr_thresh_file =  [tla_directory  'pstat_permutation_' num2str(n_smps) '_fdr_corrected_voxP' num2str(voxel_pthresh) '_clustP' num2str(fdr_q)  '.mat'];
        p_rh_inflated_fdr_thresh_file =  [tla_directory  'rh.pstat_permutation_' num2str(n_smps) '_fdr_corrected_voxP' num2str(voxel_pthresh) '_clustP' num2str(fdr_q)  '.mgz'];
        p_lh_inflated_fdr_thresh_file =  [tla_directory  'lh.pstat_permutation_' num2str(n_smps) '_fdr_corrected_voxP' num2str(voxel_pthresh) '_clustP' num2str(fdr_q)  '.mgz'];
        
    otherwise
        error('Error in tla_permutation: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"...');
        
end

% voxel-wise pvalues
load(voxelwise_p_gaussfit_file);
p_values = [G.grid_data{1}(:); G.grid_data{2}(:)];

% remove NaNs and remove sign
p_values = abs(p_values(isfinite(p_values)));

% sort
p_values = sort(p_values,'descend');
V = length(p_values);
I = (1:V)';
above_thresh = p_values >= -log10(I/V) + fdr_q;
fdr_pvalue = p_values(find( above_thresh, 1, 'last'));

% save to grid
G.grid_data{1}(G.grid_data{1} <= fdr_pvalue) = 0;
G.grid_data{2}(G.grid_data{2} <= fdr_pvalue) = 0;
save(p_fdr_thresh_file, 'G');

% save to surface file, rh
x = MRIread(voxelwise_p_gaussfit_rh_inflated_file);
x.vol(x.vol<fdr_pvalue) = 0;
MRIwrite_surface(x.vol(:), p_rh_inflated_fdr_thresh_file, 'rh');

% save to surface file, lh
x = MRIread(voxelwise_p_gaussfit_lh_inflated_file);
x.vol(x.vol<fdr_pvalue) = 0;
MRIwrite_surface(x.vol(:), p_lh_inflated_fdr_thresh_file, 'lh');

% -- plot statistic --
if optInputs(varargin, 'sigmap')
    switch volume_or_surface
        case 'volume'
            
            error('Need to setup volume analysis');
            
        case 'surface'
            
            error('Need to setup surface analysis');
            
        case 'downsampled_surface'
            
            % plot p map on the downsampled surface
            load(p_fdr_thresh_file);
            figure;
            subplot(1,2,1);
            imagesc(flipud(rot90(G.grid_data{1})), [-6 6]);
            title('Right Hemi');
            subplot(1,2,2);
            imagesc(fliplr(flipud(rot90(G.grid_data{2}))), [-6 6]); %#ok<FLUDLR>
            title('Left Hemi');
            colorbar;
            
            bounds = [voxel_pthresh, 6];
            midpoint = bounds(1) + 0.5*(bounds(2) - bounds(1));
            overlay_threshold = [fdr_pvalue, midpoint, bounds(2)];
            
            % plot map on the surface
            freeview3('fsaverage','rh','overlay',p_rh_inflated_fdr_thresh_file, 'overlay_threshold', overlay_threshold);
            freeview3('fsaverage','lh','overlay',p_lh_inflated_fdr_thresh_file, 'overlay_threshold', overlay_threshold);
            
        otherwise
            error('Error in tla_permutation: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"');
    end
end