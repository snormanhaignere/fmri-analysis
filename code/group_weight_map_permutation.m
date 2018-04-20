function [logp_gaussfit, logp_counts, Z_subject_mean, Z_null_smps] = group_weight_map_permutation(X, Y, n_perms, plot_figures)

% function [logp_gaussfit, logp_counts, Z_subject_mean, Z_null_smps] = group_weight_map_permutation(X, Y, n_perms, plot_figures)
% 
% Regresses the matrix X against voxel responses in matrix Y, converts beta weights
% to z-values and then averages these across subjects. A permutation test across the rows
% of the regression matrix, X, is used to calculate voxel weights significance values.
% 
% -- Inputs --
% 
% X: [M x P] regression matrix, typically M stimuli by P predictors
% 
% Y: [M x N x S] voxel matrix, typically M stimuli by N voxels x S subjects
% 
% n_perms: number of permutations used to calculate null
% 
% plot_figures: whether or not to plot diagnostic figures (default is false)
% 
% -- Outputs --
% 
% logp_gaussfit: [N x P] matrix of significance values (-log10[p]), calculated using 
% a Gaussian approximation to the null distribution
% 
% logp_counts: [N x P] matrix of significance values (-log10[p]), calculated using 
% direct counts, maximum significance value equals -log10(1/n_perms)
% 
% Z_subject_mean: [N x P] matrix of Z-statistics averaged across subjects
% 
% Z_null_smps: [n_perms x N x P] matrix of samples from the null
% 
% -- Example: Run on data from the natural sound experiment --
% 
% % data matrix and ICA response profiles from the natural sound experiment
% load('/mindhive/nklab/u/svnh/naturalsound-analysis/data/voxel_matrix_165x4949_formatted_grid.mat','vgrid_scan_mean_selected','grid_x','group_voxel_selection');
% load('/mindhive/nklab/u/svnh/naturalsound-analysis/analyses/ICA_best_solution_reformatted.mat', 'R', 'component_names');
% 
% % run analysis
% Y = vgrid_scan_mean_selected;
% X = R;
% n_perms = 10;
% [logp_gaussfit, logp_counts, Z_subject_mean, Z_null_smps] = group_weight_map_permutation(X, Y, n_perms);

% % plot right hemisphere significance maps for each component
% for i = 1:size(X,2)
%     % map to a flattened 2D grid
%     n_grid_voxels_rh = numel(grid_x{1});
%     p_gaussfit_rh = nan([n_grid_voxels_rh,1]);
%     xi = group_voxel_selection(1:n_grid_voxels_rh);
%     n_selected_vox_rh = sum(xi);
%     p_gaussfit_rh(xi) = logp_gaussfit(1:n_selected_vox_rh,i);
%     p_gaussfit_grid = reshape(p_gaussfit_rh, size(grid_x{1}));
%     
%     % plot
%     subplot(2,3,i);
%     imagesc(flipud(rot90(p_gaussfit_grid)));
%     colormap(parula);
%     title(component_names{i});
%     colorbar;
% end

% path to directory with useful analysis scripts
addpath([params('rootdir') 'general-analysis-code'])

if nargin < 4
    plot_figures = false;
end

% dimensions
[M,P] = size(X);
[~,N,S] = size(Y);

% unwrap subjects, M x N*S
Y_unwrapped = reshape(Y, [M, N*S]);

% Z values from regression analysis, P x N*S
[B,SE] = regress_matrix(X,Y_unwrapped);
Z = B./SE;

% rewrape, N x S x P
Z_rewrap = reshape(Z', [N, S, P]);

% average across subjects, N x P
Z_subject_mean = squeeze(mean(Z_rewrap,2));

% null distribution
Z_null_smps = nan(n_perms, N, P);
for i = 1:P
    
    % separate out shared and unique variance
    X_single_column = X(:,i);
    X_leftout_column = X(:,setdiff(1:P,i));
    X_single_column_shared = X_leftout_column * pinv(X_leftout_column) * X_single_column;
    X_single_column_unique = X_single_column - X_single_column_shared;
    
    for k = 1:n_perms
        
        if mod(k,100) == 0
            fprintf('Completed %d permutations for column %d\n',k,i); drawnow;
        end
                
        % permuted the unique variance of a single column
        X_perm = X;
        X_perm(:,i) = X_single_column_shared + X_single_column_unique(randperm(M));
        
        % Z values from regression analysis
        % just retains single row for permuted regressors 1 x N*S
        [B,SE] = regress_matrix(X_perm,Y_unwrapped);
        Z_null = B(i,:)./SE(i,:);
        
        % rewrap, N x S
        Z_null_rewrap = reshape(Z_null', [N, S]);
        
        % average across subjects, N
        Z_null_smps(k,:,i) = squeeze(mean(Z_null_rewrap,2));
       
    end
end

% mean and standard deviation of the null samples, N x P
Z_null_mean = squeeze(mean(Z_null_smps));
Z_null_std = squeeze(std(Z_null_smps));

% convert z-statistic from the regression to a true z-statistic
% using the Gaussian fit to the null
% N x P
ZZ = (Z_subject_mean - Z_null_mean) ./ Z_null_std;

% convert z-statistic to p value
p_gaussfit = 2*normcdf(-abs(ZZ), 0, 1);
logp_gaussfit = -log10(p_gaussfit) .* sign(ZZ);

% significance values based on counting the number
% of times the null exceeds the true Z-statistic
% n_perms x N x P -> N x P
Z_subject_mean_replicate_for_each_perm = repmat(shiftdim(Z_subject_mean,-1),[n_perms,1,1]); % replicate for each permutation
p_counts_right_tail = squeeze(mean(Z_null_smps > Z_subject_mean_replicate_for_each_perm));
p_counts_left_tail = squeeze(mean(Z_null_smps < Z_subject_mean_replicate_for_each_perm));
p_counts_two_tail = 2*min(p_counts_right_tail,p_counts_left_tail);
p_counts_two_tail(p_counts_two_tail==0) = 1/n_perms;
logp_counts = -log10(p_counts_two_tail);
logp_counts(p_counts_left_tail < p_counts_right_tail) = -logp_counts(p_counts_left_tail < p_counts_right_tail);

if plot_figures
    figure;
    plot(logp_counts, logp_gaussfit,'o')
    xlim(-log10(1/n_perms)*[-1 1]); ylim(-log10(1/n_perms)*[-1 1]);
    xlabel('log10(p) from counts'); ylabel('log10(p) from Gaussian fits');
end
