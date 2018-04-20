function estimate_hrf_surf_downsampled_surface(exp,usubs,runtype,input_fname,hrf_type,hrf_name,win,sr,varargin)
% function estimate_hrf_surf_downsampled_surface(exp,usubs,runtype,input_fname,roi,hrf_type,win,hrf_name,win,sr,varargin)
% 
% Estimates the HRF using an FIR model applied to responses measured on the
% downsampled cortical surface in Freesurfer.
% 
% The HRF can either be an "onset_response function", which estimates the response
% to each block, or an "impulse_response_function" HRF, which estimates the response
% to an impulse.
% 
% usubs = [45,49,51,53,55,57,59,71,73,75,171];%,48,56,46,52,54,172,50,58,74,72,70];
% usubs = [48,56,46,52,54,172,50,58,74,72,70];
% exp = 'amusia';
% input_fname = 'motcorr_smooth300mm_grid_hand-stp-stg_1.5mm';
% runtype = 'localizer';
% hrf_type = 'onset_response_function';
% hrf_name = 'amusia-controls_onset_response_function';
% exvar_threshold = 0.1;
% varargin = {};
% estimate_hrf_surf_downsampled_surface(exp,usubs,runtype,input_fname,hrf_type,hrf_name,varargin)

addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
addpath(genpath('export_fig'));

% sampling rate and hrf window
switch hrf_type
    case {'impulse_response_function'}
        sr = 1;
        win = (1/sr):(1/sr):12;
    case {'onset_response_function'}
        blockdur = read_scanparams(exp,usubs(1),runtype);
        sr = 1;
        win = (1/sr):(1/sr):blockdur+16;
end

% voxels with the highest explained variance are selected
% this threshold determines the proportion of voxels taken
exvar_threshold = 0.1;

hrf_all_subjects = nan(length(win), length(usubs));
for j = 1:length(usubs) % loop through subjects
  
  runs = read_runs(exp,usubs(j),runtype,varargin{:});
  hrf_all_runs = nan(length(win), length(runs));
  for k = 1:length(runs) % loop through runs
    
    fprintf('us %d\n',usubs(j));
    
    % timing information
    conds = read_conditions(exp,usubs(j),runtype,'run',runs(k),varargin{:});
    b = read_timings(exp,usubs(j),runtype,runs(k),varargin{:});
    TR = read_scanparams(exp,usubs(j),runtype,'TR',varargin{:},'run',runs(k));
    nTR = read_scanparams(exp,usubs(j),runtype,'nTR',varargin{:},'run',runs(k));
    
    %% Design Matrix
     
    % condition_vector has numbers which indicating the "condition" of each time point
    % block_vector has numbers which indicate the block number of each time point
    % relative to the start of the run (e.g. 1,2,3,4)    totaltime = (nTR-1)*TR;
    totaltime = (nTR-1)*TR;
    nsmps = floor(totaltime*sr)+1;
    condition_vector = zeros(nsmps,1);
    block_vector = zeros(nsmps,1);
    t = (0:nsmps-1)/sr;
    for i = 1:length(b.onsets)
        % time points of this block
        xi = t >= b.onsets(i) - (1/sr*2) & t < b.onsets(i) + b.durs(i) - (1/sr*2);

        % block number
        block_vector(xi) = i;
        
        % condition number
        if ~strcmp(b.conds{i}, 'NULL')
            condition_vector(xi) = find(strcmp(b.conds{i}, conds));
        end
    end
    
    % candlestick/boxcar regressors
    switch hrf_type
      case 'onset_response_function'
        seq_binary = zeros(size(condition_vector));
        seq_binary(1) = condition_vector(1)~=0;
        seq_binary(2:end) = condition_vector(2:end)~=0 & block_vector(2:end)~=block_vector(1:end-1);
      case 'impulse_response_function'
        seq_binary = condition_vector~=0;
      otherwise
        error('hrf_type does not match options.')
    end
    
    % fir model regressors, demeaned
    fir_model = zeros(length(seq_binary),length(win));
    for i = 1:length(win)
      fir_model(i+1:end,i) = seq_binary(1:end-i);
      fir_model(:,i) = fir_model(:,i) - mean(fir_model(:,i));
    end
    
    %% Nuissance regressors (mainly white-matter PCs)
    
    % directory
    whitematter_PCs_directory = [params('rootdir') exp '/analysis/whitematter_PCs/'];
    if ~exist(whitematter_PCs_directory,'dir')
        mkdir(whitematter_PCs_directory)
    end
    
    whitematter_mat_file = [whitematter_PCs_directory 'us' num2str(usubs(j)) '_' runtype '_r' num2str(runs(k)) '.mat'];
    if ~exist(whitematter_mat_file,'file')
        
        % version of fsl and freesurfer to use for this experiment
        fsl_version = read_fsl_version(exp);
        freesurfer_version = read_freesurfer_version(exp);
        
        % directory with the highres structural
        struct_directory = [params('rootdir') exp '/analysis/preprocess/usub' num2str(usubs(j)) '/struct_r1/'];
        
        % white-matter mask, registered to the highres structural
        register_whitematter_highres_2mmAll(exp, usubs(j));
        wm_highres_2mm = [struct_directory 'white_matter_2mm_highres.nii.gz'];
        if ~exist(wm_highres_2mm,'file')
            
            % freesurfer volume
            orig = [params('rootdir') 'freesurfer/' exp '_us' num2str(usubs(j)) '/mri/orig.mgz']; % orig is the anatomical volume in the freesurfer format
            if ~exist(orig, 'file')
                orig = [params('rootdir') 'freesurfer/' exp '_us' num2str(usubs(j)) '/mri/brain.mgz']; % orig is the anatomical volume in the freesurfer format
                fprintf('Using brain.mgz instead of orig.mgz\n'); drawnow;
            end
            
            % highres volume
            highres = [struct_directory 'brain.nii.gz'];
            highres_2mm = [struct_directory 'brain_2mm.nii.gz'];
            if ~exist(highres_2mm,'file')
                unix_fsl(fsl_version, ['fslmaths ' highres ' -subsamp2 ' highres_2mm]);
            end
            
            % orig to highres, based on headers
            orig2highres_fsl = [struct_directory 'orig2highres.mat'];
            orig2highres_freesurf = [struct_directory 'orig2highres.dat'];
            unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' exp '_us' num2str(usubs(j)) ' --mov ' orig ' --targ ' highres_2mm ' --reg ' orig2highres_freesurf ' --fslregout ' orig2highres_fsl ' --regheader --noedit ']);
            
            % resample to highres 2mm
            wm_orig = [params('rootdir') 'freesurfer/' exp '_us' num2str(usubs(j)) '/mri/wm_thresh.nii.gz'];
            unix_fsl(fsl_version, ['flirt -interp trilinear -in ' wm_orig ' -ref ' highres_2mm ' -applyxfm -init ' orig2highres_fsl ' -out ' wm_highres_2mm]);
        end
        
        % white matter mask
        wm_br = MRIread(wm_highres_2mm);
        wm_mask = wm_br.vol(:)' > 0.5;
        
        % functional volume registered to the highres anatomical
        register_func2highres_apply(exp, usubs(j), runtype, runs(k), 'bbreg', 'motcorr');
        volumefile = 'motcorr_highres_2mm';
        vol_fname = [params('rootdir') exp '/analysis/preprocess/usub' num2str(usubs(j)) '/' runtype '_r' num2str(runs(k)) '/' volumefile '.nii.gz'];
        vol_br = MRIread(vol_fname);
        vol_dims = size(vol_br.vol);
        vol_data = reshape(vol_br.vol, [prod(vol_dims(1:3)), vol_dims(4)])';
        
        % select non-zero, non-NaN voxels within the white-matter mask
        wm_matrix = vol_data(:, wm_mask & all(vol_data>max(vol_data(:))/1e6) & all(~isnan(vol_data)));
        
        % demean voxels and perform SVD
        wm_matrix = wm_matrix - ones(size(wm_matrix,1),1)*mean(wm_matrix);
        [U,S] = svd(wm_matrix, 'econ');
        exvar = cumsum(diag(S).^2 / sum(diag(S).^2));
        
        save(whitematter_mat_file, 'U', 'exvar');
        
    else
        
        load(whitematter_mat_file, 'U', 'exvar');
        
    end
    
    whitematter_figures_directory = [params('rootdir') exp '/figures/whitematter_PCs/'];
    if ~exist(whitematter_figures_directory,'dir')
        mkdir(whitematter_figures_directory)
    end
    
    % number of white-matter PCs to use as regressors
    n_whitematter_PCs = 10;
    
    % save top N PCs to X_nuissance matrix
    X_nuissance = U(:,1:n_whitematter_PCs);
    
    % add linear nuissance term
    linear_drift = (1:nTR)';
    linear_drift = zscore(linear_drift);
    X_nuissance = [X_nuissance, linear_drift];
    
    % upsample nuissance regressors
    X_nuissance_upsamp = interp1((0:nTR-1)'*TR, X_nuissance, t, 'cubic');

    %% Load and format data

    % read data from mat file
    subjid = [exp '_us' num2str(usubs(j))];
    preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(runs(k)) '/'];
    load([preprocdir input_fname '.mat']);

    % reformat data
    rh_dims = size(G.grid_data{1});
    lh_dims = size(G.grid_data{2});
    data = [reshape(G.grid_data{1}(:), [rh_dims(1)*rh_dims(2), rh_dims(3)])', reshape(G.grid_data{2}(:), [lh_dims(1)*lh_dims(2), lh_dims(3)])'];
    
    % select voxels without NaN values
    voxels_without_NaNs = find(all(~isnan(data)));
    data = data(:, voxels_without_NaNs);
    n_voxels = length(voxels_without_NaNs);
    
    % upsample the data
    data_upsamp = interp1((0:nTR-1)'*TR, data, t, 'cubic');
    
    % normalize voxel timecourses to have mean 100
    for i = 1:n_voxels
        data_upsamp(:,i) = 100*data_upsamp(:,i)/mean(data_upsamp(:,i));
    end

    % fit model
    X = [fir_model, ones(length(t), 1), X_nuissance_upsamp];
    B = pinv(X)*data_upsamp;
    B = B(1:size(fir_model,2),:);
    
    % measure of explained variance
    resid = data_upsamp-fir_model*B;
    exvar = (var(data_upsamp) - var(resid))./var(data_upsamp);
    
    % plot
    figure(1);
    x = nan(1,prod(rh_dims(1:2)) + prod(lh_dims(1:2)));
    x(voxels_without_NaNs) = exvar;
    rh_grid = reshape(x(1:prod(rh_dims(1:2))), rh_dims(1:2));
    lh_grid = reshape(x(prod(rh_dims(1:2)) + (1:prod(lh_dims(1:2)))), lh_dims(1:2));
    subplot(1,2,1);
    imagesc(flipud(rot90(rh_grid)));
    title('Right Hemi');
    subplot(1,2,2);
    imagesc(fliplr(flipud(rot90(lh_grid)))); %#ok<FLUDLR>
    title('Left Hemi');
    colorbar;
    drawnow;
    
    % number of voxels to select to compute hrf
    n_voxels_to_select = round(exvar_threshold*length(exvar));
    
    % voxels with the highest explained variance
    [~,xi] = sort(exvar,'descend');
    voxels_to_select = xi(1:n_voxels_to_select);
    
    % estimate hrf by averaging across selected voxels
    hrf = mean(B(:,voxels_to_select),2);
    hrf_all_runs(:,k) = hrf;
    
    % plot
    figure(2);
    plot(win,hrf);
    title(sprintf('us %d, run %d\n',usubs(j),runs(k)));
    drawnow;    
    
  end
  
  hrf_all_subjects(:,j) = mean(hrf_all_runs,2);
  
end

%% Save HRF data

% normalize the hrfs for each subject
hrf_all_subjects_norm = hrf_all_subjects ./ (ones(size(hrf_all_subjects,1),1)*max(abs(hrf_all_subjects)));

% average across subjects and compute standard errors
hrf_mean = mean(hrf_all_subjects_norm,2);
hrf_sem = stderr_withsub_corrected(hrf_all_subjects_norm')';

% plot average hrf
figure;
errorbar(win,hrf_mean,hrf_sem);
xlabel('Time (s)'); ylabel('Response (au)');
box off;
export_fig([params('rootdir') 'custom_hrfs/' hrf_name '.pdf'],'-nocrop','-transparent','-pdf');

figure;
% colormap(distinguishable_colors(10,1));
plot(win,hrf_all_subjects_norm);
xlabel('Time (s)'); ylabel('Response (au)');
box off;
export_fig([params('rootdir') 'custom_hrfs/' hrf_name '_allsubjects.pdf'],'-nocrop','-transparent','-pdf');

% save hrf
hrf = hrf_mean/max(abs(hrf_mean));
save([params('rootdir') 'custom_hrfs/' hrf_name '.mat'],'hrf_all_subjects','hrf','win','usubs');

