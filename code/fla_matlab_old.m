function fla_matlab(exp,us,runtype,r,hrf_name,volume_or_surface,input_fname,output_directory_name,varargin)
% Runs a standard first-level analysis in matlab
% 
% Modified on 1/9/2015 to include a permutation test
% 
% Last Modified by Sam Norman-Haignere on 1/17/2015

sr = 5;
% directory containing this script
source_directory = strrep(which('sla_permutation.m'),'sla_permutation.m','');

% directory of useful matlab scripts for dealing with Freesurfer data
addpath([source_directory 'fs/']);

% useful scripts for saving pdfs of figures
addpath([source_directory 'export_fig']);

% version of fsl and freesurfer to use
fsl_version = read_fsl_version(exp, varargin{:});
freesurfer_version = read_freesurfer_version(exp, varargin{:});

% filter cutoff
cutoff_sec = inf;
if optInputs(varargin, 'cutoff_sec')
    cutoff_sec = varargin{optInputs(varargin, 'cutoff_sec')+1};
end

% time-point exclusion cutoff based on global derivative
tp_exclusion_cutoff = inf;
if optInputs(varargin, 'tp_exclusion_cutoff') && optInputs(varargin, 'global_derivative')
    tp_exclusion_cutoff = varargin{optInputs(varargin, 'tp_exclusion_cutoff')+1};
end

% whitematter PCs to regress out
n_whitematter_PCs = 0;
if optInputs(varargin, 'n_whitematter_PCs')
    n_whitematter_PCs = varargin{optInputs(varargin, 'n_whitematter_PCs')+1};
end

% default number of permutations to use for permutation testing
n_perms = 100;
if optInputs(varargin, 'n_perms')
    n_perms = 100;
end

% number of cores to use
% additional cores can be useful for permutation analyses
n_cpus = 1;
if optInputs(varargin, 'n_cpus') 
    n_cpus = varargin{optInputs(varargin, 'n_cpus')+1};
end

% open a parallel pool if n_cpus > 1
if n_cpus > 1
    if exist('parpool','file')
        delete(gcp('nocreate'));
        max_cpus = n_cpus;
        while 1
            try
                parobj = parpool(n_cpus);
            catch
                if n_cpus < round(max_cpus/2)
                    error('Unable to start a pool with at least %d CPUs\n', n_cpus);
                else
                    n_cpus = n_cpus-1;
                    sprintf('Reducing requested number of CPUs to %d and trying again\n', n_cpus);
                end
            end
        end
    elseif exist('parcluster','file')
        pclust = parcluster;
        pclust.NumWorkers = n_cpus;
        if matlabpool('size') < n_cpus
            if matlabpool('size') > 0
                matlabpool('close');
            end
            pclust.matlabpool;
        end
    else
        if matlabpool('size') < n_cpus
            if matlabpool('size') > 0
                matlabpool('close');
            end
            matlabpool(n_cpus);
        end
    end
end

%% Creat design matrix

% condition and timing information
conds = read_conditions(exp,us,runtype,'run',r,varargin{:});
b = read_timings(exp,us,runtype,r,varargin{:});
TR = read_scanparams(exp,us,runtype,'run',r,'TR',varargin{:});
nTR = read_scanparams(exp,us,runtype,'run',r,'nTR',varargin{:});

% condition-indicator vector for upsampled timecourses
totaltime = (nTR-1)*TR;
nsmps = floor(totaltime*sr)+1;
seq_upsampled = zeros(nsmps,1);
t = (0:nsmps-1)/sr;
for i = 1:length(b.onsets)
    if ~strcmp(b.conds{i}, 'NULL')
        xi = t >= b.onsets(i) & t < b.onsets(i)+b.durs(i);
        seq_upsampled(xi) = find(strcmp(b.conds{i}, conds));
    end
end

% hrf
switch hrf_name
    case {'BOLD','MION','MION_CUSTOM1'};
        hrf_type = 'impulse_response_function';
        h = hrf_fsfast_gamma(1/sr, hrf_name, 'noplot');
        h = h/abs(sum(h));
    case {'naturalsound_block_response_function','naturalsound_monkey_block_response_function_lowpass'};
        hrf_type = 'onset_response_function';
        x = load([params('rootdir') 'custom_hrfs/' hrf_name '.mat']);
        h = interp1(x.win,x.hrf,x.win(1):(1/sr):x.win(end));
        h = h/max(h);
    otherwise
        error('Error in fla_matlab.m: ')
end

% convert indicator vector to onset vector if using "onset response function"
switch hrf_type
    case 'onset_response_function'
        block_onsets = [seq_upsampled(1)~=0; seq_upsampled(2:end)~=0 & seq_upsampled(2:end)~=seq_upsampled(1:end-1)];
        seq_upsampled(~block_onsets) = 0;
    case 'impulse_response_function'
    otherwise
        error('hrf_type does not match options.')
end

% convolve with hrf, and downsample to fMRI sampling rate
model = zeros(nsmps,length(conds));
for i = 1:length(conds)
    boxcar = zeros(nsmps,1);
    boxcar(seq_upsampled==i) = 1;
    
    x = conv(boxcar,h);
    model(:,i) = x(1:nsmps);
end
model_scans = interp1(t', model, (0:nTR-1)'*TR);

% permutation regression matrices
if optInputs(varargin, 'permute_condition_order')
    
    fprintf('Computing permuted regressors\n');
    tic;

    % condition vector for all permutations
    permuted_seq_upsampled = zeros(nsmps,n_perms);
    for q = 1:n_perms
        conds_permuted = b.conds(randperm(length(b.conds)));
        for i = 1:length(b.onsets)
            if ~strcmp(conds_permuted{i}, 'NULL')
                xi = t >= b.onsets(i) & t < b.onsets(i)+b.durs(i);
                permuted_seq_upsampled(xi,q) = find(strcmp(conds_permuted{i}, conds));
            end
        end
    end

    % convolve with hrf, and downsample to fMRI sampling rate
    % performed in parallel for each permutation
    permuted_model_scans = nan(nTR,length(conds),n_perms);
    parfor q = 1:n_perms
        permuted_model = zeros(nsmps,length(conds));
        for i = 1:length(conds)
            boxcar = zeros(nsmps,1);
            boxcar(permuted_seq_upsampled(:,q)==i) = 1;
            x = conv(boxcar,h);
            permuted_model(:,i) = x(1:nsmps);
        end
        permuted_model_scans(:,:,q) = interp1(t', permuted_model, (0:nTR-1)'*TR);
    end
    toc;
end

%% Regression analysis

% output files
switch volume_or_surface
    case 'volume'
        glmdir = [params('rootdir') exp '/analysis/fla_matlab/usub' num2str(us) '/' runtype '_r' num2str(r) '/' output_directory_name '/'];
        beta_file = [glmdir 'betas.nii.gz'];
        psc_file = [glmdir 'psc.nii.gz'];
        rvar_file = [glmdir 'rvar.nii.gz'];
        beta_white_file = [glmdir 'betas_white.nii.gz'];
        psc_white_file = [glmdir 'psc_white.nii.gz'];
        rvar_white_file = [glmdir 'rvar_white.nii.gz'];
        design_file = [glmdir 'design.mat'];
    case 'surface'
        subjid = [exp '_us' num2str(us)];
        if optInputs(varargin, 'monkey')
            glmdir = [params('rootdir') 'freesurfer/' subjid '/fla_matlab/' runtype '_r' num2str(r) '/' output_directory_name '/'];
        else
            glmdir = [params('rootdir') 'freesurfer/fsaverage/fla_matlab/' subjid '/' runtype '_r' num2str(r) '/' output_directory_name '/'];
        end
        if ~optInputs(varargin, 'hemi')
            error('Error in fit_glm: Need to specify hemisphere as optional argument for surface data.m');
        end
        hemi = varargin{optInputs(varargin, 'hemi')+1};
        beta_file = [glmdir hemi '.betas.mgz'];
        psc_file = [glmdir hemi '.psc.mgz'];
        rvar_file = [glmdir hemi '.rvar.mgz'];
        beta_white_file = [glmdir hemi '.betas_white.mgz'];
        psc_white_file = [glmdir hemi '.psc_white.mgz'];
        rvar_white_file = [glmdir hemi '.rvar_white.mgz'];
        design_file = [glmdir hemi '.design.mat'];
    otherwise
        error('Error in sla_matlab: volume_or_surface flag must be either "volume" or "surface"...');
end

if optInputs(varargin, 'permute_condition_order')
    design_file = strrep(design_file, '.mat', ['_' num2str(n_perms) '_permuted_conds.mat']);
end

% create glm directory
if ~exist(glmdir,'dir')
    mkdir(glmdir);
end

if  ~exist(design_file,'file') || ~exist(beta_file,'file') || ~exist(psc_file,'file') || ~exist(rvar_file,'file') || ~exist(beta_white_file,'file') || ~exist(psc_white_file,'file') || ~exist(rvar_white_file,'file') || optInputs(varargin, 'overwrite')
    
    % specify input file
    switch volume_or_surface
        case 'volume'
            preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
            inputfile = [preprocdir input_fname '.nii.gz'];
            copyfile([preprocdir 'example_func.nii.gz'],[glmdir 'example_func.nii.gz'],'f');
        case 'surface'
            subjid = [exp '_us' num2str(us)];
            if optInputs(varargin, 'monkey')
                preprocdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
            else
                preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
            end
            inputfile = [preprocdir hemi '.' input_fname '.mgz'];
    end
    
    try
        % read data and reformat to time x voxel matrix
        br = MRIread(inputfile);
        dims = size(br.vol);
        data = reshape(br.vol, [prod(dims(1:3)), dims(4)])';
    catch
        keyboard;
    end
    
    % select nonzero voxels and demean voxel timecourses
    n_voxels = size(data,2);
    nonzero_voxels = all(data>max(data(:))/1e6) & all(~isnan(data));
    n_nonzero_voxels = sum(nonzero_voxels);
    voxel_means = mean(data(:,nonzero_voxels));
    Y = data(:,nonzero_voxels) - ones(nTR,1)*voxel_means;
    
    % demean design matrix
    X_design = model_scans;
    X_design = X_design - ones(size(X_design,1),1)*mean(X_design);
    
    % demean permuted design matrices
    if optInputs(varargin, 'permute_condition_order')
        permuted_X_design = permuted_model_scans;
        permuted_X_design = permuted_X_design - repmat(mean(permuted_X_design), [nTR,1,1]);
    end
    
    % can optionally specify extra "nuissance" regressors
    X_nuissance = [];
    
    if n_whitematter_PCs > 0
        
        struct_directory = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/'];
        wm_highres_2mm = [struct_directory 'white_matter_2mm_highres.nii.gz'];
        if ~exist(wm_highres_2mm,'file')
                        
            % freesurfer volume
            orig = [params('rootdir') 'freesurfer/' exp '_us' num2str(us) '/mri/orig.mgz']; % orig is the anatomical volume in the freesurfer format
            if ~exist(orig, 'file')
                orig = [params('rootdir') 'freesurfer/' exp '_us' num2str(us) '/mri/brain.mgz']; % orig is the anatomical volume in the freesurfer format
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
            unix_freesurfer_version( freesurfer_version, ['tkregister2 --s ' exp '_us' num2str(us) ' --mov ' orig ' --targ ' highres_2mm ' --reg ' orig2highres_freesurf ' --fslregout ' orig2highres_fsl ' --regheader --noedit ']);
            
            % resample to highres 2mm
            wm_orig = [params('rootdir') 'freesurfer/' exp '_us' num2str(us) '/mri/wm_thresh.nii.gz'];
            unix_fsl(fsl_version, ['flirt -interp trilinear -in ' wm_orig ' -ref ' highres_2mm ' -applyxfm -init ' orig2highres_fsl ' -out ' wm_highres_2mm]);
        end
        
        % white matter mask
        wm_br = MRIread(wm_highres_2mm);
        wm_mask = wm_br.vol(:)' > 0.5;
            
        % if surface analysis need to load volume file
        if ~strcmp(volume_or_surface,'volume')
            % read in volume file in order to measure white matter voxels
            volumefile = varargin{optInputs(varargin,'n_whitematter_PCs')+2};
            vol_fname = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/' volumefile '.nii.gz'];
            vol_br = MRIread(vol_fname);
            vol_dims = size(vol_br.vol);
            vol_data = reshape(vol_br.vol, [prod(vol_dims(1:3)), vol_dims(4)])';
            wm_matrix = vol_data(:, wm_mask & all(vol_data>max(vol_data(:))/1e6) & all(~isnan(vol_data)));
        else
            wm_matrix = data(:, wm_mask & nonzero_voxels);
        end
                
        % demean voxels and perform SVD
        wm_matrix = wm_matrix - ones(size(wm_matrix,1),1)*mean(wm_matrix);
        [U,S] = svd(wm_matrix, 'econ');
        exvar = cumsum(diag(S).^2 / sum(diag(S).^2));
        
        % cumulative scree plot
        figure;
        plot(exvar,'k-o');
        hold on;
        plot(n_whitematter_PCs*[1 1],[0 1],'r--');
        xlabel('PC #'); ylabel('Cumulative Explained Variance'); ylim([0 1]);
        box off;
        % export_fig([glmdir 'whitematter_PC_explained_variance.pdf'],'-pdf','-transparent','-nocrop');
        
        % timecourses for first N PCs
        figure;
        legendlabels = cell(1,n_whitematter_PCs);
        cols = colormap(sprintf('jet(%d)',n_whitematter_PCs))*0.8;
        for i = 1:n_whitematter_PCs
            plot(U(:,i),'Color',cols(i,:));
            if i == 1
                hold on;
            end
            legendlabels{i} = sprintf('PC %d',i);
        end
        legend(legendlabels{:}); xlabel('TRs'); ylabel('Response');
        box off;
        % export_fig([glmdir 'whitematter_topPCs.pdf'],'-pdf','-transparent','-nocrop');
        
        %         sample_vox = squeeze(br.vol(44+1,35+1,58+1,:));
        %         sample_vox = sample_vox - mean(sample_vox);
        %         sample_vox_residual = sample_vox - U(:,1:n_whitematter_PCs)*pinv(U(:,1:n_whitematter_PCs))*sample_vox;
        X_nuissance = [X_nuissance, U(:,1:n_whitematter_PCs)];
        
    end
    

    % can optionally include "global derivative" as nuissance regressor
    % can also exclude time points based on global derivative
    if optInputs(varargin, 'global_derivative')
        
        % compute global derivative by convolving timecourses with laplacian
        Ypeak = nan(size(Y));
        for i = 1:size(Y,2)
            Ypeak(:,i) = conv(Y(:,i), [-1 2 -1]','same');
        end
        global_derivative = mean(abs(Ypeak),2);
        
        % normalize global derivative
        x = global_derivative;
        global_derivative_norm = (x-median(x))/(median(x)-min(x));
        
        % global derivative figures
        global_derivative_directory = [glmdir 'global_derivative_plots/'];
        if ~exist(global_derivative_directory,'dir')
            mkdir(global_derivative_directory);
        end
        figure;
        plot(global_derivative_norm);
        hold on;
        plot([1,length(global_derivative_norm)],1*[1 1],'k','LineWidth',2);
        plot([1,length(global_derivative_norm)],tp_exclusion_cutoff*[1 1],'r','LineWidth',2);
        xlabel('Scans'); ylabel('Normalized Global Derivative');
        box off;
        % output files
        switch volume_or_surface
            case 'volume'
                % export_fig([global_derivative_directory input_fname '_global_derivative_timecourse_' runtype '_r' num2str(r) '_' num2str(us) ifelse(exist('hemi','var'),'_hemi','') '.pdf'],'-pdf','-transparent','-nocrop');
            case 'surface'
                % export_fig([global_derivative_directory hemi '.' input_fname '_global_derivative_timecourse_' runtype '_r' num2str(r) '_' num2str(us) ifelse(exist('hemi','var'),'_hemi','') '.pdf'],'-pdf','-transparent','-nocrop');
        end
        
        % global derivative figures
        figure;
        hist(global_derivative_norm,length(global_derivative_norm)/20);
        hold on;
        yL = ylim;
        plot(tp_exclusion_cutoff*[1 1],yL,'r','LineWidth',2);
        plot(1*[1 1],yL,'k','LineWidth',2);
        plot(0*[1 1],yL,'k','LineWidth',2);
        xlabel('Normalized Global Derivative'); ylabel('Counts');
        box off;
        % output files
        switch volume_or_surface
            case 'volume'
                % export_fig([global_derivative_directory input_fname '_global_derivative_histogram_' runtype '_r' num2str(r) '_' num2str(us) ifelse(exist('hemi','var'),'_hemi','') '.pdf'],'-pdf','-transparent','-nocrop');
            case 'surface'
                % export_fig([global_derivative_directory hemi '.' input_fname '_global_derivative_histogram_' runtype '_r' num2str(r) '_' num2str(us) ifelse(exist('hemi','var'),'_hemi','') '.pdf'],'-pdf','-transparent','-nocrop');
        end
        
        
        % exclude voxels whose global derivative exceeds some threshold
        xi = find(global_derivative_norm > tp_exclusion_cutoff);
        exclusion_matrix = zeros(size(Y,1), length(xi));
        fprintf('%d time-points excluded (%.1f%%)\n',length(xi),100*length(xi)/length(global_derivative_norm));
        for j = 1:length(xi)
            exclusion_matrix(xi(j),j) = 1;
        end
        X_nuissance = [X_nuissance,global_derivative_norm,exclusion_matrix];
    end
    
    % optionally include polynomial regressors for each run
    if optInputs(varargin, 'polynomial')
        
        polyorder = varargin{optInputs(varargin, 'polynomial')+1};
        
        % if run is concatenated from other runs, can use this flag to specify
        % the onsets of the original runs
        if optInputs(varargin, 'concatenated_run')
            runtype_original = varargin{optInputs(varargin, 'concatenated_run')+1};
            r_original = varargin{optInputs(varargin, 'concatenated_run')+2};
            onset_index = 1;
            run_onsets = nan(1,length(r_original));
            for j = 1:length(r_original)
                [~, ~, ~, ~, ~, ~, ~, nTR_original, ~] = read_scanparams(exp,us,runtype_original,'run',r_original(j),varargin(:));
                run_onsets(j) = onset_index;
                onset_index = onset_index + nTR_original;
            end
        else
            run_onsets = 1;
        end
        
        % polynomial regressors, 1 per polynomial per run
        nruns = length(run_onsets);
        polynomial_regressors = zeros(nTR, polyorder+1, nruns);
        for j = 1:nruns
            if j == nruns;
                n = nTR-run_onsets(j)+1;
            else
                n = run_onsets(j+1)-run_onsets(j);
            end
            t = zscore((1:n)');
            t_matrix = (t*ones(1,polyorder+1)) .^ (ones(n,1)*(polyorder:-1:0));
            polynomial_regressors((1:n)+run_onsets(j)-1,:,j) = t_matrix;
        end
        
        % unwrap to nTR x nRegressors matrix, remove last boxcar to avoid
        % singular matrix
        dims = size(polynomial_regressors);
        polynomial_regressors_unwrapped = reshape(polynomial_regressors, dims(1), prod(dims(2:end)));
        polynomial_regressors_unwrapped = polynomial_regressors_unwrapped(:,1:prod(dims(2:end))-1);
        
        figure;
        plot(zscore(polynomial_regressors_unwrapped));
        xlabel('TR'); ylabel('Response');
        box off;
        % export_fig([glmdir 'polynomial_regressors.pdf'],'-pdf','-transparent','-nocrop');
        
        X_nuissance = [X_nuissance, zscore(polynomial_regressors_unwrapped)];
    end
    
    % optionally filter timecourses and design matrix
    if cutoff_sec ~= inf;
        
        cutoff_Hz = (1/cutoff_sec);
        nyq_Hz = (1/TR)/2;
        [B,A] = butter(4,cutoff_Hz/nyq_Hz,'high');
        
        % filter timecourses
        Yhp = zeros(size(Y));
        for i = 1:size(Y,2)
            Yhp(:,i) = filtfilt(B,A,Y(:,i));
        end
        
        % filter regressors
        Xhp = zeros(size(X_design));
        for i = 1:size(X_design,2)
            Xhp(:,i) = filtfilt(B,A,X_design(:,i));
        end
        
        if optInputs(varargin, 'permute_condition_order')
            permuted_Xhp = zeros(size(permuted_X_design));
            for i = 1:size(permuted_X_design,2)
                for j = 1:size(permuted_X_design,3)
                    permuted_Xhp(:,i,j) = filtfilt(B,A,permuted_X_design(:,i,j));
                end
            end
        end
        
        X_nuissance_hp = zeros(size(X_nuissance));
        for i = 1:size(X_nuissance,2)
            X_nuissance_hp(:,i) = filtfilt(B,A,X_nuissance(:,i));
        end
        
        X_design = Xhp;
        permuted_X_design = permuted_Xhp;
        Y = Yhp;
        X_nuissance = X_nuissance_hp;
        
    end
        
    
    %     vnums = find(nonzero_voxel_indices)-1;
    %     vras = nan(n_nonzero_voxels, 3);
    %     for i = 1:n_nonzero_voxels;
    %         xi = lb.vnums==vnums(i);
    %         vras(i,:) = lb.vras(xi,:);
    %     end
    
    % Xpad = [Xhp; zeros(size(Xhp))];
    
    %     % least-squares estimate without nuissance regressors
    %     X = X_design;
    %     B = pinv(X)*Y;
    %     E = Y - X*B;
    %     df = nTR-size(X,2);
    %     rvar = sum(E.^2)/df;
    %     Epad = [E; zeros(size(E))];
    %     acf_residual = real(ifft(abs(fft(Epad)).^2));
    %     acf_residual = acf_residual ./ (ones(nTR*2,1)*acf_residual(1,:));
    %     pspec_residual = real(fft(acf_residual));
    %
    %     % tukey-taper acf
    %     M = round(17/TR);
    %     maxlag = ceil((nTR*2+1)/2);
    %     tukey_window = (1/2)*(1+cos(pi*(0:M-1)/M))' * ones(1,n_nonzero_voxels);
    %     acf_onesided = acf_residual(1:maxlag,:);
    %     acf_tukey_onesided = acf_onesided;
    %     acf_tukey_onesided(M+1:end,:) = 0;
    %     acf_tukey_onesided(1:M,:) = tukey_window .* acf_onesided(1:M,:);
    %     acf_tukey_twosided = [acf_tukey_onesided; acf_tukey_onesided(end-1:-1:2,:)];
    %
    %     % invert acf in frequency domain
    %     pspec_tukey = real(fft(acf_tukey_twosided)); % real to deal with numeric error
    %     pspec_inv = 1./sqrt(pspec_tukey);
    %     K = real(ifft(pspec_inv));
    %     % acf_residual_inverse = real(ifft(abs(fft(K)).^2));
    %
    %     % prewhiten data
    %     FK = fft(K);
    %     Ypad = [Y; zeros(size(Y))];
    %     FYpad = fft(Ypad);
    %     Y_white = real(ifft(FYpad .* FK));
    %
    %     % prewhiten design matrix and estimate betas/residual
    %     Xpad = [X; zeros(size(X))];
    %     FXpad = fft(Xpad);
    %     p = size(Xpad,2);
    %     B_white = nan(p,n_nonzero_voxels);
    %     % rvar_white = nan(1,n_nonzero_voxels);
    %     norm_matrix = nan(p, p, n_nonzero_voxels);
    %     E_white = nan(nTR*2, n_nonzero_voxels);
    %     for i = 1:n_nonzero_voxels
    %         X_white = real(ifft(FXpad .* (FK(:,i)*ones(1,p))));
    %
    %         % estimate betas and residual
    %         pinv_X_white = pinv(X_white);
    %         B_white(:,i) = pinv_X_white*Y_white(:,i);
    %         E_white(:,i) = Y_white(:,i)-X_white*B_white(:,i);
    %         norm_matrix(:,:,i) = pinv_X_white*pinv_X_white';
    %     end
    %
    %     rvar_white = sum(E_white.^2)/df;
    %     acf_whitened_residual = real(ifft(abs(fft(E_white)).^2));
    %     pspec_whitened_residual = real(fft(acf_whitened_residual));
    %
    %
    %     cope_white = c*B_white;
    %     % cope_var_white = c*pinv(X_white)*Vi*pinv(X_white)'*c'*rvar_white_allvoxels; %#ok<*MINV>
    %     cope_var_white = nan(size(cope_white));
    %     for j = 1:n_nonzero_voxels
    %         cope_var_white(j) = c*norm_matrix(:,:,j)*c'*rvar_white(j);
    %     end
    %
    %     % cope_var_white_accurate = c*inv(Xpad'*Vi*Xpad)*c'*rvar_white_allvoxels;
    %
    %     % t-stat, df, and p-stat
    %     tstat_white = cope_white ./ sqrt(cope_var_white);
    %     % pstat_white = 2*tpvalue_copy(-abs(tstat_white), df);
    %     pstat_white = -sign(tstat_white).*log10(2*tpvalue_copy(-abs(tstat_white), df)); % two-tailed
    %     tstat_white(cope_var_white == 0) = 0;
    %     tstat_white(cope_var_white == 0) = 0;
        
    % least-squares estimate without nuissance regressors
    X = X_design;
    B = pinv(X)*Y;
    E = Y - X*B;
    df = nTR-size(X,2);
    rvar = sum(E.^2)/df;
    Epad = [E; zeros(size(E))];
    acf_residual = mean(real(ifft(abs(fft(Epad)).^2)),2); 
    acf_residual = acf_residual/acf_residual(1);
    pspec_residual = real(fft(acf_residual));
    
    % design matrix spectrum
    pspec_design = mean(abs(fft(X,nTR*2)).^2,2);
    maxlag = ceil((nTR*2+1)/2);
    f = 0.5*(1/TR)*(0:maxlag-1)/nTR;
    plot(f,pspec_design(1:maxlag)/max(pspec_design));
    xlabel('Frequency (Hz)'); ylabel('Power'); xlim([0 f(end)]);
    % export_fig([glmdir 'spectrum_design_matrix.pdf'],'-pdf','-nocrop','-transparent');
    
    % least-squares estimate with nuissance regressors if present
    if size(X_nuissance,2) > 0
        X = [X_design,zscore(X_nuissance)];
        B = pinv(X)*Y;
        E = Y - X*B;
        df = nTR-size(X,2);
        rvar = sum(E.^2)/df;
        Epad = [E; zeros(size(E))];
        acf_nuissance_residual = mean(real(ifft(abs(fft(Epad)).^2)),2);
        acf_nuissance_residual = acf_nuissance_residual/acf_nuissance_residual(1);
        pspec_nuissance_residual = real(fft(acf_nuissance_residual));
        
        figure;
        subplot(1,2,1);
        maxlag = ceil((nTR*2+1)/2);
        plot((0:maxlag-1)*TR, [acf_residual(1:maxlag)/acf_residual(1), acf_nuissance_residual(1:maxlag)/acf_nuissance_residual(1)]);
        xlabel('Lag (sec)'); ylabel('Power'); xlim([0 maxlag]);
        legend('Without Nuissance','With Nuissance','Location','Best');
        subplot(1,2,2);
        f = 0.5*(1/TR)*(0:maxlag-1)/nTR;
        plot(f,[pspec_residual(1:maxlag)/max(pspec_residual),  pspec_nuissance_residual(1:maxlag)/max(pspec_nuissance_residual)]);
        xlabel('Frequency (Hz)'); ylabel('Power'); xlim([0 f(end)]);
        legend('Without Nuissance','With Nuissance','Location','Best'); 
        % output files
        switch volume_or_surface
            case 'volume'
                % export_fig([glmdir 'acf_nuissance.pdf'],'-pdf','-nocrop','-transparent');
            case 'surface'
                % export_fig([glmdir hemi '.acf_nuissance.pdf'],'-pdf','-nocrop','-transparent');
        end
        
        % design matrix spectrum
        figure;
        pspec_design = mean(abs(fft(X_nuissance,nTR*2)).^2,2);
        maxlag = ceil((nTR*2+1)/2);
        f = 0.5*(1/TR)*(0:maxlag-1)/nTR;
        plot(f,pspec_design(1:maxlag)/max(pspec_design));
        xlabel('Frequency (Hz)'); ylabel('Power'); xlim([0 f(end)]);
        % export_fig([glmdir 'spectrum_nuissance.pdf'],'-pdf','-nocrop','-transparent');
        
        acf_residual = acf_nuissance_residual;
        pspec_residual = pspec_nuissance_residual;
    end
    
    % tukey-taper acf
    M = round(100/TR);
    maxlag = ceil((nTR*2+1)/2);
    tukey_window = (1/2)*(1+cos(pi*(0:M-1)/M))';
    acf_onesided = acf_residual(1:maxlag);
    acf_tukey_onesided = acf_onesided;
    acf_tukey_onesided(M+1:end) = 0;
    acf_tukey_onesided(1:M) = tukey_window .* acf_onesided(1:M);
    acf_tukey_twosided = [acf_tukey_onesided; acf_tukey_onesided(end-1:-1:2)]; 

    % invert acf in frequency domain
    pspec_tukey = real(fft(acf_tukey_twosided)); % real to deal with numeric error
    pspec_inv = 1./sqrt(pspec_tukey);
    K = real(ifft(pspec_inv));
    acf_residual_inverse = real(ifft(abs(fft(K)).^2));
    
    % prewhiten data
    FK = fft(K);
    Ypad = [Y; zeros(size(Y))];
    FYpad = fft(Ypad);
    Y_white = real(ifft(FYpad .* (FK*ones(1,n_nonzero_voxels))));

    % prewhiten design matrix
    Xpad = [X; zeros(size(X))];
    FXpad = fft(Xpad);
    p = size(Xpad,2);
    X_white = real(ifft(FXpad .* (FK*ones(1,p))));
    
    % estimate betas and residual
    B_white = pinv(X_white)*Y_white;
    E_white = Y_white-X_white*B_white;
    rvar_white = sum(E_white.^2)/df;
    acf_whitened_residual = mean(real(ifft(abs(fft(E_white)).^2)),2);
    pspec_whitened_residual = real(fft(acf_whitened_residual));
    
    % plot
    figure;
    subplot(1,2,1);
    maxlag = ceil((nTR*2+1)/2);
    plot((0:maxlag-1)*TR, [acf_residual(1:maxlag)/acf_residual(1), acf_whitened_residual(1:maxlag)/acf_whitened_residual(1)]);
    xlabel('Lag (sec)'); ylabel('Power'); xlim([0 maxlag]);
    legend('No Prewhitening','Prewhitening','Location','Best');
    subplot(1,2,2);
    f = 0.5*(1/TR)*(0:maxlag-1)/nTR;
    plot(f,[pspec_residual(1:maxlag)/max(pspec_residual),  pspec_whitened_residual(1:maxlag)/max(pspec_whitened_residual)]);
    xlabel('Frequency (Hz)'); ylabel('Power'); xlim([0 f(end)]);
    legend('No Prewhitening','Prewhitening','Location','Best');
    switch volume_or_surface
        case 'volume'
            % export_fig([glmdir 'acf_prewhitening.pdf'],'-pdf','-nocrop','-transparent');
        case 'surface'
            % export_fig([glmdir hemi '.acf_prewhitening.pdf'],'-pdf','-nocrop','-transparent');
    end
    

    p = size(X,2);
    B_allvoxels = zeros(p, n_voxels);
    B_allvoxels(:,nonzero_voxels) = B;
    
    p = size(X,2);
    B_white_allvoxels = zeros(p, n_voxels);
    B_white_allvoxels(:,nonzero_voxels) = B_white;
    
    rvar_allvoxels = zeros(1, n_voxels);
    rvar_allvoxels(:,nonzero_voxels) = rvar;
    
    rvar_white_allvoxels = zeros(1, n_voxels);
    rvar_white_allvoxels(:,nonzero_voxels) = rvar_white;
    
    psc = zeros(p, n_voxels);
    psc(:,nonzero_voxels) = B./ (ones(p,1)*voxel_means);
    psc(:,any(isnan(psc))) = 0;
    
    psc_white = zeros(p, n_voxels);
    psc_white(:,nonzero_voxels) = B_white./ (ones(p,1)*voxel_means);
    psc_white(:,any(isnan(psc_white))) = 0;
    
    %     norm_matrix_allvoxels = zeros(p, p, n_voxels);
    %     norm_matrix_allvoxels(:,:,nonzero_voxels) = norm_matrix;
    %
    % residual variance, and df
    % see Woolrich et al., 2001 and Worsley and Frison, 2005
    %     R = eye(size(E,1)) - Xnuissance*pinv(Xnuissance);
    %     rvar = zeros(1,size(data,2));
    %     rvar(nonzero_voxels) = sum(E.^2)/trace(R);
    %     df = trace(R*Ecorr).^2 / trace(R*Ecorr*R*Ecorr);
    
    % wrap data in old volume format
    dims_data = size(br.vol);
    B_wrap = reshape(B_allvoxels',[dims_data(1:3), size(B_allvoxels,1)]);
    psc_wrap = reshape(psc',[dims_data(1:3), size(psc,1)]);
    rvar_wrap = reshape(rvar_allvoxels',[dims_data(1:3),1]);
    B_white_wrap = reshape(B_white_allvoxels',[dims_data(1:3), size(B_white_allvoxels,1)]);
    psc_white_wrap = reshape(psc_white',[dims_data(1:3), size(psc_white,1)]);
    rvar_white_wrap = reshape(rvar_white_allvoxels',[dims_data(1:3),1]);
    
    % write betas
    br_beta = br;
    br_beta.vol = B_wrap;
    br_beta.nframes = size(br_beta.vol,4);
    br_beta.fspec = beta_file;
    MRIwrite(br_beta, beta_file);
    
    % write psc
    br_psc = br;
    br_psc.vol = psc_wrap;
    br_psc.nframes = size(br_psc.vol,4);
    br_psc.fspec = psc_file;
    MRIwrite(br_psc, psc_file);
    
    % variances
    br_rvar = br;
    br_rvar.vol = rvar_wrap;
    br_rvar.nframes = size(br_rvar.vol,4);
    br_rvar.fspec = rvar_file;
    MRIwrite(br_rvar, rvar_file);
    
    % write betas, prewhitened
    br_beta_white = br;
    br_beta_white.vol = B_white_wrap;
    br_beta_white.nframes = size(br_beta_white.vol,4);
    br_beta_white.fspec = beta_white_file;
    MRIwrite(br_beta_white, beta_white_file);
    
    % write psc, prewhitened
    br_white_psc = br;
    br_white_psc.vol = psc_white_wrap;
    br_white_psc.nframes = size(br_white_psc.vol,4);
    br_white_psc.fspec = psc_white_file;
    MRIwrite(br_white_psc, psc_white_file);
    
    % variances, prewhitened
    br_rvar_white = br;
    br_rvar_white.vol = rvar_white_wrap;
    br_rvar_white.nframes = size(br_rvar_white.vol,4);
    br_rvar_white.fspec = rvar_white_file;
    MRIwrite(br_rvar_white, rvar_white_file);
    
    if exist('permuted_X_design','var')
        save(design_file,'permuted_X_design','X_design','X_white','X_nuissance','X','acf_residual','acf_whitened_residual','acf_residual_inverse','df');
    else
        save(design_file,'X_design','X_nuissance','X','X_white','acf_residual','acf_whitened_residual','acf_residual_inverse','df');
    end
end


%% register statistics

if optInputs(varargin, 'reg2highres') && strcmp(volume_or_surface,'volume')
    
    % preprocessing directory
    if optInputs(varargin, 'featreg')
        regdir = [params('rootdir') exp '/analysis/fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_block_3mm.feat/reg/'];
    else
        regdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/reg_handtune/'];
    end
    highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz']; % brain extracted structural image
    highres_2mm = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain_2mm.nii.gz']; % brain extracted structural image
    if ~exist(highres_2mm,'file') || optInputs(varargin, 'overwrite')
        unix_fsl(fsl_version, ['fslmaths ' highres ' -subsamp2 ' highres_2mm]);
    end
    
    if ~exist([glmdir 'highres_2mm.nii.gz'],'file') || optInputs(varargin, 'overwrite')
        copyfile(highres_2mm, [glmdir 'highres_2mm.nii.gz']);
    end
    
    exfunc2highres_regmat = [regdir 'example_func2highres.mat']; % functional to anatomical registration matrix
    
    brainmask = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/brain_mask.nii.gz'];
    brainmask_2mm = [glmdir 'brain_mask_highres_2mm.nii.gz'];
    if ~exist(brainmask_2mm,'file') || optInputs(varargin, 'overwrite')
        unix_fsl(fsl_version, ['flirt -interp trilinear -in ' brainmask ' -ref ' highres_2mm ' -applyxfm -init ' exfunc2highres_regmat ' -out ' brainmask_2mm]);
        unix_fsl(fsl_version, ['fslmaths ' brainmask_2mm ' -thr 0.25 -bin ' brainmask_2mm]);
    end
    
    statistics_to_register = {psc_file, beta_file, psc_white_file, beta_white_file, rvar_file};
    for i = 1:length(statistics_to_register)
        statistic_highres = strrep(statistics_to_register{i}, '.nii.gz','_highres_2mm.nii.gz');
        if ~exist(statistic_highres,'file') || optInputs(varargin, 'overwrite')
            unix_fsl(fsl_version, ['flirt -interp trilinear -in ' statistics_to_register{i} ' -ref ' highres_2mm ' -applyxfm -init ' exfunc2highres_regmat ' -out ' statistic_highres]);
            unix_fsl(fsl_version, ['fslmaths ' statistic_highres ' -mas ' brainmask_2mm ' ' statistic_highres]);
        end
    end
end

if optInputs(varargin, 'nocontrasts')
    return;
end

%% setup contrasts

load(design_file);
evname = conds;
nevs = length(evname);
if optInputs(varargin, 'contrasts')
    contrastnames = varargin{optInputs(varargin,'contrasts')+1};
else
    contrastnames = read_contrasts(exp,us,runtype,model,'run',r,varargin{:});
end
contrasts = zeros(length(contrastnames),nevs);

% create contrast matrix
for c = 1:length(contrastnames);
    
    % an ev group name is a name for a set of conditions that can be used in a contrast
    evgroup = regexp(contrastnames{c},'_vs_','split');
    
    if length(evgroup)==1
        
        evinds = strmatchMANY( evgroup2evname(exp,us,evgroup{1},model,'run',r,varargin{:}),  evname );
        if isempty(evinds); error(['Failed to match contrasts to evs: ' contrastnames{c}]); end
        contrasts(c,evinds) = 1;
        
    else
        
        evinds1 = strmatchMANY( evgroup2evname(exp,us,evgroup{1},model,'run',r,varargin{:}),  evname );
        evinds2 = strmatchMANY( evgroup2evname(exp,us,evgroup{2},model,'run',r,varargin{:}),  evname );
        
        if isempty(evinds1) || isempty(evinds2); error(['Failed to match contrasts to evs: ' contrastnames{c}]); end
        
        contrasts(c,evinds1) = length(evinds2);
        contrasts(c,evinds2) = -length(evinds1);
        
    end
end

% "complicated" contrasts
[~,~,cc] = read_contrasts(exp,us,runtype,model,'run',r);
if ~isempty(cc)
    ccmat = nan(size(cc,1),nevs);
    
    for i = 1:size(cc,1)
        weightcell = cc{i,2};
        weightmat = zeros(1,nevs);
        for j = 1:size(weightcell,1)
            x = strcmp(weightcell{j,1},evname);
            weightmat(x) = weightmat(x) + weightcell{j,2};
        end
        ccmat(i,:) = weightmat;
    end
    
    contrasts = [contrasts; ccmat];
    contrastnames = [contrastnames, cc(:,1)'];
end

%% View beta/psc/variance maps using freeview

if strcmp(volume_or_surface,'surface') && optInputs(varargin, 'plot_betas')
    if optInputs(varargin, 'monkey')
        freeview3(subjid, hemi, 'overlay', beta_file);
    else
        freeview3('fsaverage', hemi, 'overlay', beta_file);
    end
end

if strcmp(volume_or_surface,'surface') && optInputs(varargin, 'plot_psc')
    if optInputs(varargin, 'monkey')
        freeview3(subjid, hemi, 'overlay', psc_file);
    else
        freeview3('fsaverage', hemi, 'overlay', psc_file);
    end
end

if strcmp(volume_or_surface,'surface') && optInputs(varargin, 'plot_rvar')
    if optInputs(varargin, 'monkey')
        freeview3(subjid, hemi, 'overlay', rvar_file);
    else
        freeview3('fsaverage', hemi, 'overlay', rvar_file);
    end
end

%% compute t statistic

statistics_to_register = {};
for i = 1:size(contrasts,1)
    
    contrast_directory = [glmdir 'contrasts/'];
    if ~exist(contrast_directory,'dir')
        mkdir(contrast_directory);
    end
    
    switch volume_or_surface
        case 'volume'
            t_file =  [contrast_directory 'tstat_' contrastnames{i} '.nii.gz'];
            p_file =  [contrast_directory 'pstat_' contrastnames{i} '.nii.gz'];
            cope_file =  [contrast_directory 'cope_' contrastnames{i} '.nii.gz'];
            cope_var_file =  [contrast_directory 'cope_var_' contrastnames{i} '.nii.gz'];
            t_white_file =  [contrast_directory 'tstat_white_' contrastnames{i} '.nii.gz'];
            p_white_file =  [contrast_directory 'pstat_white_' contrastnames{i} '.nii.gz'];
            cope_white_file =  [contrast_directory 'cope_white_' contrastnames{i} '.nii.gz'];
            cope_var_white_file =  [contrast_directory 'cope_var_white_' contrastnames{i} '.nii.gz'];
        case 'surface'
            t_file =  [contrast_directory hemi '.tstat_' contrastnames{i} '.mgz'];
            p_file =  [contrast_directory hemi '.pstat_' contrastnames{i} '.mgz'];
            cope_file =  [contrast_directory hemi '.cope_' contrastnames{i} '.mgz'];
            cope_var_file =  [contrast_directory hemi '.cope_var_' contrastnames{i} '.mgz'];
            t_white_file =  [contrast_directory hemi '.tstat_white_' contrastnames{i} '.mgz'];
            p_white_file =  [contrast_directory hemi '.pstat_white_' contrastnames{i} '.mgz'];
            cope_white_file =  [contrast_directory hemi '.cope_white_' contrastnames{i} '.mgz'];
            cope_var_white_file =  [contrast_directory hemi '.cope_var_white_' contrastnames{i} '.mgz'];
    end
    
    statistics_to_register = [statistics_to_register, t_file, p_file, cope_file, cope_var_file,  t_white_file, p_white_file, cope_white_file, cope_var_white_file]; %#ok<AGROW>

    if ~exist(t_file,'file') || ~exist(p_file,'file') || ~exist(cope_file,'file') || ~exist(cope_var_file,'file') || optInputs(varargin, 'overwrite')
        
        if ~exist('B_allvoxels','var')
            br_beta = MRIread(beta_file);
            dims = size(br_beta.vol);
            B_allvoxels = reshape(br_beta.vol, [prod(dims(1:3)), dims(4)])';
        end
        
        if ~exist('B_white_allvoxels','var')
            br_beta_white = MRIread(beta_white_file);
            dims = size(br_beta_white.vol);
            B_white_allvoxels = reshape(br_beta_white.vol, [prod(dims(1:3)), dims(4)])';
        end
        
        if ~exist('rvar_allvoxels','var')
            br_rvar = MRIread(rvar_file);
            rvar_allvoxels = br_rvar.vol(:)';
        end
        
        if ~exist('rvar_white_allvoxels','var')
            br_rvar_white = MRIread(rvar_white_file);
            rvar_white_allvoxels = br_rvar_white.vol(:)';
        end
                
        c = [contrasts(i,:), zeros(1,size(X_nuissance,2))];        
        
        % contrast magnitude and variance
        cope = c*B_allvoxels;
        cope_var = c*pinv(X)*pinv(X)'*c'*rvar_allvoxels;
        % cope_var = c*pinv(Xpad)*toeplitz(acf_residual)*pinv(Xpad)'*c'*rvar;
        
        % t-stat, df, and p-stat
        tstat = cope ./ sqrt(cope_var);
        % pstat = 2*tpvalue_copy(-abs(tstat), df);
        pstat = -sign(tstat).*log10(2*tpvalue_copy(-abs(tstat), df)); % two-tailed
        tstat(cope_var == 0) = 0;
        pstat(cope_var == 0) = 0;

        % Xpad = [X; zeros(size(X))];
        
        % contrast magnitude and variance for whitened data
        % Vi = toeplitz(acf_residual_inverse);
        cope_white = c*B_white_allvoxels;
        % cope_var_white = c*pinv(X_white)*Vi*pinv(X_white)'*c'*rvar_white_allvoxels; %#ok<*MINV>
        cope_var_white = c*pinv(X_white)*pinv(X_white)'*c'*rvar_white_allvoxels;
        % cope_var_white_accurate = c*inv(Xpad'*Vi*Xpad)*c'*rvar_white_allvoxels;
        
        % t-stat, df, and p-stat
        tstat_white = cope_white ./ sqrt(cope_var_white);
        % pstat_white = 2*tpvalue_copy(-abs(tstat_white), df);
        pstat_white = -sign(tstat_white).*log10(2*tpvalue_copy(-abs(tstat_white), df)); % two-tailed
        tstat_white(cope_var_white == 0) = 0;
        tstat_white(cope_var_white == 0) = 0;
        
        % t-stat, df, and p-stat
        %         tstat_white_accurate = cope_white ./ sqrt(cope_var_white_accurate);
        %         % pstat_white = 2*tpvalue_copy(-abs(tstat_white), df);
        %         pstat_white_accurate = -sign(tstat_white_accurate).*log10(2*tpvalue_copy(-abs(tstat_white_accurate), df)); % two-tailed
        %         pstat_white_accurate(cope_var_white == 0) = 0;
        %         pstat_white_accurate(cope_var_white == 0) = 0;
        
        % save degrees of freedom for higher-level analysis
        save([contrast_directory 'df.mat'],'df');        
        
        % template file
        br = MRIread(rvar_file);    
        
        % write files
        br_tstat = br;
        br_tstat.vol = reshape(tstat(:),size(br.vol));
        br_tstat.fspec = t_file;
        MRIwrite(br_tstat, t_file);
        
        br_pstat = br;
        br_pstat.vol = reshape(pstat(:),size(br.vol));
        br_pstat.fspec = p_file;
        MRIwrite(br_pstat, p_file);
        
        br_cope = br;
        br_cope.vol = reshape(cope(:),size(br.vol));
        br_cope.fspec = cope_file;
        MRIwrite(br_cope, cope_file);
        
        br_cope_var = br;
        br_cope_var.vol = reshape(cope_var(:),size(br.vol));
        br_cope_var.fspec = cope_var_file;
        MRIwrite(br_cope_var, cope_var_file);
        
        br_tstat_white = br;
        br_tstat_white.vol = reshape(tstat_white(:),size(br.vol));
        br_tstat_white.fspec = t_white_file;
        MRIwrite(br_tstat, t_white_file);
        
        br_pstat_white = br;
        br_pstat_white.vol = reshape(pstat_white(:),size(br.vol));
        br_pstat_white.fspec = p_white_file;
        MRIwrite(br_pstat_white, p_white_file);
        
        br_cope_white = br;
        br_cope_white.vol = reshape(cope_white(:),size(br.vol));
        br_cope_white.fspec = cope_white_file;
        MRIwrite(br_cope_white, cope_white_file);
        
        br_cope_var_white = br;
        br_cope_var_white.vol = reshape(cope_var_white(:),size(br.vol));
        br_cope_var_white.fspec = cope_var_white_file;
        MRIwrite(br_cope_var_white, cope_var_white_file);
    end
        
    if optInputs(varargin, 'permute_condition_order')
        switch volume_or_surface
            case 'volume'
                perm_cond_cope_file = [contrast_directory 'perm_condition_order_' num2str(n_perms) 'smps_' contrastnames{i} '.mat'];
            case 'surface'
                perm_cond_cope_file = [contrast_directory hemi '.perm_condition_order_' num2str(n_perms) 'smps_' contrastnames{i} '.mat'];
        end
    end
    
    if optInputs(varargin, 'permute_condition_order') && (~exist(perm_cond_cope_file,'file') || optInputs(varargin,'overwrite'))
                
        if ~exist('B_allvoxels','var')
            br_beta = MRIread(beta_file);
            dims = size(br_beta.vol);
            B_allvoxels = reshape(br_beta.vol, [prod(dims(1:3)), dims(4)])';
        end
       
        if ~exist('rvar_allvoxels','var')
            br_rvar = MRIread(rvar_file);
            rvar_allvoxels = br_rvar.vol(:)';
        end
        
        % GLM permutation analyses
        % not saved because it would consume considerable memory
        if ~exist('permuted_B','var') || ~exist('permuted_rvar','var')
            
            load(design_file);
            if ~exist('Y','var')
                % specify input file
                switch volume_or_surface
                    case 'volume'
                        preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
                        inputfile = [preprocdir input_fname '.nii.gz'];
                        copyfile([preprocdir 'example_func.nii.gz'],[glmdir 'example_func.nii.gz'],'f');
                    case 'surface'
                        subjid = [exp '_us' num2str(us)];
                        if optInputs(varargin, 'monkey')
                            preprocdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
                        else
                            preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
                        end
                        inputfile = [preprocdir hemi '.' input_fname '.mgz'];
                end
                try
                    % read data and reformat to time x voxel matrix
                    br = MRIread(inputfile);
                    dims = size(br.vol);
                    data = reshape(br.vol, [prod(dims(1:3)), dims(4)])';
                catch
                    keyboard;
                end
                % select nonzero voxels and demean voxel timecourses
                nonzero_voxels = all(data>max(data(:))/1e6) & all(~isnan(data));
                voxel_means = mean(data(:,nonzero_voxels));
                Y = data(:,nonzero_voxels) - ones(nTR,1)*voxel_means;
            end
            
            tic;
            fprintf('Computing betas and residual for design permutations.\n');
            permuted_B = nan(size(X_design,2), size(Y,2), n_perms);
            permuted_rvar = nan(size(Y,2),n_perms);
            df = nTR-size(X_design,2);
            parfor q = 1:n_perms
                tic;
                if size(X_nuissance,2) > 0
                    permuted_X = [permuted_X_design(:,:,q),zscore(X_nuissance)];
                else
                    permuted_X = permuted_X_design(:,:,q);
                end
                single_perm_B = pinv(permuted_X)*Y;
                permuted_B(:,:,q) = single_perm_B(1:size(X_design,2),:);
                permuted_E = Y - permuted_X*single_perm_B;
                permuted_rvar(:,q) = sum(permuted_E.^2)/df;
            end
            toc;
        end
        
        tic;
        fprintf('Contrast statistics for permutation test.\n'); drawnow;
        c = contrasts(i,:);
        c_with_nuissance = [contrasts(i,:), zeros(1,size(X_nuissance,2))];        
        permuted_cope = nan(n_perms, sum(nonzero_voxels));
        permuted_cope_var = nan(n_perms, sum(nonzero_voxels));
        parfor q = 1:n_perms
            if size(X_nuissance,2) > 0;
                permuted_X = [permuted_X_design(:,:,q),zscore(X_nuissance)];
            else
                permuted_X = permuted_X_design(:,:,q);
            end
            permuted_cope(q,:) = c*permuted_B(:,:,q);
            permuted_cope_var(q,:) = (c_with_nuissance*pinv(permuted_X)*pinv(permuted_X)'*c_with_nuissance'*permuted_rvar(:,q));
        end
        permuted_tstat = permuted_cope ./ sqrt(permuted_cope_var);
        
        % compare with measured tstat
        cope = c_with_nuissance*B_allvoxels(:,nonzero_voxels);
        cope_var = c_with_nuissance*pinv(X)*pinv(X)'*c_with_nuissance'*rvar_allvoxels(:,nonzero_voxels);
        tstat = cope ./ sqrt(cope_var);        
        permuted_pstat = mean(abs(permuted_tstat) > (ones(n_perms,1)*abs(tstat))); % two-tailed
        toc; 
        
        % save
        save(perm_cond_cope_file, 'permuted_cope','permuted_cope_var','permuted_tstat','permuted_pstat','nonzero_voxels');
        
        % plot useful statistics
        figure;
        subplot(1,2,1);
        x = linspace(0,1,100);
        nx = hist(permuted_pstat,x);
        plot(x,nx/sum(nx));
        xlabel('P-Value'); ylabel('Proportion of Voxels');
        title('PDF of Permutation-Based P-values');
        subplot(1,2,2);
        plot(sort(permuted_pstat), cumsum(ones(1,length(permuted_pstat)))/length(permuted_pstat), 'k-','LineWidth',2); hold on;
        plot([0 1],[0 1],'r--','LineWidth',2);
        xlabel('P-Value'); ylabel('Cumulative Proportion of Voxels');
        title('CDF of Permutation-Based P-values');
        switch volume_or_surface
            case 'volume'
                % export_fig([contrast_directory 'pstat_permuted_conditions_' contrastnames{i} '.pdf'],'-pdf','-transparent','-nocrop');
            case 'surface'
                % export_fig([contrast_directory hemi '.pstat_permuted_conditions_' contrastnames{i} '.pdf'],'-pdf','-transparent','-nocrop');
        end
        
    end
end

%% register contrasts
if optInputs(varargin, 'reg2highres') && strcmp(volume_or_surface,'volume')
    for i = 1:length(statistics_to_register)
        
        % preprocessing directory
        if optInputs(varargin, 'featreg')
            regdir = [params('rootdir') exp '/analysis/fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_block_3mm.feat/reg/'];
        else
            regdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/reg_handtune/'];
        end
        highres = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz']; % brain extracted structural image
        highres_2mm = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain_2mm.nii.gz']; % brain extracted structural image
        if ~exist(highres_2mm,'file') || optInputs(varargin, 'overwrite')
            unix_fsl(fsl_version, ['fslmaths ' highres ' -subsamp2 ' highres_2mm]);
        end
        
        if ~exist([glmdir 'highres_2mm.nii.gz'],'file') || optInputs(varargin, 'overwrite')
            copyfile(highres_2mm, [glmdir 'highres_2mm.nii.gz'],'f');
        end
        
        brainmask = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/brain_mask.nii.gz'];
        brainmask_2mm = [glmdir 'brain_mask_highres_2mm.nii.gz'];
        if ~exist(brainmask_2mm,'file') || optInputs(varargin, 'overwrite')
            unix_fsl(fsl_version, ['flirt -interp trilinear -in ' brainmask ' -ref ' highres_2mm ' -applyxfm -init ' exfunc2highres_regmat ' -out ' brainmask_2mm]);
            unix_fsl(fsl_version, ['fslmaths ' brainmask_2mm ' -thr 0.25 -bin ' brainmask_2mm]);
        end
        
        exfunc2highres_regmat = [regdir 'example_func2highres.mat']; % functional to anatomical registration matrix
        statistic_highres = strrep(statistics_to_register{i}, '.nii.gz','_highres_2mm.nii.gz');
        if ~exist(statistic_highres,'file') || optInputs(varargin, 'overwrite')
            unix_fsl(fsl_version, ['flirt -interp trilinear -in ' statistics_to_register{i} ' -ref ' highres_2mm ' -applyxfm -init ' exfunc2highres_regmat ' -out ' statistic_highres]);
            unix_fsl(fsl_version, ['fslmaths ' statistic_highres ' -mas ' brainmask_2mm ' ' statistic_highres]);
        end
    end
end

close all;
%% Old Code


% output files
% x = [];
% if tp_exclusion_cutoff < inf
%   x = [x, '_tpexclusion' num2str(tp_exclusion_cutoff)];
% end
%
% if cutoff_sec < inf
%   x = [x, '_hpfilt-cutoff' num2str(cutoff_sec)];
% end


%   plot(global_derivative_norm);

%   hist(global_derivative,100);

%   tp_exclusion_cutoff = 3;
%   figure
%   plot(sort(global_derivative_norm));
%   hold on;
%   plot([1,length(global_mean)],tp_exclusion_cutoff*[1 1],'r');
%   plot([1,length(global_mean)],0*[1 1],'k--');
%   plot([1,length(global_mean)],1*[1 1],'k');

%   if tp_exclusion_cutoff < inf
%       figure;
%       plot(abs(global_mean));
%       hold on;
%       plot([1,length(global_mean)],tp_exclusion_cutoff*std(global_mean)*[1 1],'r');
%       xlabel('TRs'); ylabel('abs(global mean)');
%       drawnow;
%   end

%% Prewhitening scraps


%     % tapered acf, effectively smoothing spectrum
%     max_lag = ceil((nTR*2+1)/2);    
%     M = round(100/TR); % parameter of tukey filter
%     tukey_window = (1/2)*(1+cos(pi*(0:M-1)/M))' * ones(1,n_nonzero_voxels);
%     acf_onesided = acf(1:max_lag,:);
%     acf_tukey_onesided = acf_onesided;
%     acf_tukey_onesided(M+1:end,:) = 0;
%     acf_tukey_onesided(1:M,:) = tukey_window .* acf_onesided(1:M,:);
%     if mod(size(acf,1),2) == 0
%         acf_tukey_twosided = [acf_tukey_onesided; acf_tukey_onesided(end-1:-1:2,:)];
%     else
%         acf_tukey_twosided = [acf_tukey_onesided; acf_tukey_onesided(end:-1:2,:)];
%     end
%     
%     % invert acf in frequency domain
%     pspec_orig = real(fft(acf));
%     pspec = real(fft(acf_tukey_twosided)); % real to deal with numeric error
%     pspec_inv = pspec;
%     pspec_inv(1,:) = 0;
%     pspec_inv(2:end,:) = 1./sqrt(pspec(2:end,:));
%     K = real(ifft(pspec_inv));
%     K = K ./ (ones(nTR*2,1)*K(1,:));
%     
%     % prewhiten data
%     FK = fft(K);
%     Ypad = [Yhp; zeros(size(Yhp))];
%     FYpad = fft(Ypad);
%     Y_white = real(ifft(FYpad .* FK));
%     
%     % estimate betas
%     tic
%     Xnuissance = [Xhp,zscore(nuissance_regressors)];
%     Xpad = [Xnuissance; zeros(size(Xnuissance))];
% %     FXpad = fft(Xpad);
% %     p = size(Xpad,2);
% %     B2 = zeros(p, n_nonzero_voxels);
% %     E_white = nan(nTR*2,n_nonzero_voxels);
% %     rvar = nan(1,n_nonzero_voxels);
% %     df = nTR-p;
% %     for i = 1:n_nonzero_voxels
% %         X = real(ifft(FXpad .* (FK(:,i)*ones(1,p)))); % apply whitening matrix to regression matrix
% %         B2(:,i) = pinv(X)*Y_white(:,i);
% %         E_white(:,i) = Y_white(:,i)-X*B2(:,i);
% %         rvar(i) = sum(E_white(:,i).^2)/df;
% %     end
% %     toc
% %     
% %     acf_white2 = real(ifft(abs(fft(E_white)).^2));
% %     spec_white2 = abs(fft(E_white)).^2;
% %     spec_raw = abs(fft(Epad)).^2;
% %     acf_raw = ifft(abs(fft(Epad)).^2);
% %     

% 
% 
%     
%     x = rand(3,1);
%     n = size(x,1);
%     xpad = [zeros(n,1); x; zeros(n,1)];
%     pspec = abs(fft(xpad)).^2;
%     acf = ifft(pspec);
% 
%     acf_white = 
%     
%     % residual ACF
%     zero_lag = size(E,1);
%     acf = nan(size(E,1)*2-1,size(E,2));
%     for i = 1:size(E,2);
%         acf(:,i) = xcorr(E(:,i),'unbiased');
%     end
%     
% 
%     % tukey windowing
%     tukey_window = (1/2)*(1+cos(pi*(0:M-1)/M))' * ones(1,n_nonzero_voxels);
%     acf_onesided = acf(zero_lag:end,:);
%     acf_tukey_onesided = acf_onesided;
%     acf_tukey_onesided(M+1:end,:) = 0;
%     acf_tukey_onesided(1:M,:) = tukey_window .* acf_onesided(1:M,:);
%     acf_tukey = [acf_tukey_onesided(end:-1:2,:); acf_tukey_onesided];
%     
%     % invert in spectral domain
%     power_spectrum = real(fft(ifftshift(acf_tukey))); % imaginary part is numerical error
%     inverse_spectrum = 1./sqrt(power_spectrum);
%     inverse_kernel = fftshift(ifft(inverse_spectrum));
%     
%     new_residual = conv(inverse_kernel,E,'same');
% 
%     
%     % acf_white = conv(acf_tukey(:,1), acf_inverse(:,1),'same'); check whitened acf is a delta function
%     % plot(acf_white)
%     K = convmtx(inverse_kernel)
%     
%     
%     Vi = toeplitz(ifftshift(acf_inverse(:,1)));
%     [X,S2,~] = svd(Vi);
%     [Y,~,~] = svd(Vi');
%     Ki = X*sqrt(S2)*Y;
%     S = Ki;
%     
%     % fit
%     Xnuissance = [Xhp,zscore(nuissance_regressors)];
%     B2 = zeros(size(Xnuissance,2), size(data,2));
%     B2(:,nonzero_voxels) = pinv(S*Xnuissance)*S*Yhp;
% 
%     % PSC values
%     psc = zeros(size(Xnuissance,2), size(data,2));
%     psc(:,nonzero_voxels) = B(:,nonzero_voxels)./ (ones(size(B,1),1)*voxel_means);
%     psc(:,any(isnan(psc))) = 0;
%     
%     
%     acf_mean = mean(acf,2);
%     acf_onesided = acf_mean(size(E,1):end);
%     Ecorr = toeplitz(acf_onesided);
%     
%     % plot and save autocorrelation function figure
%     figure;
%     plot((0:size(E,1)-1)*TR,acf_onesided,'k-o','LineWidth',2);  
%     xlabel('Time (s)')
%     ylabel('Average correlation (r)');
%     box off;
%     export_fig([glmdir 'acf.pdf'],'-pdf','-nocrop','-transparent');
%     