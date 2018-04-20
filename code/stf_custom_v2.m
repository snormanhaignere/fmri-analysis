function stf_custom_v2(exp,us,runtype,r,model,varargin)

% function stf_custom(exp,us,runtype,r,model,varargin)
%
% convolves boxcar with custom hrf to generate timing files for FSL

% sampling rate for convolution
sr = 5;
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
addpath(genpath('export_fig'));
fsl_version = read_fsl_version(exp, varargin{:});
freesurfer_version = read_freesurfer_version(exp, varargin{:});

% type of hrf function to use
% current options are 'BOLD', 'MION', 'MION_CUSTOM1'
hrf_name = read_scanparams(exp,us,runtype,'run',r,'hrf_type',varargin{:});
if optInputs(varargin, 'hrf_name')
    hrf_name = varargin{optInputs(varargin, 'hrf_name')+1};
end

% filter cutoff
cutoff_sec = inf;
if optInputs(varargin, 'cutoff_sec')
    cutoff_sec = varargin{optInputs(varargin, 'cutoff_sec')+1};
end

% whitematter PCs to regress out
n_whitematter_PCs = 0;
if optInputs(varargin, 'n_whitematter_PCs')
    n_whitematter_PCs = varargin{optInputs(varargin, 'n_whitematter_PCs')+1};
end

% stf directory
stfdir = [params('rootdir') exp '/analysis/stf/usub' num2str(us) '/'];
if ~exist(stfdir, 'dir');
    mkdir(stfdir);
end

%% Design matrix

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

% demean design matrix
X_design = model_scans;
X_design = X_design - ones(size(X_design,1),1)*mean(X_design);

% plot
if ~optInputs(varargin, 'noplot')
    figure;
    set(gcf, 'Position',[0 0 1440 700]/2);
    plot((0:nTR-1)'*TR, X_design);
    ylim([min(X_design(:))-0.1,max(X_design(:))+0.1]);
    ylabel('Response');
    xlabel('Time (s)');
    legend(conds,'Location','NorthEastOutside','Orientation','Vertical');
    export_fig([stfdir runtype '_r' num2str(r) '.pdf'],'-pdf','-nocrop');
end

%% nuissance regressors

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
    % read in volume file in order to measure white matter voxels
    volumefile = varargin{optInputs(varargin,'n_whitematter_PCs')+2};
    vol_fname = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/' volumefile '.nii.gz'];
    vol_br = MRIread(vol_fname);
    vol_dims = size(vol_br.vol);
    vol_data = reshape(vol_br.vol, [prod(vol_dims(1:3)), vol_dims(4)])';
    wm_matrix = vol_data(:, wm_mask & all(vol_data>max(vol_data(:))/1e6) & all(~isnan(vol_data)));
    
        
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
    export_fig([stfdir runtype '_r' num2str(r) '_whitematter_PC_explained_variance.pdf'],'-pdf','-transparent','-nocrop');

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
    export_fig([stfdir runtype '_r' num2str(r) '_whitematter_topPCs.pdf'],'-pdf','-transparent','-nocrop');
    
    %         sample_vox = squeeze(br.vol(44+1,35+1,58+1,:));
    %         sample_vox = sample_vox - mean(sample_vox);
    %         sample_vox_residual = sample_vox - U(:,1:n_whitematter_PCs)*pinv(U(:,1:n_whitematter_PCs))*sample_vox;
    X_nuissance = [X_nuissance, zscore(U(:,1:n_whitematter_PCs))];
    
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
    export_fig([stfdir runtype '_r' num2str(r) '_polynomial_regressors.pdf'],'-pdf','-transparent','-nocrop');
    
    X_nuissance = [X_nuissance, zscore(polynomial_regressors_unwrapped)];
end

%% optionally design matrix
if cutoff_sec ~= inf;
    
    cutoff_Hz = (1/cutoff_sec);
    nyq_Hz = (1/TR)/2;
    [B,A] = butter(4,cutoff_Hz/nyq_Hz,'high');
    
    % filter regressors
    Xhp = zeros(size(X_design));
    for i = 1:size(X_design,2)
        Xhp(:,i) = filtfilt(B,A,X_design(:,i));
    end
    
    X_nuissance_hp = zeros(size(X_nuissance));
    for i = 1:size(X_nuissance,2)
        X_nuissance_hp(:,i) = filtfilt(B,A,X_nuissance(:,i));
    end
    
    X_design = Xhp;
    X_nuissance = X_nuissance_hp;
    
end

%% Write design matrix to stf file using one-column fsl format

for i = 1:length(conds)        
    regfile = [stfdir runtype '_r'  num2str(r) '_' conds{i} '_custom.stf'];
    fid = fopen(regfile,'w');
    fprintf(fid,repmat('%.6f\n',1,nTR),X_design(:,i));
    fclose(fid);
end

%% Write nuissance regressors to file

% optionally add polynomial regressors to model
p = size(X_nuissance,2);
if p > 0
  regfile = [stfdir runtype '_r'  num2str(r) '_confound_regressors.stf'];
  fid = fopen(regfile,'w');
  for i = 1:nTR
    for j = 1:p
      fprintf(fid,'%.6f ',X_nuissance(i,j));
    end
    fprintf(fid,'\n');
  end
  fclose(fid);
end

close all;
