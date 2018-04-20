function estimate_hrf_surf(exp,usubs,runtype,input_fname,roi,hrf_type,win,hrf_name,varargin)

% usubs = [157 158];
% sr = 5;
% % input_fname = 'brain_thresh_detrend1_demean_smooth300mm';
% input_fname = 'motcorr_smooth571mm_detrend1';
% roi = 'stp';
% hrf_type = 'onset_response_function';
% % win = 0:(1/sr):40;
% win = 0:(1/sr):100;
% % exvar_threshold = 1;
% exvar_threshold = 0.1;
% % hrf_name = 'naturalsound_block_response_function';
% % hrf_name = 'naturalsound_monkey_block_response_function';
% hrf_name = 'pitch_localizer_monkey_block_response_function';
% varargin = {'monkey'};
% addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
% addpath(genpath('export_fig'));

usubs_amusia  = [45,49,51,53,55,57,59,71,73,75,171]; % 47
usubs_control = [48,56,46,52,54,172,50,58,74,72,70]; % 44
usubs = [usubs_amusia, usubs_control];
exp = 'amusia';
sr = 5;
input_fname = 'motcorr_smooth300mm';
roi = 'stp';
hrf_type = 'onset_response_function';
% win = 0:(1/sr):40;
win = 0:(1/sr):100;
% exvar_threshold = 1;
exvar_threshold = 0.1;
runtype = 'localizer';
% hrf_name = 'naturalsound_block_response_function';
% hrf_name = 'naturalsound_monkey_block_response_function';
hrf_name = 'amusia_irf';
varargin = {};
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
addpath(genpath('export_fig'));

%%
hrf_all_subjects = nan(length(win), length(usubs));
for j = 1:length(usubs)
  
  runs = read_runs(exp,usubs(j),runtype,varargin{:});
  hrf_all_runs = nan(length(win), length(runs));
  for k = 1:length(runs)
    
    %% design matrix
    fprintf('us %d\n',usubs(j));
    conds = read_conditions(exp,usubs(j),runtype,'run',runs(k),varargin{:});
    b = read_timings(exp,usubs(j),runtype,runs(k),varargin{:});
    % hrf_type = read_scanparams(exp,usubs(j),runtype,'hrf_type',varargin{:});
    TR = read_scanparams(exp,usubs(j),runtype,'TR',varargin{:},'run',runs(k));
    nTR = read_scanparams(exp,usubs(j),runtype,'nTR',varargin{:},'run',runs(k));
    
    % upsampled block indicator variable
    totaltime = (nTR-1)*TR;
    nsmps = floor(totaltime*sr);
    seq_upsampled = zeros(nsmps,1);
    t = (0:nsmps-1)/sr;
    for i = 1:length(b.onsets)
      if ~strcmp(b.conds{i}, 'NULL')
        xi = t >= b.onsets(i) & t < b.onsets(i)+b.durs(i);
        seq_upsampled(xi) = find(strcmp(b.conds{i}, conds));
      end
    end
    
    % candlestick/boxcar regressors
    switch hrf_type
      case 'onset_response_function'
        seq_binary = zeros(size(seq_upsampled));
        seq_binary(1) = seq_upsampled(1)~=0;
        seq_binary(2:end) = seq_upsampled(2:end)~=0 & seq_upsampled(2:end)~=seq_upsampled(1:end-1);
      case 'impulse_response_function'
        seq_binary = seq_upsampled~=0;
      otherwise
        error('hrf_type does not match options.')
    end
    
    % fir model regressors
    fir_model = zeros(length(seq_binary),length(win));
    for i = 1:length(win)
      fir_model(i:end,i) = seq_binary(1:end-(i-1));
    end
    
    %% regression
    
    subjid = [exp '_us' num2str(usubs(j))];
    if optInputs(varargin, 'monkey')
      surfdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(usubs(j)) '/' runtype '_r' num2str(runs(k)) '/'];
    else
      surfdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(runs(k)) '/'];
    end
    
    % read data
    hemis = {'rh','lh'};
    roi_data = cell(1,2);
    for q = 1:2
      
      % read in roi label
      if optInputs(varargin, 'monkey')
        roifile = [params('rootdir') 'freesurfer/' subjid '/label/' hemis{q} '.' roi '.label'];
      else
        roifile = [params('rootdir') 'freesurfer/fsaverage/label/' hemis{q} '.' roi '.label'];
      end
      lb = read_label(roifile);
      
      inputfile = [surfdir hemis{q} '.' input_fname '.mgz'];
      br = MRIread(inputfile);
      roi_data{q} = squeeze(br.vol(:,lb.vnums+1,:,:))';
      
    end
    
    % concatenate two hemispheres and upsample
    roi_data_both_hemi = cat(2,roi_data{:});
    roi_data_upsamp = interp1((0:nTR-1)'*TR, roi_data_both_hemi, t, 'cubic');
    
    % demean timepoints
    voxel_means = mean(roi_data_upsamp);
    Y = roi_data_upsamp - ones(size(roi_data_upsamp,1),1)*voxel_means;
    
    % deman design matrix
    X = fir_model - ones(size(fir_model,1),1)*mean(fir_model);
    
    B = pinv(X)*Y;
    
    resid = Y-X*B;
    exvar = (var(Y) - var(resid))./var(Y);
    [~,xi] = sort(exvar,'descend');
    xi = xi(1:round(exvar_threshold*length(exvar)));
    hrf = mean(B(:,xi),2);
    figure(1);
    plot(win,hrf);
    title(sprintf('us %d, run %d\n',usubs(j),runs(k)));
    drawnow;
    hrf_all_runs(:,k) = hrf;
  end
  
  hrf_all_subjects(:,j) = mean(hrf_all_runs,2);
  
end

%% Save HRF data
hrf_all_subjects_norm = hrf_all_subjects ./ (ones(size(hrf_all_subjects,1),1)*max(abs(hrf_all_subjects)));
hrf_mean = mean(hrf_all_subjects_norm,2);
hrf_sem = stderr_withsub_corrected(hrf_all_subjects_norm')';

figure;
% errorbar(win(1:sr:end),hrf_mean(1:sr:end),hrf_sem(1:sr:end));
plot(win,hrf_mean);
xlabel('Time (s)'); ylabel('Response (au)');
box off;
export_fig([params('rootdir') 'custom_hrfs/' hrf_name '.pdf'],'-nocrop','-transparent','-pdf');

figure;
colormap(distinguishable_colors(10,1));
plot(win,hrf_all_subjects_norm);
xlabel('Time (s)'); ylabel('Response (au)');
box off;
export_fig([params('rootdir') 'custom_hrfs/' hrf_name '_allsubjects.pdf'],'-nocrop','-transparent','-pdf');

hrf = hrf_mean/max(abs(hrf_mean));
hrf(1) = 0;
save([params('rootdir') 'custom_hrfs/' hrf_name '.mat'],'hrf_all_subjects','hrf','win');

[B,A] = butter(4,0.1/(sr/2),'low');
hrf_lowpass = filtfilt(B,A,hrf);
hrf_lowpass = hrf_lowpass/max(abs(hrf_lowpass));
hrf_lowpass(1) = 0;
plot(win,[hrf,hrf_lowpass]);
xlabel('Time (s)'); ylabel('Response (au)');
box off;
export_fig([params('rootdir') 'custom_hrfs/' hrf_name '_lowpass.pdf'],'-nocrop','-transparent','-pdf');

hrf = hrf_lowpass;
save([params('rootdir') 'custom_hrfs/' hrf_name '_lowpass.mat'],'hrf_all_subjects','hrf','win');

% hrf_error = std(B,[],2) / sqrt(size(B,2)-1);
% errorbar(win,hrf);
  
%% Old Code
%   if win(i) < 0
%     fir_model(1:end-abs(win(i)),i) = seq_binary( 1+abs(win(i)):end );
%   elseif win(i) > 0
%     fir_model(1+abs(win(i)):end,i) = seq_binary( 1:end-abs(win(i)) );
%   else
%     fir_model(:,i) = seq_binary;
%   end

% model_scans = model_scans/max(model_scans(:));
