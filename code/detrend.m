function detrend(exp,us,runtype,r,varargin)

% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

polyorder = 1; % linear detrend default
fwhm = read_smooth(exp, varargin{:});
smooth_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm.nii.gz'];
intnorm_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm_intnorm.nii.gz'];
detrend_file = [preprocdir 'smooth' num2str(fwhm*100,'%.0f') 'mm_intnorm_detrend' num2str(polyorder) '.nii.gz'];

% brightness threshold
mode_file = [preprocdir 'mode.txt'];
fid = fopen(mode_file,'r');
x = textscan(fid,'%f'); fclose(fid);
factor = 1e4/x{1};

fsl_version = read_fsl_version(exp);

if ~exist([detrend_file],'file') || optInputs(varargin, 'overwrite')

    % sigma
    [~, ~, TR] = read_scanparams(exp,us,runtype,'run',r,varargin{:});
    sigma = params('hpcutoff')/(2*TR);
    
    % normalize mode
    unix_fsl(fsl_version, ['fslmaths ' smooth_file ' -mul ' num2str(factor) ' '  intnorm_file]);
        
    
    %% Detrend data
    fprintf('Detrending smoothed data: run %d\n', r);
    drawnow;
    
    % data matrix
    % unsmoothedbr = readmr([intnorm_file '.nii.gz'],'NOPROGRESSBAR');
    smoothedbr = MRIread([intnorm_file]);
    dims = size(smoothedbr.vol);
    data_cols = reshape(shiftdim(smoothedbr.vol,3), dims(4), prod(dims(1:3)));
    nonzero_cols = ~all(data_cols < 1e-3);
    data_nonzero_cols = data_cols(:,nonzero_cols);
    
    % polynomial fit
    dims = size(smoothedbr.vol);
    t = (1:dims(4))' - floor(dims(4)/2);
    t_matrix = (t*ones(1,polyorder+1)) .^ (ones(dims(4),1)*[0:polyorder]);
    predicted = t_matrix * pinv(t_matrix) * data_nonzero_cols;
    
    % residual with mean added back in
    resid = zeros(size(data_cols));
    resid(:,nonzero_cols) = data_nonzero_cols - predicted + ones(dims(4),1)*mean(data_nonzero_cols);
    
    % write new data file
    detrendbr = smoothedbr;
    detrendbr.vol = shiftdim(reshape(resid, [dims(4), dims(1:3)]),1);
    % bxhwrite(detrendbr, detrend_file);
    detrendbr.fspec = detrend_file;
    MRIwrite(detrendbr, detrend_file);

    %     %% Detrend data
    %     inormbr = readmr([intnorm_file '.nii.gz'],'NOPROGRESSBAR');
    %     detrendbr = inormbr;
    %     detrendbr.data = zeros(size(inormbr.data));
    %     dims = size(inormbr.data);
    %     t = (1:dims(4))' - floor(dims(4)/2);
    %     fprintf('Detrending slice...');
    %     drawnow;
    %     for i = 1:dims(1)
    %         fprintf('%d ',i);
    %         drawnow;
    %         for j = 1:dims(2)
    %             for k = 1:dims(3)
    %                 data = squeeze(inormbr.data(i,j,k,:));
    %                 if ~all(data < 1e-3)
    %                     [p,S] = polyfit(t,data,polyorder);
    %                     predicted = 0;
    %                     for q = 0:polyorder
    %                         predicted = predicted + p(polyorder+1-q) * t.^q;
    %                     end
    %                     detrendbr.data(i,j,k,:) = p(polyorder+1) + (data - predicted);
    %
    %                     %                     plot(t,[squeeze(inormbr.data(i,j,k,:)),squeeze(detrendbr.data(i,j,k,:))])
    %                     %                     drawnow;
    %                     %                     pause(0.5);
    %                 end
    %             end
    %         end
    %     end
    %     fprintf('\n');
    %     bxhwrite(detrendbr, detrend_file);
    
end

unsmoothed_file = [preprocdir 'brain_thresh.nii.gz'];
detrend_unsmoothed_file = [preprocdir 'brain_thresh_detrend' num2str(polyorder) '.nii.gz'];

if ~exist([detrend_unsmoothed_file],'file') || optInputs(varargin, 'overwrite')

 
    %% Detrend data
    fprintf('Detrending unsmoothed data: run %d\n', r);
    drawnow;
    
    % data matrix
    % unsmoothedbr = readmr([unsmoothed_file '.nii.gz'],'NOPROGRESSBAR');
    unsmoothedbr = MRIread([unsmoothed_file]);
    dims = size(unsmoothedbr.vol);
    data_cols = reshape(shiftdim(unsmoothedbr.vol,3), dims(4), prod(dims(1:3)));
    nonzero_cols = ~all(data_cols < 1e-3);
    data_nonzero_cols = data_cols(:,nonzero_cols);
    
    % polynomial fit
    dims = size(unsmoothedbr.vol);
    t = (1:dims(4))' - floor(dims(4)/2);
    t_matrix = (t*ones(1,polyorder+1)) .^ (ones(dims(4),1)*[0:polyorder]);
    predicted = t_matrix * pinv(t_matrix) * data_nonzero_cols;
    
    % residual with mean added back in
    resid = zeros(size(data_cols));
    resid(:,nonzero_cols) = data_nonzero_cols - predicted + ones(dims(4),1)*mean(data_nonzero_cols);
    
    % write new data file
    detrendbr = unsmoothedbr;
    detrendbr.vol = shiftdim(reshape(resid, [dims(4), dims(1:3)]),1);
    detrendbr.fspec = detrend_unsmoothed_file;
    % bxhwrite(detrendbr, detrend_unsmoothed_file);
    MRIwrite(detrendbr, detrend_unsmoothed_file);

    
end
