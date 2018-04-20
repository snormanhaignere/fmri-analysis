function hpfilt_surf(exp,us,runtype,r,hemi,varargin)

subjid = [exp '_us' num2str(us)];

filter_order = 2;
if optInputs(varargin, 'filter_order')
    filter_order = varargin{optInputs(varargin, 'filter_order')+1};
end

cutoff = inf;
if optInputs(varargin, 'cutoff')
    cutoff = varargin{optInputs(varargin, 'cutoff')+1};
end

demean_flag = '';
if optInputs(varargin,'demean')
    demean_flag = '_demean';
end
if optInputs(varargin, 'nulldemean')
    demean_flag = '_nulldemean';
end

[~,fwhm] = read_smooth(exp, varargin{:});
if optInputs(varargin, 'fwhm');
    fwhm = varargin{optInputs(varargin, 'fwhm')+1};
end

[~, ~, TR] = read_scanparams(exp, us, runtype, 'run', r, varargin{:});
sr = 1/TR;

if optInputs(varargin, 'monkey')
    preprocdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
else
    preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype '_r' num2str(r) '/'];
end
unfiltered_file = [preprocdir hemi '.brain_thresh_detrend1' demean_flag  '_smooth' num2str(100*fwhm, '%.0f') 'mm.mgz'];
filtered_file = [preprocdir hemi '.brain_thresh_detrend1' demean_flag  '_smooth' num2str(100*fwhm, '%.0f') 'mm_hpfilt-' num2str(cutoff) 'sec-order' num2str(filter_order) '.mgz'];

if ~exist(filtered_file,'file') || optInputs(varargin,'overwrite')
    if cutoff == inf
        copyfile(unfiltered_file, filtered_file, 'f');
    else
        f = MRIread(unfiltered_file);
        original_data = squeeze(f.vol)';
        freq_cutoff = (1/cutoff);
        [B,A] = butter(2,freq_cutoff/(sr/2),'high');
        filtered_data = filtfilt(B,A,original_data) + ones(size(original_data,1),1)*mean(original_data);
        %         plot([original_data(:,10000),filtered_data(:,10000)],'LineWidth',1);
        fn = f;
        fn.vol = nan(size(f.vol));
        fn.vol(1,:,1,:) = filtered_data';
        fn.fspec = filtered_file;
        MRIwrite(fn, filtered_file);
    end
end

