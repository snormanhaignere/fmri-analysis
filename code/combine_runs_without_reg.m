function combine_runs_without_reg(exp, us, runtype_original, r_original, runtype_combined, r_combined, fname, varargin)

% function combine_runs_without_reg(exp, us, runtype_original, r_original, fname, runtype_combined, r_combined, varargin)
% concatenates a set of runs all of the same runtype into a larger file
% with a new runtype and run number
% runs are concatenated without any initial registration
% intended to be used mainly with monkey data

% freesurfer matlab scrips
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));

% preprocessing directory for new combined runtype
combined_preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_combined '_r' num2str(r_combined) '/'];
if ~exist(combined_preprocdir,'dir');
  mkdir(combined_preprocdir)
end

% output file
outputfile = [combined_preprocdir fname '.nii.gz'];

if ~exist(outputfile, 'file') || optInputs(varargin, 'overwrite')
  cat_data = [];
  for j = 1:length(r_original)
    
    fprintf('Run %d\n',r_original(j)); drawnow;
    
    % read file
    %     if optInputs(varargin, 'monkey')
    %       run_specific_file = [params('rootdir') exp '/data/brain/nifti/usub' num2str(us) '/' runtype_original '_r' num2str(r_original(j)) '_falseheader.nii.gz'];
    %     else
    %       run_specific_file = [params('rootdir') exp '/data/brain/nifti/usub' num2str(us) '/' runtype_original '_r' num2str(r_original(j)) '_LAS.nii.gz'];
    %     end
    
    
    % input file
    run_specific_file = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_original '_r' num2str(r_original(j)) '/' fname '.nii.gz'];
    
    try
      br = MRIread(run_specific_file);
    catch
      keyboard;
    end
    
    % check number of TRs is correct
    [~, ~, ~, ~, ~, ~, ~, nTR, ~] = read_scanparams(exp,us,runtype_original,'run',r_original(j),varargin(:));
    if nTR ~= size(br.vol,4)
        fprintf('comine_runs_without_reg: Incorrect number of TRs in run %d\n',r_original(j)); drawnow;
        keyboard;
    end
    
    % add data to matrix
    if j == 1;
      cat_data = br.vol;
    else
      cat_data = cat(4, cat_data, br.vol);
    end
    
  end
  
  % write to file
  br.vol = cat_data;
  br.nframes = size(cat_data,4);
  br.fspec = outputfile;
  MRIwrite(br, outputfile);
  
end

%% create new para file for combined data

new_datafile = [params('rootdir') exp '/data/experiment/usub' num2str(us) '/seq/' runtype_combined '_r' num2str(r_combined) '_block.par'];
fid = fopen(new_datafile, 'w');
onset_index = 0;
for j = 1:length(r_original)
  
  % read timing information
  b = read_timings(exp,us,runtype_original,r_original(j),varargin{:});
  
  % write timing info to new para file
  for k = 1:length(b.onsets)
    if length(b.conds{k}) > 48
      error('Condition string too long');
    end
    fprintf(fid, '%8.2f%5d%8.2f%5d%50s\n', b.onsets(k) + onset_index, b.condition_indices(k), b.durs(k), 1, b.conds{k});
  end
  
  % update onset for next run
  [~, ~, TR, ~, ~, ~, ~, nTR, ~] = read_scanparams(exp,us,runtype_original,'run',r_original(j),varargin(:));
  onset_index = onset_index + nTR*TR;
end
fclose(fid);

%% Old Code

% end

%     clear br;
%     cat_data = [];
%     for j = 1:length(old_runtypes)
%         br = readmr([params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' old_runtypes{j} '_r' num2str(runs(r_original)) '/smooth' params('smooth') 'mm_intnorm_hpfilt' num2str(params('hpcutoff')) '.nii.gz'],'NOPROGRESSBAR');
%         cat_data = cat(4, cat_data, br.data);
%     end
%
%     br_new = br;
%     [blockdur, nulldur, TR, TA, stimdur, stim2scan, win, nTR, disdaqs] = read_scanparams(exp,us,new_runtype,varargin(:));
%     br_new.info.dimensions(4).size = nTR;
%     br_new.data = cat_data;