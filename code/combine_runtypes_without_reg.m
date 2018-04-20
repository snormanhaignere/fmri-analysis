function combine_runtypes_without_reg(exp, us, runtypes_original, runtype_combined, runs, fname, varargin)

% Combines corresponding runs from different runtypes. Primarily useful when
% there are too many conditions to be included in a single run, and the
% conditions are thus split-up between a few different runtypes.
%
% 2016-020-15: Created, Sam NH

% freesurfer matlab scrips
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));

for j = 1:length(runs)
    
    % preprocessing directory for new combined runtype
    combined_preprocdir = [params('rootdir') exp '/analysis/preprocess' ...
        '/usub' num2str(us) '/' runtype_combined '_r' num2str(runs(j)) '/'];
    if ~exist(combined_preprocdir,'dir');
        mkdir(combined_preprocdir)
    end
    
    % output file
    outputfile = [combined_preprocdir fname '.nii.gz'];
    
    if ~exist(outputfile, 'file') || optInputs(varargin, 'overwrite')
        
        fprintf('Run %d\n',runs(j)); drawnow;
        
        cat_data = [];
        for i = 1:length(runtypes_original)
            
            % input file
            run_specific_file = [...
                params('rootdir') exp '/analysis/preprocess/usub' num2str(us) ...
                '/' runtypes_original{i} '_r' num2str(runs(j)) '/'...
                fname '.nii.gz'];
            
            try
                br = MRIread(run_specific_file);
            catch me
                print_error_message(me);
                keyboard;
            end
            
            % check number of TRs is correct
            [~, ~, ~, ~, ~, ~, ~, nTR, ~] = read_scanparams(...
                exp,us,runtypes_original{i},'run',runs(j),varargin{:});
            if nTR ~= size(br.vol,4)
                fprintf(['comine_runs_without_reg: ' ...
                    'Incorrect number of TRs in run %d\n'],runs(j));
                drawnow;
                keyboard;
            end
            
            % add data to matrix
            if i == 1;
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
    
end

%% create new para file for combined data

for j = 1:length(runs)
    
    new_datafile = [params('rootdir') exp '/data/experiment/usub' num2str(us) ...
        '/seq/' runtype_combined '_r' num2str(runs(j)) '_block.par'];
    fid = fopen(new_datafile, 'w');
    onset_index = 0;
    
    for i = 1:length(runtypes_original)
                
        % read timing information
        b = read_timings(exp,us,runtypes_original{i},runs(j),varargin{:});
        
        % write timing info to new para file
        for k = 1:length(b.onsets)
            if length(b.conds{k}) > 48
                error('Condition string too long');
            end
            fprintf(fid, '%8.2f%5d%8.2f%5d%50s\n', b.onsets(k) + onset_index, ...
                b.condition_indices(k), b.durs(k), 1, b.conds{k});
        end
        
        % update onset for next run
        [~, ~, TR, ~, ~, ~, ~, nTR, ~] = read_scanparams(...
            exp,us,runtypes_original{i},'run',runs(j),varargin(:));
        onset_index = onset_index + nTR*TR;
    end
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