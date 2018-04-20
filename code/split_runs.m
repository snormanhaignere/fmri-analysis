function split_runs(exp, us, runtype_combined, r_combined, runtype_split, r_split, fname, volume_or_surface, varargin)

% function split_runs(exp, us, runtype_combined, r_combined, runtype_split, r_split, fname, nTR_per_run, varargin)
% splits larger run into several smaller runs
% 
% 2017-03-06: Fixed bug that caused onset index not to be updated

addpath(genpath('/software/Freesurfer/5.3.0/matlab'));

% combine functional data
onset_index = 0;
for j = 1:length(r_split)
        
    switch volume_or_surface
        case 'volume'
            split_preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_split '_r' num2str(r_split(j)) '/'];
            if ~exist(split_preprocdir,'dir')
                mkdir(split_preprocdir);
            end
            outputfile = [split_preprocdir fname '.nii.gz'];
        case 'surface'
            subjid = [exp '_us' num2str(us)];
            if optInputs(varargin, 'monkey')
                split_preprocdir = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype_split '_r' num2str(r_split(j)) '/'];
            else
                split_preprocdir = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype_split '_r' num2str(r_split(j)) '/'];
            end
            if ~exist(split_preprocdir,'dir')
                mkdir(split_preprocdir);
            end
            if ~optInputs(varargin, 'hemi')
                error('Error in fit_glm: Need to specify hemisphere as optional argument for surface data.m');
            end
            hemi = varargin{optInputs(varargin, 'hemi')+1};
            outputfile = [split_preprocdir hemi '.' fname '.mgz'];
    end
    
    [~, ~, ~, ~, ~, ~, ~, nTR, ~] = read_scanparams(exp,us,runtype_split,'run',r_split(j),varargin(:));
    onset_index
    
    if ~exist(outputfile, 'file') || optInputs(varargin, 'overwrite')
        
        if ~exist('br','var')
            switch volume_or_surface
                case 'volume'
                    inputfile = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_combined '_r' num2str(r_combined) '/' fname '.nii.gz'];
                case 'surface'
                    subjid = [exp '_us' num2str(us)];
                    if optInputs(varargin, 'monkey')
                        inputfile = [params('rootdir') 'freesurfer/' subjid '/preprocess/usub' num2str(us) '/' runtype_combined '_r' num2str(r_combined) '/' hemi '.' fname '.mgz'];
                    else
                        inputfile = [params('rootdir') 'freesurfer/fsaverage/preprocess/' subjid '/' runtype_combined '_r' num2str(r_combined) '/' hemi '.' fname '.mgz'];
                    end
            end
            br = MRIread(inputfile);
        end
        fprintf('Run %d\n',r_split(j)); drawnow;
        
        drawnow;
        splitbr = br;
        splitbr.vol = br.vol(:,:,:,(1:nTR)+onset_index);
        splitbr.nframes = nTR;
        splitbr.fspec = outputfile;
        MRIwrite(splitbr, outputfile);
        
    end
    
    % update onset counter
    onset_index = onset_index + nTR;
end

%% Create new split para files

% read timing information
b = read_timings(exp,us,runtype_combined,r_combined,varargin{:});

onset_index = 0;
for j = 1:length(r_split)
    
    new_datafile = [params('rootdir') exp '/data/experiment/usub' num2str(us) '/seq/' runtype_split '_r' num2str(r_split(j)) '_block.par'];
    fid = fopen(new_datafile, 'w');
    
    % timing information
    [~, ~, TR, ~, ~, ~, ~, nTR, ~] = read_scanparams(exp,us,runtype_split,'run',r_split(j),varargin(:));
    
    % write timing info to new para file
    for k = find(b.onsets >= onset_index & b.onsets < onset_index + TR*nTR)'
        if length(b.conds{k}) > 48
            error('Condition string too long');
        end
        fprintf(fid, '%8.2f%5d%8.2f%5d%50s\n', b.onsets(k) - onset_index, b.condition_indices(k), b.durs(k), 1, b.conds{k});
    end
    
    onset_index = onset_index + nTR*TR;
    fclose(fid);

end

%% Copy registration files
% registration directories to copy over
reg_directories = {'reg_flirt','reg_handtune','reg_bbreg'};
for q = 1:length(reg_directories)
    
    combined_flirt_reg_directory = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_combined '_r' num2str(r_combined) '/' reg_directories{q} '/'];
    if exist(combined_flirt_reg_directory,'dir')
        files = mydir(combined_flirt_reg_directory);
        for j = 1:length(r_split)
            split_flirt_reg_directory = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_split '_r' num2str(r_split(j)) '/' reg_directories{q} '/'];
            if ~exist(split_flirt_reg_directory,'dir')
                mkdir(split_flirt_reg_directory);
            end
            for i = 1:length(files)
                if ~exist([split_flirt_reg_directory files{i}],'file') || optInputs(varargin, 'overwrite')
                    copyfile([combined_flirt_reg_directory files{i}],[split_flirt_reg_directory files{i}],'f')
                end
            end
        end
    end
end

%% Copy over other useful files
misc_files = {'example_func.nii.gz','example_func_bet.nii.gz','brain_mask.nii.gz'};
for j = 1:length(r_split)
    for i = 1:length(misc_files)
        combined_preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_combined '_r' num2str(r_combined) '/'];
        split_preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype_split '_r' num2str(r_split(j)) '/'];
        if ~exist([split_preprocdir misc_files{i}],'file') || optInputs(varargin, 'overwrite')
            copyfile([combined_preprocdir misc_files{i}],[split_preprocdir misc_files{i}],'f')
        end
    end
end

%% Runorders file for new split runtype
runorders_file = [params('rootdir') exp '/data/brain/runorders/' runtype_split '_us' num2str(us) '.txt'];
if exist(runorders_file,'file')
    [runs,~,~] = textread(runorders_file, '%d%d%d'); %#ok<*NASGU>
else
    runs = [];
end

% add new runs
if any(~ismember(r_split,runs))
    fid = fopen(runorders_file,'a');
    try
        [~,dicom,scan] = read_runs(exp, us, runtype_combined, varargin{:}, 'runnum', r_combined);
    catch
        fprintf('Problem reading run information in split_runs.m\n'); drawnow;
        keyboard;
    end
    for j = 1:length(r_split)
        if ~ismember(r_split(j),runs)
            fprintf(fid, '%d %d %d\n', r_split(j), dicom, scan);
        end
    end
    fclose(fid);
end

% sort runs
[runs,dicoms,scans] = textread(runorders_file, '%d%d%d'); %#ok<*NASGU>
[~,xi] = sort(runs,'ascend');
if ~isequal(xi,(1:length(runs))')
    fid = fopen(runorders_file,'w');
    for j = 1:length(runs)
        fprintf(fid, '%d %d %d\n', runs(xi(j)), dicoms(xi(j)), scans(xi(j)));
    end
    fclose(fid);
end

