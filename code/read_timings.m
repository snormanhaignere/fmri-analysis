function b = read_timings(exp,us,runtype,r,varargin)

% b = read_timings(s,r,varargin)
% 
% reads in timing information from the experimental output file
% inputs are the subject number, runtype, and run number

%% read in data
datadir = [params('rootdir') exp '/data/experiment/'];

if strcmp(exp,'naturalsound') && strcmp(runtype, 'texture_combined')
    load('texture_experiment_stimulus_order.mat');
    uncombined_runs = (1:12) + (r-1)*12;
    b.onsets = nan(18*12, 1);
    b.condition_indices = nan(18*12,1);
    b.conds = cell(18*12,1);
    b.durs = 17*ones(18*12,1);
    [~, ~, TR, ~, ~, ~, ~, nTR, ~] = read_scanparams('naturalsound',us,'texture',varargin(:));
    for i = 1:length(uncombined_runs)
        event_file = [params('rootdir') 'naturalsound/data/experiment/usub' num2str(us) '/texture_r' num2str(uncombined_runs(i)) '_event.txt'];
        block_file = [params('rootdir') 'naturalsound/data/experiment/usub' num2str(us) '/texture_r' num2str(uncombined_runs(i)) '_block.txt'];
        [rb.onsets, rb.condition_indices, rb.synthesis_type] = textread(block_file, '%f%d%s', 'headerlines', 1); %#ok<*NASGU>
        [re.onsets, re.condition_indices, re.synthesis_type, ~, ~, re.stim, ~, ~, ~, ~] = textread(event_file, '%f%d%s%d%s%f%f%f%f%f', 'headerlines', 1); %#ok<*NASGU>
        rb.conds = cell(length(rb.onsets),1);
        for j = 1:length(rb.onsets)
            if rb.condition_indices(j) == 0
                rb.conds{j} = 'NULL';
            else
                xi = rb.onsets(j) == re.onsets;
                if sum(xi) ~= 1
                    error('Error in read_timings.m: Should only be one onset match between block and event file.');
                end
                rb.conds{j} = [rb.synthesis_type{j} '-' strrep(texture_stimulus_order{re.stim(xi)}, '.wav', '')];
            end
        end
        b.conds((1:18) + (i-1)*18) = rb.conds;
        b.onsets((1:18) + (i-1)*18) = rb.onsets + (i-1)*(TR*nTR);
        b.condition_indices((1:18) + (i-1)*18) = rb.condition_indices;
    end
elseif strcmp(exp,'naturalsound') && strcmp(runtype, 'spectrotemporal_combined')
    load([params('rootdir') 'naturalsound/scripts_exp/spectrotemporal_stims.mat']); % loads cell array with stimuli in the same order as the experiment
    uncombined_runs = (1:10) + (r-1)*10;
    b.onsets = nan(21*10, 1);
    b.condition_indices = nan(21*10,1);
    b.conds = cell(21*10,1);
    b.durs = 17*ones(21*10,1);
    [~, ~, TR, ~, ~, ~, ~, nTR, ~] = read_scanparams('naturalsound',us,'spectrotemporal',varargin(:));
    for i = 1:length(uncombined_runs)
        event_file = [params('rootdir') 'naturalsound/data/experiment/usub' num2str(us) '/spectrotemporal_r' num2str(uncombined_runs(i)) '_event.txt'];
        block_file = [params('rootdir') 'naturalsound/data/experiment/usub' num2str(us) '/spectrotemporal_r' num2str(uncombined_runs(i)) '_block.txt'];
        [rb.onsets, rb.condition_indices, rb.synthesis_type] = textread(block_file, '%f%d%s', 'headerlines', 1); %#ok<*NASGU>
        [re.onsets, re.condition_indices, re.synthesis_type, ~, ~, re.stim, ~, ~, ~, ~] = textread(event_file, '%f%d%s%d%s%f%f%f%f%f', 'headerlines', 1); %#ok<*NASGU>
        rb.conds = cell(length(rb.onsets),1);
        if length(rb.onsets) ~= 21;
            error('Error in read_timings: wrong number of blocks');
        end
        for j = 1:length(rb.onsets)
            if rb.condition_indices(j) == 0
                rb.conds{j} = 'NULL';
            else
                xi = rb.onsets(j) == re.onsets;
                if sum(xi) ~= 1
                    error('Error in read_timings.m: Should only be one onset match between block and event file.');
                end
                if ~strcmp(rb.synthesis_type{j}, re.synthesis_type{xi})
                    error('Error in read_timings.m: Block condition names in event and block file do not match.');
                end
                rb.conds{j} = [rb.synthesis_type{j} '-si' num2str(re.stim(xi)) '-' strrep(stims{re.stim(xi)}, '.wav', '')];
            end
        end
        b.conds((1:21) + (i-1)*21) = rb.conds;
        b.onsets((1:21) + (i-1)*21) = rb.onsets + (i-1)*(TR*nTR);
        b.condition_indices((1:21) + (i-1)*21) = rb.condition_indices;
    end
elseif strcmp(exp,'naturalsound') && strcmp(runtype, 'spectrotemporal_v2_combined')
    load([params('rootdir') 'naturalsound/scripts_exp/spectrotemporal_stims_v4.mat']); % loads cell array with stimuli in the same order as the experiment
    nr = 12;
    nb = 18;
    uncombined_runs = (1:nr) + (r-1)*nr;
    b.onsets = nan(nb*nr, 1);
    b.condition_indices = nan(nb*nr,1);
    b.conds = cell(nb*nr,1);
    b.durs = 17*ones(nb*nr,1);
    [~, ~, TR, ~, ~, ~, ~, nTR, ~] = read_scanparams('naturalsound',us,'spectrotemporal_v2',varargin(:));
    for i = 1:length(uncombined_runs)
        event_file = [params('rootdir') 'naturalsound/data/experiment/usub' num2str(us) '/spectrotemporal_v2_r' num2str(uncombined_runs(i)) '_event.txt'];
        block_file = [params('rootdir') 'naturalsound/data/experiment/usub' num2str(us) '/spectrotemporal_v2_r' num2str(uncombined_runs(i)) '_block.txt'];
        [rb.onsets, rb.condition_indices, rb.synthesis_type] = textread(block_file, '%f%d%s', 'headerlines', 1); %#ok<*NASGU>
        [re.onsets, re.condition_indices, re.synthesis_type, ~, ~, re.stim,~, ~, ~, ~, ~] = textread(event_file, '%f%d%s%d%s%f%s%f%f%f%f', 'headerlines', 1); %#ok<*NASGU>
        rb.conds = cell(length(rb.onsets),1);
        if length(rb.onsets) ~= nb;
            error('Error in read_timings: wrong number of blocks');
        end
        for j = 1:nb
            if rb.condition_indices(j) == 0
                rb.conds{j} = 'NULL';
            else
                xi = rb.onsets(j) == re.onsets;
                if sum(xi) ~= 1
                    error('Error in read_timings.m: Should only be one onset match between block and event file.');
                end
                if ~strcmp(rb.synthesis_type{j}, re.synthesis_type{xi})
                    error('Error in read_timings.m: Block condition names in event and block file do not match.');
                end
                rb.conds{j} = [rb.synthesis_type{j} '-si' num2str(re.stim(xi)) '-' strrep(stims{re.stim(xi)}, '.wav', '')];
            end
        end
        b.conds((1:nb) + (i-1)*nb) = rb.conds;
        b.onsets((1:nb) + (i-1)*nb) = rb.onsets + (i-1)*(TR*nTR);
        b.condition_indices((1:nb) + (i-1)*nb) = rb.condition_indices;
    end
    
elseif strcmp(exp,'naturalsound-infant') && strcmp(runtype, 'main_v1')
    
    load([params('rootdir') exp '/data/experiment/usub' num2str(us) '/' 'r' num2str(r) '.mat'],'seq','measured_onsets','condnames','starting_block','last_completed_block'); % loads cell array with stimuli in the same order as the experiment
    b.onsets = measured_onsets(starting_block:last_completed_block);
    b.condition_indices = seq(starting_block:last_completed_block);
    b.durs = 12*ones(size(b.onsets));
    b.conds = cell(size(b.condition_indices));
    for j = 1:length(b.condition_indices);
        if b.condition_indices(j)==0
            b.conds{j} = 'NULL';
        else
            b.conds{j} = condnames{b.condition_indices(j)};
        end
    end
    
else
    try
        x1 = [datadir 'usub' num2str(us) '/seq/' runtype '_r' num2str(r) '_block.par'];
        x2 = [datadir 'usub' num2str(us) '/seq/' runtype '_standard_r' num2str(r) '_block.par'];
        x3 = [datadir 'usub' num2str(us) '/seq/' runtype '_r' num2str(r) '.par'];
        x4 = [datadir 'usub' num2str(us) '/seq/' runtype '__r' num2str(r) '_block.par'];
        if exist(x1,'file')
            [b.onsets, b.condition_indices, b.durs, ~, b.conds] = textread(x1, '%f%d%f%f%s'); %#ok<*NASGU>
        elseif exist(x2,'file')
            [b.onsets, b.condition_indices, b.durs, ~, b.conds] = textread(x2, '%f%d%f%f%s'); %#ok<*NASGU>
        elseif exist(x3,'file')
            [b.onsets, b.condition_indices, b.durs, ~, b.conds] = textread(x3, '%f%d%f%f%s'); %#ok<*NASGU>
        else
            [b.onsets, b.condition_indices, b.durs, ~, b.conds] = textread(x4, '%f%d%f%f%s'); %#ok<*NASGU>
        end
    catch me
        print_error_message(me);
        keyboard
    end
end

%% format

% use inds to remove null trials if specifie
if optInputs(varargin,'nofix');
    nullevents = strcmp('NULL',b.conds);
    b.onsets = b.onsets(~nullevents);
    b.conds = b.conds(~nullevents);
    b.durs = b.durs(~nullevents);
    b.condition_indices = b.condition_indices(~nullevents);
end

% b.stims = data{strmatch('stimname',fields,'exact')}(inds);
% b.stims = strrep(b.stims,'.wav','');
% 
% b.condstims = cell(size(b.stims));
% for i = 1:length(b.stims)
%     b.condstims{i} = [b.conds{i} '_s_' b.stims{i}]; 
% end
