function tla_matlab(P,volume_or_surface,tla_directory_name,statistical_test,varargin)

% Example Data (from Amusia experiment)
clear P;
usubs_amusia = [45,49,51,53,55,57,59,71,73,75,171];
for i = 1:11;
    P(i).exp = 'amusia';
    P(i).us = usubs_amusia(i);
    P(i).runtype = 'localizer';
    P(i).runs = 1;
    P(i).contrast = 'harm_vs_noise';
    P(i).lower_level_directory_name = 'smooth300mm_grid_hand-stp-stg_1.5mm_10whitematterPCs';
end
tla_directory_name = 'amusia_group_permutation_harm_vs_noise';
volume_or_surface = 'downsampled_surface';
statistical_test = 'random';
varargin = {};

% scripts directories
source_directory = strrep(which('fla_matlab.m'),'fla_matlab.m','');
addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
addpath([source_directory 'export_fig']);

% output files
switch volume_or_surface
    case 'volume'
        
        error('Need to setup volume analysis');
        
    case 'surface'
        
        error('Need to setup surface analysis');

    case 'downsampled_surface'
        
        tla_directory = [params('rootdir') 'freesurfer/fsaverage/tla_matlab/' tla_directory_name '_downsampled_hash' DataHash(P) '/'];
        if ~exist(tla_directory,'dir');
            mkdir(tla_directory);
        end
        
        t_file =  [tla_directory  'tstat_' statistical_test '.mat'];
        p_file =  [tla_directory  'pstat_' statistical_test '.mat'];
        tla_cope_file =  [tla_directory  'cope_' statistical_test '.mat'];
        tla_cope_var_file =  [tla_directory  'cope_var_' statistical_test '.mat'];
        
    otherwise
        error('Error in sla_matlab: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"...');
        
end

if ~exist(t_file,'file') || ~exist(p_file,'file') || ~exist(tla_cope_file,'file') || ~exist(tla_cope_var_file,'file') || optInputs(varargin, 'overwrite')
    
    % read in data
    n_subjects = length(P);
    for i = 1:n_subjects
        
        % input files
        switch volume_or_surface
            case 'volume'

                error('Need to setup volume analysis');

            case 'surface'
                
                error('Need to setup surface analysis');

            case 'downsampled_surface'
                
                % subject id
                subjid = [P(i).exp '_us' num2str(P(i).us)];
                
                % first or second level analysis directory to use as input
                if length(P(i).runs) == 1
                    lower_level_directory = [params('rootdir') 'freesurfer/fsaverage/fla_matlab/' subjid '/' P(i).runtype '_r' num2str(P(i).runs) '/' P(i).lower_level_directory_name '_downsampled/contrasts/'];
                else
                    if length(P(i).runs) > 50
                        run_string = ['_' num2str(length(P(i).runs)) 'r-' num2str(P(i).runs(1)) '-' num2str(P(i).runs(end)) '_' DataHash(P(i).runs)];
                    else
                        run_string = ['_r' sprintf('%d',P(i).runs)];
                    end
                    lower_level_directory = [params('rootdir') 'freesurfer/fsaverage/sla_matlab/' subjid '/' P(i).runtype run_string '/' P(i).lower_level_directory_name '_downsampled/contrasts/'];
                end
                
                switch statistical_test
                    case 'fixed'
                        cope_file = [lower_level_directory  'cope_white_' P(i).contrast '.mat'];
                        cope_var_file = [lower_level_directory  'cope_var_white_' P(i).contrast '.mat'];
                    case {'random','fixed_without_whitening'}
                        cope_file = [lower_level_directory  'cope_' P(i).contrast '.mat'];
                        cope_var_file = [lower_level_directory 'cope_var_' P(i).contrast '.mat'];
                    otherwise
                        error('"fixed_or_random_effects" must be "fixed", "random", or "fixed_without_whitening"');
                end
                
            otherwise
                error('Error in sla_matlab: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"');
        end
        
        switch volume_or_surface
            
            case {'volume','surface'}
                % actual contrasts and variances
                x = MRIread(cope_file);
                if i == 1
                    br = x;
                    lower_level_cope = nan(n_subjects,numel(br.vol));
                    lower_level_cope_var = nan(n_subjects,numel(br.vol),1);
                    lower_level_fla_df = nan(n_subjects,1);
                end
                lower_level_cope(i,:) = x.vol(:);
                
                x = MRIread(cope_var_file);
                lower_level_cope_var(i,:) = x.vol(:);
                
                x = load([lower_level_directory 'df.mat']);
                lower_level_fla_df(i) = x.df;
                
            case 'downsampled_surface'
                
                % cope
                load(cope_file);
                if i == 1
                    n_voxels = numel(G.grid_data{1}) + numel(G.grid_data{2});
                    lower_level_cope = nan(n_subjects,n_voxels);
                    lower_level_cope_var = nan(n_subjects,n_voxels);
                    lower_level_fla_df = nan(n_subjects,1);
                end
                lower_level_cope(i,:) = [G.grid_data{1}(:)', G.grid_data{2}(:)'];
                
                % var cope
                load(cope_var_file);
                lower_level_cope_var(i,:) = [G.grid_data{1}(:)', G.grid_data{2}(:)'];
                
                % degrees of freedom
                x = load([lower_level_directory 'df.mat']);
                lower_level_fla_df(i) = x.df;
                
            otherwise
                error('Error in sla_matlab: volume_or_surface flag must be either "volume", "surface" or "downsampled_surface"');
        end
        
    end
    
    % divide by variances
    xi = lower_level_cope ~= 0 & lower_level_cope_var ~=0;
    lower_level_cope(xi) = lower_level_cope(xi) ./ sqrt(lower_level_cope_var(xi));
    lower_level_cope_var(xi) = 1;
    
    switch statistical_test
        case {'fixed', 'fixed_without_whitening'}
            tla_cope = mean(lower_level_cope);
            tla_cope_var = mean(lower_level_cope_var) / n_subjects;
            tla_df = sum(lower_level_fla_df);
        case 'random'
            tla_cope = mean(lower_level_cope);
            tla_cope_var = var(lower_level_cope,1) / (n_subjects-1);
            tla_df = n_subjects-1;
        otherwise
            error('"fixed_or_random_effects" must be "fixed", "random", or "fixed_without_whitening"');
    end
    
    tla_cope(any(lower_level_cope==0)) = 0;
    tla_cope_var(any(lower_level_cope_var==0)) = 0;

    % t-stat, df, and p-stat
    tstat = tla_cope ./ sqrt(tla_cope_var);
    pstat = -sign(tstat).*log10(2*tpvalue_copy(-abs(tstat), tla_df)); % two-tailed
    tstat(tla_cope_var==0) = 0;
    pstat(tla_cope_var==0) = 0;
    
    switch volume_or_surface
        
        case {'volume','surface'}
            
            br_tstat = br;
            br_tstat.vol = reshape(tstat(:),size(br.vol));
            br_tstat.fspec = t_file;
            MRIwrite(br_tstat, t_file);
            
            br_pstat = br;
            br_pstat.vol = reshape(pstat(:),size(br.vol));
            br_pstat.fspec = p_file;
            MRIwrite(br_pstat, p_file);
            
            br_cope = br;
            br_cope.vol = reshape(tla_cope(:),size(br.vol));
            br_cope.fspec = tla_cope_file;
            MRIwrite(br_cope, tla_cope_file);
            
            br_cope_var = br;
            br_cope_var.vol = reshape(tla_cope_var(:),size(br.vol));
            br_cope_var.fspec = tla_cope_var_file;
            MRIwrite(br_cope_var, tla_cope_var_file);
            
        case 'downsampled_surface'
            
            G.grid_data{1}(:) = pstat(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = pstat(numel(G.grid_data{1})+1:end);
            save(p_file, 'G');
            
            G.grid_data{1}(:) = tstat(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = tstat(numel(G.grid_data{1})+1:end);
            save(t_file, 'G');
            
            G.grid_data{1}(:) = tla_cope(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = tla_cope(numel(G.grid_data{1})+1:end);
            save(tla_cope_file, 'G');
            
            G.grid_data{1}(:) = tla_cope_var(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = tla_cope_var(numel(G.grid_data{1})+1:end);
            save(tla_cope_var_file, 'G');
            
            
        otherwise
            error('Error in sla_matlab: volume_or_surface flag must be either "volume", "surface" or "downsampled_surface"');
    end    
        
end

stat_range = [3,4.5,6];
if optInputs(varargin, 'stat_range')
    stat_range = varargin{optInputs(varargin, 'stat_range')+1};
end

stat_file = p_file;
if optInputs(varargin, 'tstat')
    stat_file = t_file;
end

if optInputs(varargin, 'sigmap')
    switch volume_or_surface
        case 'volume'
            if optInputs(varargin, 'monkey')
                unix_freesurfer_version('5.3.0',['freeview ' highres_2mm ':grayscale=0,150 ' stat_file ':colormap=heat:heatscale=' num2str(stat_range(1)) ',' num2str(stat_range(2)) ',' num2str(stat_range(3)) '  &']);
            else
                unix_freesurfer_version('5.3.0',['freeview ' highres_2mm ':grayscale=0,1000 ' stat_file ':colormap=heat:heatscale=' num2str(stat_range(1)) ',' num2str(stat_range(2)) ',' num2str(stat_range(3)) '  &']);
            end
        case 'surface'
            if optInputs(varargin, 'truncate')
                br = MRIread(stat_file);
                br.vol(br.vol<0) = 0;
                br.fspec = strrep(stat_file, '.mgz', '_truncate.mgz');
                MRIwrite(br, br.fspec);
                stat_file = br.fspec;
            end
            subjid = [exp '_us' num2str(us)];
            if optInputs(varargin, 'monkey')
                freeview3(subjid,hemi,'overlay',stat_file,'overlay_threshold',stat_range,varargin{:});
            else
                freeview3('fsaverage',hemi,'overlay',stat_file,'overlay_threshold',stat_range,varargin{:});
            end
        case 'downsampled_surface'
            load(p_file);
            figure;
            subplot(1,2,1);
            imagesc(flipud(rot90(G.grid_data{1})));
            title('Right Hemi');
            subplot(1,2,2);
            imagesc(fliplr(flipud(rot90(G.grid_data{2})))); %#ok<FLUDLR>
            title('Left Hemi');
            colorbar;
            
        otherwise
            error('Error in sla_matlab: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"');
    end
end

