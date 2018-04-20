function p_file = sla_matlab(exp,us,runtype,contrast,runs,volume_or_surface,fla_directory_name,sla_directory_name,fixed_or_random_effects,varargin)

if length(runs) > 50
    run_string = ['_' num2str(length(runs)) 'r-' num2str(runs(1)) '-' num2str(runs(end)) '_' DataHash(runs)];
else
    run_string = ['_r' sprintf('%d',runs)];
end

% output files
switch volume_or_surface
    case 'volume'
        
        sla_directory = [params('rootdir') exp '/analysis/sla_matlab/usub' num2str(us) '/' runtype run_string '/' sla_directory_name '/'];
        if ~exist(sla_directory,'dir');
            mkdir(sla_directory);
        end
        t_file =  [sla_directory 'tstat_' contrast '_' fixed_or_random_effects '.nii.gz'];
        p_file =  [sla_directory 'pstat_' contrast '_' fixed_or_random_effects '.nii.gz'];
        sla_cope_file =  [sla_directory 'cope_' contrast '_' fixed_or_random_effects '.nii.gz'];
        sla_cope_var_file =  [sla_directory 'cope_var_' contrast '_' fixed_or_random_effects '.nii.gz'];
        
        highres_2mm = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/struct_r1/brain_2mm.nii.gz']; % brain extracted structural image
        if ~exist([sla_directory 'highres_2mm.nii.gz'],'file') || optInputs(varargin, 'overwrite')
            copyfile(highres_2mm, [sla_directory 'highres_2mm.nii.gz'],'f');
        end
        
    case 'surface'
        subjid = [exp '_us' num2str(us)];
        if optInputs(varargin, 'monkey')
            sla_directory = [params('rootdir') 'freesurfer/' subjid '/sla_matlab/' runtype run_string '/' sla_directory_name '/'];
        else
            sla_directory = [params('rootdir') 'freesurfer/fsaverage/sla_matlab/' subjid '/' runtype run_string '/' sla_directory_name '/'];
        end
        if ~exist(sla_directory,'dir');
            mkdir(sla_directory);
        end
        if ~optInputs(varargin, 'hemi')
            error('Error in fit_glm: Need to specify hemisphere as optional argument for surface data.m');
        end
        hemi = varargin{optInputs(varargin, 'hemi')+1};
        
        t_file =  [sla_directory hemi '.tstat_' contrast '_' fixed_or_random_effects '.mgz'];
        p_file =  [sla_directory hemi '.pstat_' contrast '_' fixed_or_random_effects '.mgz'];
        sla_cope_file =  [sla_directory hemi '.cope_' contrast '_' fixed_or_random_effects '.mgz'];
        sla_cope_var_file =  [sla_directory hemi '.cope_var_' contrast '_' fixed_or_random_effects '.mgz'];
        
    case 'downsampled_surface'
        subjid = [exp '_us' num2str(us)];
        if optInputs(varargin, 'monkey')
            sla_directory = [params('rootdir') 'freesurfer/' subjid '/sla_matlab/' runtype run_string '/' sla_directory_name '_downsampled/'];
        else
            sla_directory = [params('rootdir') 'freesurfer/fsaverage/sla_matlab/' subjid '/' runtype run_string '/' sla_directory_name '_downsampled/'];
        end
        if ~exist(sla_directory,'dir');
            mkdir(sla_directory);
        end
        
        t_file =  [sla_directory  'tstat_' contrast '_' fixed_or_random_effects '.mat'];
        p_file =  [sla_directory  'pstat_' contrast '_' fixed_or_random_effects '.mat'];
        sla_cope_file =  [sla_directory  'cope_' contrast '_' fixed_or_random_effects '.mat'];
        sla_cope_var_file =  [sla_directory  'cope_var_' contrast '_' fixed_or_random_effects '.mat'];
        
    otherwise
        error('Error in sla_matlab: volume_or_surface flag must be either "volume" or "surface"...');
end

if ~exist(t_file,'file') || ~exist(p_file,'file') || ~exist(sla_cope_file,'file') || ~exist(sla_cope_var_file,'file') || optInputs(varargin, 'overwrite')
    
    % read in data
    
    
    source_directory = strrep(which('fla_matlab.m'),'fla_matlab.m','');
    addpath(genpath('/software/Freesurfer/5.3.0/matlab'));
    addpath([source_directory 'export_fig']);
    
    for i = 1:length(runs)
        % input files
        switch volume_or_surface
            case 'volume'
                fla_directory = [params('rootdir') exp '/analysis/fla_matlab/usub' num2str(us) '/' runtype '_r' num2str(runs(i)) '/' fla_directory_name '/contrasts/'];
                %                 fla_cope_file = [fla_directory 'cope_' contrast '_highres_2mm.nii.gz'];
                %                 fla_cope_var_file = [fla_directory 'cope_var_' contrast '_highres_2mm.nii.gz'];
                switch fixed_or_random_effects
                    case 'fixed'
                        fla_cope_file = [fla_directory 'cope_white_' contrast '.nii.gz'];
                        fla_cope_var_file = [fla_directory 'cope_var_white_' contrast '.nii.gz'];
                    case {'random','fixed_without_whitening'}
                        fla_cope_file = [fla_directory 'cope_' contrast '.nii.gz'];
                        fla_cope_var_file = [fla_directory 'cope_var_' contrast '.nii.gz'];
                    otherwise
                        error('"fixed_or_random_effects" must be "fixed", "random", or "fixed_without_whitening"');
                end
                
            case 'surface'
                subjid = [exp '_us' num2str(us)];
                if optInputs(varargin, 'monkey')
                    fla_directory = [params('rootdir') 'freesurfer/' subjid '/fla_matlab/' runtype '_r' num2str(runs(i)) '/' fla_directory_name '/contrasts/'];
                else
                    fla_directory = [params('rootdir') 'freesurfer/fsaverage/fla_matlab/' subjid '/' runtype '_r' num2str(runs(i)) '/' fla_directory_name '/contrasts/'];
                end
                
                switch fixed_or_random_effects
                    case 'fixed'
                        fla_cope_file = [fla_directory hemi '.cope_white_' contrast '.mgz'];
                        fla_cope_var_file = [fla_directory hemi '.cope_var_white_' contrast '.mgz'];
                    case {'random','fixed_without_whitening'}
                        fla_cope_file = [fla_directory hemi '.cope_' contrast '.mgz'];
                        fla_cope_var_file = [fla_directory hemi '.cope_var_' contrast '.mgz'];
                    otherwise
                        error('"fixed_or_random_effects" must be "fixed", "random", or "fixed_without_whitening"');
                end
            case 'downsampled_surface'
                subjid = [exp '_us' num2str(us)];
                if optInputs(varargin, 'monkey')
                    fla_directory = [params('rootdir') 'freesurfer/' subjid '/fla_matlab/' runtype '_r' num2str(runs(i)) '/' fla_directory_name '_downsampled/contrasts/'];
                else
                    fla_directory = [params('rootdir') 'freesurfer/fsaverage/fla_matlab/' subjid '/' runtype '_r' num2str(runs(i)) '/' fla_directory_name '_downsampled/contrasts/'];
                end
                
                switch fixed_or_random_effects
                    case 'fixed'
                        fla_cope_file = [fla_directory  'cope_white_' contrast '.mat'];
                        fla_cope_var_file = [fla_directory  'cope_var_white_' contrast '.mat'];
                    case {'random','fixed_without_whitening'}
                        fla_cope_file = [fla_directory  'cope_' contrast '.mat'];
                        fla_cope_var_file = [fla_directory 'cope_var_' contrast '.mat'];
                    otherwise
                        error('"fixed_or_random_effects" must be "fixed", "random", or "fixed_without_whitening"');
                end
            otherwise
                error('Error in sla_matlab: volume_or_surface flag must be either "volume", "surface", or "downsampled_surface"');
        end
        
        switch volume_or_surface
            
            case {'volume','surface'}
                % actual contrasts and variances
                x = MRIread(fla_cope_file);
                if i == 1
                    br = x;
                    fla_cope = nan(length(runs),numel(br.vol));
                    fla_cope_var = nan(length(runs),numel(br.vol),1);
                    fla_df = nan(length(runs),1);
                end
                fla_cope(i,:) = x.vol(:);
                
                x = MRIread(fla_cope_var_file);
                fla_cope_var(i,:) = x.vol(:);
                
                x = load([fla_directory 'df.mat']);
                fla_df(i) = x.df;
                
            case 'downsampled_surface'
                
                % cope
                load(fla_cope_file);
                if i == 1
                    n_voxels = numel(G.grid_data{1}) + numel(G.grid_data{2});
                    fla_cope = nan(length(runs),n_voxels);
                    fla_cope_var = nan(length(runs),n_voxels);
                    fla_df = nan(length(runs),1);
                end
                fla_cope(i,:) = [G.grid_data{1}(:)', G.grid_data{2}(:)'];
                
                % var cope
                load(fla_cope_var_file);
                fla_cope_var(i,:) = [G.grid_data{1}(:)', G.grid_data{2}(:)'];
                
                % degrees of freedom
                x = load([fla_directory 'df.mat']);
                fla_df(i) = x.df;
                
            otherwise
                error('Error in sla_matlab: volume_or_surface flag must be either "volume", "surface" or "downsampled_surface"');
        end
        
    end
    
    switch fixed_or_random_effects
        case {'fixed', 'fixed_without_whitening'}
            sla_cope = mean(fla_cope);
            sla_cope_var = mean(fla_cope_var) / length(runs);
            sla_df = sum(fla_df);
        case 'random'
            sla_cope = mean(fla_cope);
            sla_cope_var = var(fla_cope,1) / (length(runs)-1);
            sla_df = length(runs)-1;
        otherwise
            error('"fixed_or_random_effects" must be "fixed", "random", or "fixed_without_whitening"');
    end
    
    
    % t-stat, df, and p-stat
    tstat = sla_cope ./ sqrt(sla_cope_var);
    pstat = -sign(tstat).*log10(2*tpvalue_copy(-abs(tstat), sla_df)); % two-tailed
    tstat(sla_cope_var==0) = 0;
    pstat(sla_cope_var==0) = 0;
    
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
            br_cope.vol = reshape(sla_cope(:),size(br.vol));
            br_cope.fspec = sla_cope_file;
            MRIwrite(br_cope, sla_cope_file);
            
            br_cope_var = br;
            br_cope_var.vol = reshape(sla_cope_var(:),size(br.vol));
            br_cope_var.fspec = sla_cope_var_file;
            MRIwrite(br_cope_var, sla_cope_var_file);
            
        case 'downsampled_surface'
            
            G.grid_data{1}(:) = pstat(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = pstat(numel(G.grid_data{1})+1:end);
            save(p_file, 'G');
            
            G.grid_data{1}(:) = tstat(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = tstat(numel(G.grid_data{1})+1:end);
            save(t_file, 'G');
            
            G.grid_data{1}(:) = sla_cope(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = sla_cope(numel(G.grid_data{1})+1:end);
            save(sla_cope_file, 'G');
            
            G.grid_data{1}(:) = sla_cope_var(1:numel(G.grid_data{1}));
            G.grid_data{2}(:) = sla_cope_var(numel(G.grid_data{1})+1:end);
            save(sla_cope_var_file, 'G');
            
            
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

