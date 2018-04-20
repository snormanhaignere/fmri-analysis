function sla_matlabAll(exp,usubs,volume_or_surface,fla_directory_name,sla_directory_name,varargin)

% function slaAll(usubs,model,varargin)
%
% wrapper function for sla.m
model = 'block';
for i = 1:length(usubs);
    
    runtypes = read_runtypes(exp,usubs(i),'func');
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(runtypes)
        
        if strcmp(runtypes{j},'main_v3')
            contrasts = {'all_main_v3_runspec'};
            loc_contrasts = {};
        else
            [contrasts,loc_contrasts] = read_contrasts(exp,usubs(i),runtypes{j},model,varargin{:});
        end
        
        if optInputs(varargin,'contrasts');
            contrasts = varargin{optInputs(varargin,'contrasts')+1};
        end
        
        for k = 1:length(contrasts)
            
            allruns = read_runs(exp,usubs(i),runtypes{j},varargin{:});
            runset = {allruns};
            
            % leave one out analysis for localizer runs
            if ~optInputs(varargin, 'skip_loc_runs') && any(strcmp(contrasts{k},loc_contrasts));
                for l = 1:length(allruns);
                    runset = [runset, setdiff(allruns,allruns(l))]; %#ok<AGROW>
                end
            end
            
            % loop through sets of runs
            for l = 1:length(runset)
                if isempty(runset{l});
                    warning(['Skipping at least 1 run set of ' contrasts{k}]); %#ok<WNTAG>
                    continue;
                end
                fprintf('%s, sub %d, %s, runs %s, %s\n',exp,usubs(i),runtypes{j},sprintf('%d',runset{l}),contrasts{k});
                switch volume_or_surface
                    case 'volume'
                        if optInputs(varargin, 'fixed')
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed',varargin{:});
                        elseif optInputs(varargin, 'random')
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random',varargin{:});
                        else
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed',varargin{:});
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random',varargin{:});
                        end
                    case 'surface'
                        if optInputs(varargin, 'permutation')
                            sla_permutation(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'hemi','rh',varargin{:});
                            sla_permutation(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'hemi','lh',varargin{:});
                        elseif optInputs(varargin, 'fixed')
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed','hemi','rh',varargin{:});
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed','hemi','lh',varargin{:});
                        elseif optInputs(varargin, 'random')
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random','hemi','rh',varargin{:});
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random','hemi','lh',varargin{:});
                        else
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed','hemi','rh',varargin{:});
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed','hemi','lh',varargin{:});
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random','hemi','rh',varargin{:});
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random','hemi','lh',varargin{:});
                        end
                    case 'downsampled_surface'
                        if optInputs(varargin, 'permutation')
                            sla_permutation(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,varargin{:});
                        elseif optInputs(varargin, 'fixed_without_whitening')
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed_without_whitening',varargin{:});
                        elseif optInputs(varargin, 'fixed')
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed',varargin{:});
                        elseif optInputs(varargin, 'random')
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random',varargin{:});
                        else
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'fixed_without_whitening',varargin{:});
                            sla_matlab(exp,usubs(i),runtypes{j},contrasts{k},runset{l},volume_or_surface,fla_directory_name,sla_directory_name,'random',varargin{:});
                        end
                        
                    otherwise
                        error('"volume_or_surface" must be "volume", "surface", or "downsampled_surface"')
                end
            end
        end
    end
end