function slaAll(exp,usubs,model,varargin)

% function slaAll(usubs,model,varargin)
%
% wrapper function for sla.m

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
            
            allruns = read_runs(exp,usubs(i),runtypes{j},'contrast',contrasts{k},model,varargin{:});
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
                sla(exp,usubs(i),runtypes{j},runset{l},contrasts{k},model,varargin{:});
            end
        end
    end
end
