function preprocAll(exp, usubs, varargin)

% if optInputs(varargin, 'monkey')
%     convert_monkey_func(exp, usubs, varargin{:});
%     convert_monkey_struct(exp, usubs, varargin{:})
% else
%     convertraw(exp,usubs);
%     reorientBXH(exp,usubs);
% end

for i = 1:length(usubs)
    
    %     runnum = read_runs(exp,usubs(i),'struct',varargin{:});
    %     if ~optInputs(varargin, 'monkey') && ~isempty(runnum)
    %         bet_struct(exp,usubs(i),varargin{:});
    %     end
    
    switch exp
        case {'pitch_overlap_v2'}
            combine_runtypes(exp, usubs(i), 'overlap_v2', {'overlap_v2_lowconds','overlap_v2_highconds'}, varargin{:});
        case {'pitch_overlap_v3'}
            combine_runtypes(exp, usubs(i), 'overlap_v3', {'overlap_v3_lowconds','overlap_v3_highconds'}, varargin{:});
        case {'naturalsound'}
            if usubs(i) == 1
                combine_runs(exp, usubs(i), 'main_v2', {1:4,5:8}, varargin{:});
                combine_runs(exp, usubs(i), 'main', {1:5,6:10}, varargin{:});
            end
            runnum = read_runs(exp, usubs(i), 'main_v3');
            if all(ismember(1:33, runnum))
                combine_runs(exp, usubs(i), 'main_v3', {1:11,12:22,23:33}, varargin{:});
                %                 combine_runs_raw(exp, usubs(i), 'main_v3', {1:11,12:22,23:33}, varargin{:});
                x = {1:11,12:22,23:33};
                for j = 1:3
                    combine_runs_without_reg(exp, usubs(i), 'main_v3', x{j}, 'main_v3_combined_raw', j, 'raw', varargin{:})
                end
            elseif all(ismember(1:22, runnum))
                combine_runs(exp, usubs(i), 'main_v3', {1:11,12:22}, varargin{:});
                %                 combine_runs_raw(exp, usubs(i), 'main_v3', {1:11,12:22}, varargin{:});
                x = {1:11,12:22};
                for j = 1:2
                    combine_runs_without_reg(exp, usubs(i), 'main_v3', x{j}, 'main_v3_combined_raw', j, 'raw', varargin{:})
                end
            elseif all(ismember(1:11, runnum))
                combine_runs(exp, usubs(i), 'main_v3', {1:11}, varargin{:});
                %                 combine_runs_raw(exp, usubs(i), 'main_v3', {1:11}, varargin{:});
                x = {1:11};
                for j = 1:1
                    combine_runs_without_reg(exp, usubs(i), 'main_v3', x{j}, 'main_v3_combined_raw', j, 'raw', varargin{:})
                end
            elseif usubs(i) == 136
            else
                fprintf('Error in preprocAll: not enough runs found.\n');
                drawnow;
                keyboard;
            end
            runnum = read_runs(exp, usubs(i), 'texture');
            if all(ismember(1:12, runnum))
                combine_runs(exp, usubs(i), 'texture', {1:12}, varargin{:});
                %                 combine_runs_raw(exp, usubs(i), 'texture', {1:12}, varargin{:});
                combine_runs_without_reg(exp, usubs(i), 'texture', 1:12, 'texture_combined_raw', 1, 'raw', varargin{:})
            end
            runnum = read_runs(exp, usubs(i), 'spectrotemporal');
            if all(ismember(1:10, runnum))
                combine_runs(exp, usubs(i), 'spectrotemporal', {1:10}, varargin{:});
                %                 combine_runs_raw(exp, usubs(i), 'texture', {1:12}, varargin{:});
                combine_runs_without_reg(exp, usubs(i), 'spectrotemporal', 1:10, 'spectrotemporal_combined_raw', 1, 'raw', varargin{:})
            end
            runnum = read_runs(exp, usubs(i), 'spectrotemporal_v2');
            if all(ismember(1:12, runnum))
                combine_runs(exp, usubs(i), 'spectrotemporal_v2', {1:12}, varargin{:});
                %                 combine_runs_raw(exp, usubs(i), 'texture', {1:12}, varargin{:});
                combine_runs_without_reg(exp, usubs(i), 'spectrotemporal_v2', 1:12, 'spectrotemporal_v2_combined_raw', 1, 'raw', varargin{:})
            end
            
    end
    
    runtypes = read_runtypes(exp, usubs(i), 'func', varargin{:});
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(runtypes)
        
        runnum = read_runs(exp,usubs(i),runtypes{j},varargin{:});
        if optInputs(varargin,'runnum')
            runnum = varargin{optInputs(varargin,'runnum')+1};
        end
        
        for k = 1:length(runnum)
            
            if ~strcmp('struct',runtypes{j})
                preproc(exp,usubs(i),runtypes{j},runnum(k),varargin{:});
            end
        end
    end
end