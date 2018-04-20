function stfAll(exp,usubs,varargin)

% stfAll(usubs,varargin)
%
% wrapper for stfAll.m

for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp,usubs(i),'func',varargin{:});
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    for j = 1:length(runtypes)
        
        runnum = read_runs(exp,usubs(i),runtypes{j},varargin{:});
        if optInputs(varargin, 'runnum');
            runnum = varargin{optInputs(varargin,'runnum')+1};
        end

        models = read_models(exp,usubs(i),runtypes{j});
        
        for k = 1:length(runnum)
            for l = 1:length(models)
                
                fprintf('%s, usub %d, %s, run %d, %s\n',exp,usubs(i),runtypes{j},runnum(k),models{l});
                if optInputs(varargin, 'custom_stf')
                    stf_custom(exp,usubs(i),runtypes{j},runnum(k),models{l},varargin{:});
                elseif optInputs(varargin, 'custom_stf_v2')
                    stf_custom_v2(exp,usubs(i),runtypes{j},runnum(k),models{l},varargin{:});
                else
                    stf(exp,usubs(i),runtypes{j},runnum(k),models{l},varargin{:});
                end
            end
        end
    end
end