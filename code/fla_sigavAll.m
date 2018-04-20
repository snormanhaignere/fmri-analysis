function fla_sigavAll(exp,usubs,varargin)

% stfAll(exp,usubs,varargin{:});
% keyboard;
for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp,usubs(i),'func');
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(runtypes)
        
        runnum = read_runs(exp,usubs(i),runtypes{j},varargin{:});
        if optInputs(varargin, 'runnum');
            runnum = varargin{optInputs(varargin,'runnum')+1};
        end
        
        for k = 1:length(runnum)
            
            fprintf('%s, subject %d, %s, run %d\n',exp,usubs(i),runtypes{j},runnum(k));
            fla_sigav(exp,usubs(i),runtypes{j},runnum(k),varargin{:});
            
        end
    end
end