function bbregisterAll(exp,usubs,model,varargin)

% function bbregisterAll(usubs,model,varargin)
%
% wrapper for bbregister.m


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

        for r = 1:length(runnum)
            fprintf('%s, subject %d, %s, run %d\n',exp, usubs(i), runtypes{j}, runnum(r));
            bbregister(exp,usubs(i),runtypes{j},runnum(r),model,varargin{:});
        end
    end
end