function bbregister_flirtfilesAll(subjnum,model,varargin)

% function bbregister_flirtfilesAll(subjnum,model,varargin)
% 
% wrapper for bbregister_flirtfiles.m

runset = false;
if optInputs(varargin, 'runnum');
    runnum = varargin{optInputs(varargin,'runnum')+1};
    runset = true;
end

runtypes = read_runtypes('func');
if optInputs(varargin,'runtypes')
    runtypes = varargin{optInputs(varargin,'runtypes')+1};
end

for s = 1:length(subjnum)
    for t = 1:length(runtypes)
        
        if ~runset; runnum = runsread(subjnum(s),'main'); end
        for r = 1:length(runnum)
            fprintf('subject %d, %s, run %d\n',subjnum(s), runtypes{t}, runnum(r));
            bbregister_flirtfiles(subjnum(s),runtypes{t},runnum(r),model,varargin{:});
        end
        
    end
end
