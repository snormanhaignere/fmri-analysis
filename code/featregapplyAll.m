function featregapplyAll(subjnum,model,varargin)

mypath;
runset = false;
if optInputs(varargin, 'runnum');
    runnum = varargin{optInputs(varargin,'runnum')+1};
    runset = true;
end

for s = 1:length(subjnum)
    if ~runset; runnum = runsread(subjnum(s),'main'); end
    for r = 1:length(runnum)
        fprintf('subject %d, run %d\n',subjnum(s), runnum(r));
        featregapply(subjnum(s),runnum(r),model);
    end
end
