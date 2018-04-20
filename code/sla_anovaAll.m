function sla_anovaAll(exp,usubs,runtype,condlevels,model,varargin)

% function slaAll(usubs,model,varargin)
%
% wrapper function for sla.m

for i = 1:length(usubs);
    
    allruns = read_runs(exp,usubs(i),runtype,'contrast',condlevels{1},model,varargin{:});
    runset = {allruns};
    
    % leave one out analysis for localizer runs
    %     for l = 1:length(allruns);
    %         runset = [runset, setdiff(allruns,allruns(l))]; %#ok<AGROW>
    %     end
    
    for l = 1:length(runset)
        fprintf('%s, sub %d, %s, runs %s\n',exp,usubs(i),runtype,sprintf('%d',runset{l}));
        sla_anova(exp,usubs(i),runtype,runset{l},condlevels,model,varargin{:});
    end
end