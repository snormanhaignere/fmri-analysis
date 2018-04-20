function stf(exp,us,runtype,r,model,varargin)

% function stf(us,runtype,r,model,varargin)
%
% generates an stf file for a particular subject, runtype, run, and model type
% this is based on experiment file output from main experiment (read by read_timings.m)

stfdir = [params('rootdir') exp '/analysis/stf/usub' num2str(us) '/'];
if ~exist(stfdir, 'dir');
    mkdir(stfdir);
end

conds = read_conditions(exp,us,runtype,'run',r,varargin{:});
b = read_timings(exp,us,runtype,r,'nofix');

[~, ~, TR, TA, ~, ~, ~] = read_scanparams(exp,us,runtype,'run',r,varargin{:});
delay = TR/2 - TA/2;

for i = 1:length(conds)
    ci = strcmp(conds{i},b.conds);
    regfile = [stfdir runtype '_r'  num2str(r) '_' conds{i} '.stf'];
    stfmat = [b.onsets(ci)+delay, b.durs(ci), ones(sum(ci),1)]; %b.durs(ci)

    fid = fopen(regfile,'w');
    fprintf(fid,repmat('%8.2f %8.2f %5d\n',1,size(stfmat,1)),stfmat');
    fclose(fid);
end

% elseif strcmp(model,'event')
%
%     for j = 1:length(b.condstims)
%
%         regfile = [stfdir 'sub' num2str(s) '/' runtype '_r' num2str(r) '_' b.condstims{j} '.stf'];
%         stfmat = [round(b.onsets(j)) triallength 1];
%
%         fid = fopen(regfile,'w');
%         fprintf(fid,repmat('%5d %5d %5d\n',1,size(stfmat,1)),stfmat');
%         fclose(fid);
%
%     end

% conditions = {'intact','pitch','rhythm','both'};
% for c = 1:length(conditions)
%     condorder = strmatch(conditions{c},data{6});
%     stfmat = [onsets(condorder) durations(condorder) ones(length(condorder),1)];
%     if ~exist([stfdir 'sub' num2str(s)], 'dir'); mkdir([stfdir 'sub' num2str(s)]); end
%     fid = fopen([stfdir 'sub' num2str(s) '/' 'timings_mus_s' num2str(s) '_r' num2str(runnum) '_' conditions{c} '.stf'],'w');
%     fprintf(fid,repmat('%5d %5d %5d\n',1,size(stfmat,1)),stfmat');
%     fclose(fid);
% end

