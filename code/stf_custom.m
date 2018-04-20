function stf_custom(exp,us,runtype,r,model,varargin)

% function stf_custom(exp,us,runtype,r,model,varargin)
%
% convolves boxcar with custom hrf to generate timing files for FSL

% sampling rate for convolution
sr = 5;

% type of hrf function to use
% current options are 'BOLD', 'MION', 'MION_CUSTOM1'
hrf_type = read_scanparams(exp,us,runtype,'run',r,'hrf_type',varargin{:});

% stf directory
stfdir = [params('rootdir') exp '/analysis/stf/usub' num2str(us) '/'];
if ~exist(stfdir, 'dir');
    mkdir(stfdir);
end

% misc parameters
conds = read_conditions(exp,us,runtype,'run',r,varargin{:});
b = read_timings(exp,us,runtype,r,varargin{:});
TR = read_scanparams(exp,us,runtype,'run',r,'TR',varargin{:});
nTR = read_scanparams(exp,us,runtype,'run',r,'nTR',varargin{:});
[h, t] = hrf_fsfast_gamma(1/sr, hrf_type, 'noplot');

% if strcmp(runtype,'pitchloc2_combined_raw') && r == 2
%   keyboard;
% end

% upsampled condition vector
totaltime = nTR*TR;
nsmps = ceil(totaltime*sr) + 1;
seq_upsampled = zeros(nsmps,1);
t = (0:nsmps-1)/sr;
for i = 1:length(b.onsets)
    if ~strcmp(b.conds{i}, 'NULL')
        xi = t >= b.onsets(i) & t < b.onsets(i)+b.durs(i);
        seq_upsampled(xi) = find(strcmp(b.conds{i}, conds));
    end
end

% upsampled hrf, one per condition
model = zeros(nsmps,length(conds));
for i = 1:length(conds)
    boxcar = zeros(nsmps,1);
    boxcar(seq_upsampled==i) = 1;
    
    x = conv(boxcar,h);
    model(:,i) = x(1:nsmps)/abs(sum(h));
end

% interpolate upsampled HRF
model_scans = interp1(t', model, (0:nTR-1)'*TR);

% demean
model_scans = model_scans - ones(nTR,1)*mean(model_scans);

try
    % plot
    if ~optInputs(varargin, 'noplot')
        figure;
        set(gcf, 'Position',[0 0 1440 700]/2);
        plot((0:nTR-1)'*TR, model_scans);
        ylim([min(model_scans(:)), max(model_scans(:))]);
        ylabel('Response');
        xlabel('Time (s)');
        legend(conds,'Location','NorthEastOutside','Orientation','Vertical');
        export_fig([stfdir runtype '_r' num2str(r) '.pdf'],'-pdf','-nocrop');
    end
catch
    keyboard;
end

% write to file
for i = 1:length(conds)
    regfile = [stfdir runtype '_r'  num2str(r) '_' conds{i} '_custom.stf'];
    fid = fopen(regfile,'w');
    fprintf(fid,repmat('%.6f\n',1,nTR),model_scans(:,i));
    fclose(fid);
end

% optionally add polynomial regressors to model
if optInputs(varargin, 'polynomial_confound_regressors')
    polyorder = varargin{optInputs(varargin, 'polynomial_confound_regressors')+1};
    t = zscore((1:nTR)');
    t_matrix = (t*ones(1,polyorder)) .^ (ones(nTR,1)*(1:polyorder));
    regfile = [stfdir runtype '_r'  num2str(r) '_confound_regressors.stf'];
    fid = fopen(regfile,'w');
    for i = 1:nTR
        for j = 1:polyorder
            fprintf(fid,'%.6f ',t_matrix(i,j));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end

%% Old Code
%
%     plot(t,model,'Color',colors{mod(j-1,length(colors))+1}, 'LineWidth',2);
%



% [~, ~, TR, TA, ~, ~, ~] = read_scanparams(exp,us,runtype,'run',r,varargin{:});
% delay = TR/2 - TA/2;




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

