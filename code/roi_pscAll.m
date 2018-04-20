function [psc_mean psc_sem psc plotconds] = roi_pscAll(usubs, roi_name, test_exp, test_runtype, varargin)%loc_runtype, loc, loc_sign, loc_op, loc_opval, test_runtype, varargin)

% [psc_mean runs_used mainstr psc] = roi_plotconds(s, test_runtype, testmodel, masktype, roi, varargin)
% For anatomical + contrast: roi_plotconds(s, test_runtype, testmodel, masktype, roi, loc_runtype, locmodel, contrast)
%
% Basically a wrapper/plotting function for roi_psc to be used in single subjects
%
% Specificy a subject, a test runtype (e.g. 'main', 'localizer'), a test model (e.g. 'block' or 'event'), and then a type of mask: anat, anat_contrast, and cluster
%
% For anat you only need to specify an roi
% For cluster and anat_contrast, you need to specify an anatomical roi, a localizer contrast, a localizer runtype, & a localizer model (next two arguments respectively)
%
% Also returns the the mean psc value which is used in by roi_plotconds
%
% When test and localizer runs overlap, analyses are performed in a leave-one-out style
% See roi_psc.m for more details/options
%
% Sample:
% [psc_mean runs_used mainstr psc] = roi_plotconds(1, 'main', 'block', 'anat_contrast', 'destrieux_pp', 'harm-highfreq_vs_noise-highfreq', 'localizer', 'block')

%% paramaters, directories, etc

addpath(genpath('export_fig'));

% test parameters
clear test;
test.exp = test_exp;
test.runtype = test_runtype;
test.conds = read_conditions(test.exp, usubs(1), test.runtype);
test.model = 'block';
test.fwhm = read_smooth(test.exp, usubs(1), test.runtype, varargin{:});

% roi parameters
clear roi;
roi.name = roi_name;
roi.op = 'thr';
roi.opval = 0.5;
if optInputs(varargin, 'roi_op')
    roi.op = varargin{optInputs(varargin, 'roi_op') + 1};
    roi.opval = varargin{optInputs(varargin, 'roi_op') + 2};
end

% localizer parameters
loc = [];
if optInputs(varargin, 'loc')
    clear loc;
    inds = optInputs(varargin, 'loc');
    for i = 1:length(inds)
        loc(i).exp = varargin{inds(i)+1}; %#ok<*AGROW>
        loc(i).runtype = varargin{inds(i)+2};
        loc(i).con = varargin{inds(i)+3};
        loc(i).sign = varargin{inds(i)+4};
        loc(i).op = varargin{inds(i)+5};
        loc(i).opval = varargin{inds(i)+6};
        loc(i).model = 'block';
        loc(i).fwhm = read_smooth(loc(i).exp, usubs(1), loc(i).runtype, varargin{:});%read_smooth(loc(i).exp, usubs(1), loc(i).runtype, varargin{:});
    end
end


% parameters for signal averaging analysis
if optInputs(varargin,'sigav_plateau') || optInputs(varargin,'sigav_peaktime') || optInputs(varargin,'sigav_timecourse')
    
    % parameters and analysis window
    [blockdur, nulldur, TR, TA, stimdur, stim2scan, win] = read_scanparams(test.exp, usubs(1), test.runtype); %#ok<ASGLU>
    
    if optInputs(varargin,'sigav_peaktime');
        
        peak_time = varargin{optInputs(varargin,'sigav_peaktime')+1};
        
        % peak time point
        x = abs(win - peak_time);
        [zz, peak_smp] = min(x); peak_smp = peak_smp(1); %#ok<ASGLU>
        peak_time = win(peak_smp);
        
    elseif optInputs(varargin,'sigav_plateau');
        
        start_time = varargin{optInputs(varargin,'sigav_plateau')+1};
        if optInputs(varargin, 'blockdur')
            blockdur = varargin{optInputs(varargin,'blockdur')+1};
        end
        
        % plateau samples
        plat_smp = find(win > start_time - 1e-3 & win < blockdur-TR+stim2scan + 1e-3);
        plat_time = win(plat_smp);
        plat_time
        
    end
end

scripst_directory = [pwd '/'];
figure_directory = strrep(scripst_directory,'scripts/',[test.exp '/figures/roi/']);
if ~exist(figure_directory,'dir');
    mkdir(figure_directory);
end

%% main

if optInputs(varargin,'sigav_timecourse')
    psc_allsub = nan(length(usubs), length(win), length(test.conds));
else
    psc_allsub = nan(length(usubs), length(test.conds));
end

% loop through subjects
skipped_sub = [];
for i = 1:length(usubs)
    
    
    % runs to test on
    if optInputs(varargin, 'testruns') % specify test runs
        
        testruns = varargin{optInputs(varargin, 'testruns')+1};
        
    elseif optInputs(varargin, 'low-tsnr') % runs with low tsnr
        
        % tsnr values and threshold
        thresh = varargin{optInputs(varargin, 'low-tsnr')+1};
        tsnr = roi_tsnr(test.exp, usubs(i), roi.name, test.runtype, 'raw', 'nomask', 'allruns');
        [~,ti] = sort(tsnr,'ascend');
        
        % select runs
        x = read_runs(test.exp, usubs(i), test.runtype, varargin{:});
        testruns = x( ti( 1:round( thresh*length(x) ) ) );
        
    elseif optInputs(varargin, 'high-tsnr') % runs with high tsnr
        
        % tsnr values and threshold
        thresh = varargin{optInputs(varargin, 'low-tsnr')+1};
        tsnr = roi_tsnr(test.exp, usubs(i), roi.name, test.runtype, 'raw', 'nomask', 'allruns');
        [~,ti] = sort(tsnr,'descend');
        
        % select runs
        x = read_runs(test.exp, usubs(i), test.runtype, varargin{:});
        testruns = x( ti( 1:round( thresh*length(x) ) ) );
        
    elseif optInputs(varargin, 'tsnr-thresh') % runs that pass an absolute tsnr threshold
        
        % tsnr values and threshold
        thresh = varargin{optInputs(varargin, 'tsnr-thresh')+1};
        tsnr = roi_tsnr(test.exp, usubs(i), roi.name, test.runtype, 'raw', 'nomask', 'allruns');
        
        % select runs
        x = read_runs(test.exp, usubs(i), test.runtype, varargin{:});
        testruns = x(tsnr>thresh);
        
    else
        
        testruns = read_runs(test.exp, usubs(i), test.runtype, varargin{:}, 'roi_test_runs'); % all possible test runs
        
    end
    
    if optInputs(varargin, 'within_run_analysis')
        
        for k = 1:length(loc)
            loc(k).runs = read_runs(loc(k).exp, usubs(i), loc(k).runtype, varargin{:});
            test.runs = testruns;
            x = roi_psc_v2(usubs(i), roi, test, loc, varargin{:});
        end
        
        if any(isnan(x(:)));
            error('Bad run: NaNs found.');
        elseif optInputs(varargin,'sigav_timecourse')
            psc_sub = x;
        elseif optInputs(varargin,'sigav_peaktime')
            psc_sub = squeeze(x(:,peak_smp,:));
        elseif optInputs(varargin,'sigav_plateau')
            psc_sub = squeeze(mean(x(:,plat_smp,:),2));
        end
        
        if optInputs(varargin,'sigav_timecourse')
            psc_allsub(i,:,:) = squeeze(mean(psc_sub,1));
        else
            psc_allsub(i,:) = squeeze(mean(psc_sub,1));
        end
        
    else
        
        % initialize matrix to store psc values for a given subject
        if optInputs(varargin,'sigav_timecourse')
            psc_sub = nan(length(testruns), length(win), length(test.conds));
        else
            psc_sub = nan(length(testruns), length(test.conds));
        end
        
        % loop through test runs
        skipped_runs = [];
        for j = 1:length(testruns)
            
            test.run = testruns(j);
            
            % select nonindependent test runs
            if optInputs(varargin, 'nonindependent')
                warning('MATLAB:CUSTOM','Nonindependent Analysis!');
            end
            
            for k = 1:length(loc)
                if optInputs(varargin, 'locruns')
                    loc(k).runs = varargin{optInputs(varargin, 'locruns')+1};
                elseif optInputs(varargin, 'locruns_pool')
                    loc(k).runs = varargin{optInputs(varargin, 'locruns_pool')+1};
                    loc(k).runs = setdiff(loc(k).runs, testruns(j));
                elseif strcmp(loc(k).exp, test.exp) && strcmp(loc(k).runtype, test.runtype) && ~optInputs(varargin,'nonindependent')
                    x = read_runs(loc(k).exp, usubs(i), loc(k).runtype, varargin{:});
                    loc(k).runs = setdiff(x, testruns(j));
                else
                    loc(k).runs = read_runs(loc(k).exp, usubs(i), loc(k).runtype, varargin{:});
                end
            end
            
            if optInputs(varargin, 'surf')
                x = roi_psc_surf(usubs(i), roi, test, loc, varargin{:});
            elseif optInputs(varargin, 'matlab_volume')
                x = roi_psc_matlab_volume(usubs(i), roi, test, loc, varargin{:});
            elseif optInputs(varargin, 'volume')
                x = roi_psc_volume(usubs(i), roi, test, loc, varargin{:});
            else
                x = roi_psc(usubs(i), roi, test, loc, varargin{:});
            end
            
            if any(isnan(x(:)));
                skipped_runs = [skipped_runs j]; %#ok<AGROW>
                fprintf('Skipping run %d\n',testruns(j)); drawnow;
            elseif optInputs(varargin,'sigav_timecourse')
                psc_sub(j,:,:) = x;
            elseif optInputs(varargin,'sigav_peaktime')
                psc_sub(j,:) = x(peak_smp,:);
            elseif optInputs(varargin,'sigav_plateau')
                psc_sub(j,:) = mean(x(plat_smp,:),1);
            else
                psc_sub(j,:) = x;
            end
            
        end
        
        % remove skipped runs
        if ~isempty(skipped_runs) && optInputs(varargin,'sigav_timecourse')
            psc_sub(skipped_runs,:,:) = [];
            error('\nSkipped %d Runs\n', length(skipped_runs));
        elseif ~isempty(skipped_runs)
            psc_sub(skipped_runs,:) = [];
            error('\nSkipped %d Runs\n', length(skipped_runs));
        end
        
        if ~isempty(skipped_runs)
            skipped_sub = [skipped_sub i]; %#ok<AGROW>
        elseif optInputs(varargin,'sigav_timecourse')
            psc_allsub(i,:,:) = squeeze(mean(psc_sub,1));
        else
            psc_allsub(i,:) = squeeze(mean(psc_sub,1));
        end
    end

    
end

% remove skipped subjects
if ~isempty(skipped_sub) && optInputs(varargin,'sigav_timecourse')
    psc_allsub(skipped_sub,:,:) = [];
    error('\nSkipped Subjects: %s\n', sprintf('%d ',usubs(skipped_sub)));
elseif ~isempty(skipped_sub)
    psc_allsub(skipped_sub,:) = [];
    error('\nSkipped Subjects: %s\n', sprintf('%d ',usubs(skipped_sub)));
end

%% Standard Errors

if length(usubs) > 1
    psc = psc_allsub;
else
    psc = psc_sub;
end

if optInputs(varargin,'mean-norm');
    psc = psc ./ repmat( mean(psc,2), 1, size(psc,2) );
end

if isempty(psc)
    error('No valid subjects.');
    return; %#ok<UNRCH>
end


%% T-tests Direct Contrasts

test.contrasts = read_maincontrasts(test_exp, test_runtype, 'us', usubs(1), varargin{:});

if ~isempty(test.contrasts)
    
    if optInputs(varargin,'sigav_timecourse')
        psc_contrasts = nan(size(psc,1), length(win), length(test.contrasts));
    else
        psc_contrasts = nan(size(psc,1), length(test.contrasts));
    end
    
    for i = 1:length(test.contrasts)
        
        if optInputs(varargin,'sigav_timecourse')
            for j = 1:length(win)
                psc_contrasts(:,j,i) = psc_contr(squeeze(psc(:,j,:)), test_exp, test.contrasts{i}, test.conds, usubs, test.model, varargin{:});
            end
        else
            psc_contrasts(:,i) = psc_contr(psc, test_exp, test.contrasts{i}, test.conds, usubs, test.model, varargin{:});
            [~, pval,~,stats] = ttest(psc_contrasts(:,i));
            fprintf('%s: %8.6f, %8.6f\n',test.contrasts{i},pval,stats.tstat);
        end
    end
end

%% anovas

if optInputs(varargin,'anova2x2')
    switch test.exp
        case {'pitch_f0inharm', 'pitch_f0smallfreq'}
            
            data = psc_contrasts(:,[2 3 5 6]);%- repmat(mean(psc(:,[1 5 7 11]),2),1,4);
            f0fac = repmat(1:2,size(psc_contrasts,1),2);
            harmfac = repmat([1 1 2 2],size(psc_contrasts,1),1);
            
            subORrunfac = repmat((1:size(psc_contrasts,1))',1,4);
            anovan(data(:),{f0fac(:) harmfac(:) subORrunfac(:)},'model',2,'random',3,'varnames',{'F0 Adapt','Harmonicity','Runs/Subjects'});
            rm_anova2(data,subORrunfac,f0fac,harmfac,{'F0 Adapt','Harmonicity'})
    end
end


if optInputs(varargin,'anova2x3')
    switch test.exp
        case {'pitch_overlap_v3'}
            if optInputs(varargin, 'collapse_frequency')
                data = psc_contrasts(:,4:9);%- repmat(mean(psc(:,[1 5 7 11]),2),1,4);
                idfac = repmat(1:3,size(psc_contrasts,1),2);
                varfac = repmat([1 1 1 2 2 2],size(psc_contrasts,1),1);
                subORrunfac = repmat((1:size(psc_contrasts,1))',1,6);
                
                anovan(data(:),{idfac(:) varfac(:) subORrunfac(:)},'model',2,'random',3,'varnames',{'Note Identity','Note Variation','Runs/Subjects'});
                rm_anova2(data,subORrunfac,idfac,varfac,{'Note Identity','Note Variation'})
            end
    end
end

if optInputs(varargin,'anova3x3')
    switch test.exp
        case {'pitch_overlap_v3'}
            if optInputs(varargin, 'collapse_frequency')
                data = psc_contrasts;%- repmat(mean(psc(:,[1 5 7 11]),2),1,4);
                idfac = repmat(1:3,size(psc_contrasts,1),3);
                varfac = repmat([1 1 1 2 2 2 3 3 3],size(psc_contrasts,1),1);
                subORrunfac = repmat((1:size(psc_contrasts,1))',1,9);
                
                anovan(data(:),{idfac(:) varfac(:) subORrunfac(:)},'model',2,'random',3,'varnames',{'Note Identity','Note Variation','Runs/Subjects'});
                rm_anova2(data,subORrunfac,idfac,varfac,{'Note Identity','Note Variation'})
            end
    end
end

%% T-tests 2nd order contrasts
%
contrasts_2ndorder = read_maincontrasts_2ndorder(test_exp, test_runtype);
psc_contrasts_2ndorder = nan(size(psc,1), length(contrasts_2ndorder));

if ~optInputs(varargin, 'sigav_timecourse')
    
    for i = 1:length(contrasts_2ndorder)
        
        contrasts_1storder = regexp(contrasts_2ndorder{i},'_vs2_','split');
        
        contr1 = psc_contr(psc, test_exp, contrasts_1storder{1}, test.conds, usubs, test.model, varargin{:});
        contr2 = psc_contr(psc, test_exp, contrasts_1storder{2}, test.conds, usubs, test.model, varargin{:});
        psc_contrasts_2ndorder(:,i) = contr1-contr2;
        
        [~, pval] = ttest(psc_contrasts_2ndorder(:,i));
        fprintf('%s: %8.6f\n',contrasts_2ndorder{i},pval);
        
    end
end

%%

if optInputs(varargin, 'plotcontrasts')
    
    psc = psc_contrasts;
    plotconds = test.contrasts;
    plotorder = test.contrasts;
    
elseif optInputs(varargin, 'plotcontrasts_2ndorder')
    
    psc = psc_contrasts_2ndorder;
    plotconds = contrasts_2ndorder;
    plotorder = contrasts_2ndorder;
    
else
    
    plotconds = test.conds;
    plotorder = read_plotorder(test.exp, usubs(1), test.runtype);
    
end

if optInputs(varargin, 'flipdata')
    psc = -psc;
end

%% Title

if optInputs(varargin, 'withsub_sem')
    semstring = 'Within-Subject SEM';
else
    semstring = 'SEM';
end

if isempty(loc)
    titlestring = [sprintf('\n') strrep(test_exp, '_','') ' sub ' sprintf('%d',usubs) sprintf('\n') 'Atlas: ' strrep(roi.name,'_',' ') sprintf('\n') semstring];
else
    
    locstring = [];
    
    for i = 1:length(loc)
        
        if i > 1
            locstring = [locstring sprintf('\n')];
        end
        
        a = regexp(loc(i).con,'_vs_','split');
        if length(a) == 1;
            a = [a, 'NULL'];
        end
        
        b = cell(1,2);
        for j = 1:length(a)
            
            x = read_plotformat(a{j});
            b{j} = x.sym;
            
        end
        
        switch loc(i).sign
            case 'pos'
                locstring = [locstring b{1} ' > ' b{2}];
            case 'neg'
                locstring = [locstring b{2} ' > ' b{1}];
            case 'abs'
                locstring = [locstring b{1} ' > or < ' b{2}];
            otherwise
                error('Bad sign');
        end
        
        
        if strcmp(loc(i).op,'nvox')
            locstring = [locstring ', Vox ' num2str(loc(i).opval(1)) ' to ' num2str(loc(i).opval(2))];
        elseif strcmp(loc(i).op,'perc')
            locstring = [locstring, sprintf(', Top %.2f to %.2f%% of voxels', 100*loc(i).opval(1), 100*loc(i).opval(2))];
        elseif strcmp(loc(i).op,'thr')
            locstring = [locstring sprintf(', P < %.3f', 1-normcdf(loc(i).opval))];
        elseif strcmp(loc(i).op,'negthr')
            locstring = [locstring sprintf(', P > %.3f', 1-normcdf(loc(i).opval))];
        else
            error('Need to add title string for this loc_op');
        end
        
    end
    
    titlestring = [sprintf('\n') strrep(test_exp, '_','') ' s' sprintf('%d',usubs) sprintf('\n') 'Atlas: ' strrep(roi.name,'_',' ') sprintf('\n') locstring sprintf('\n') semstring];% strrep(test_runtype, '_', ' ')  ', ' strrep(masktype,'_',' ') ', '
end

%% Plot

% Mean PSC value
psc_mean = squeeze(mean(psc,1));

% Standard Error
if optInputs(varargin,'withsub_sem')
    [psc_sem psc_norm] = stderr_withsub_corrected(psc);
elseif optInputs(varargin,'bootstrap')
    [lsem usem] = bootstrap_mean(psc);
elseif optInputs(varargin,'bootstrap_withsub')
    [~, psc_norm] = stderr_withsub_corrected(psc);
    [x y] = bootstrap_mean(psc_norm);
    lsem = x + mean(psc(:));
    usem = y + mean(psc(:));
else
    psc_sem = squeeze(std(psc,[],1)/sqrt(size(psc,1)-1));
    psc_norm = zeros(size(psc_sem));
end

% Return if Don't want to plot
if optInputs(varargin,'noplot')
    return;
end

%%

% Plot
if optInputs(varargin,'sigav_timecourse')
    if size(psc_mean,1) == 1
        psc_mean = psc_mean';
        psc_sem = psc_sem';
    end
    myplot(win',psc_mean,plotconds,plotorder,'errorbar',psc_sem);
elseif optInputs(varargin, 'scatter')
    %     myscatter(mean(psc(:)) + psc_norm, plotconds, plotorder);
    %     keyboard;
    myscatter(psc, plotconds, plotorder);
elseif optInputs(varargin, 'bootstrap') || optInputs(varargin, 'bootstrap_withsub')
    mybar(psc_mean,plotconds,plotorder,'errorbar2',lsem,usem,varargin{:});
else
    mybar(psc_mean,plotconds,plotorder,'errorbar',psc_sem,varargin{:});
end

if ~optInputs(varargin, 'simplot')
    title(titlestring);
end
% drawnow;
% keyboard;
runs_string = '';
if optInputs(varargin, 'testruns')
    x = varargin{optInputs(varargin, 'testruns')+1};
    runs_string = [runs_string, '_testruns' num2str(min(testruns)) '-' num2str(max(testruns))];
end

if optInputs(varargin, 'scans')
    x = varargin{optInputs(varargin, 'scans')+1};
    runs_string = [runs_string, '_scans' sprintf('%d', x)];
end

set(gcf, 'Position', [0 0 1000, 900]);
% keyboard;
if ~isempty(loc)
    export_fig([figure_directory 'us' sprintf('%d',usubs) '_'  roi.name '_' loc(1).con '_' loc(1).sign runs_string  '_' DataHash(varargin)  '.pdf'], '-pdf', '-nocrop')
else
    export_fig([figure_directory 'us' sprintf('%d',usubs) '_'  roi.name '_' runs_string  '_' DataHash(varargin)  '.pdf'], '-pdf', '-nocrop')
end
% titlestring = ;

%%


% contrasts = read_maincontrasts;
% psc_contrasts = nan(size(psc,1), length(contrasts));
% psc_contrasts_weights = nan(size(psc,1), length(contrasts));
%
% if ~optInputs(varargin, 'sigav_timecourse')
%     for i = 1:length(contrasts)
%         psc_contrasts(:,i) = psc_contr(psc, contrasts{i}, conds, usubs, testmodel, varargin{:});
%
%         if optInputs(varargin, 'weightmean')
%             for j = 1:length(usubs)
%                 x = psc_contr(psc_allsub_allrun{j}, contrasts{i}, conds, usubs(j), testmodel, varargin{:});
%                 sem = std(x)/sqrt(size(x,1)-1);
%                 psc_contrasts_weights(j,i) = 1/(sem.^2);
%             end
%             [~, pval] = ttest(psc_contrasts(:,i).*psc_contrasts_weights(:,i));
%             fprintf('%s: %8.6f\n',contrasts{i},pval);
%         else
%             [~, pval] = ttest(psc_contrasts(:,i));
%             fprintf('%s: %8.6f\n',contrasts{i},pval);
%         end
%     end
% end

%% T-tests 2nd order contrasts

% contrasts_2ndorder = read_maincontrasts_2ndorder;
% psc_contrasts_2ndorder = nan(size(psc,1), length(contrasts_2ndorder));
% psc_contrasts_2ndorder_weights = nan(size(psc,1), length(contrasts));
%
% if ~optInputs(varargin, 'sigav_timecourse')
%
%     for i = 1:length(contrasts_2ndorder)
%
%         contrasts_1storder = regexp(contrasts_2ndorder{i},'_vs2_','split');
%
%         contr1 = psc_contr(psc, contrasts_1storder{1}, conds, usubs, testmodel, varargin{:});
%         contr2 = psc_contr(psc, contrasts_1storder{2}, conds, usubs, testmodel, varargin{:});
%         psc_contrasts_2ndorder(:,i) = contr1-contr2;
%
%         if optInputs(varargin, 'weightmean')
%             for j = 1:length(usubs)
%                 x1 = psc_contr(psc_allsub_allrun{j},  contrasts_1storder{1}, conds, usubs(j), testmodel, varargin{:});
%                 x2 = psc_contr(psc_allsub_allrun{j},  contrasts_1storder{2}, conds, usubs(j), testmodel, varargin{:});
%                 sem = std(x1-x2)/sqrt(size(x1,1)-1);
%                 psc_contrasts_2ndorder_weights(j,i) = 1/(sem.^2);
%             end
%             [~, pval] = ttest(psc_contrasts_2ndorder(:,i).*psc_contrasts_2ndorder_weights(:,i));
%             fprintf('%s: %8.6f\n',contrasts_2ndorder{i},pval);
%         else
%             [~, pval] = ttest(psc_contrasts_2ndorder(:,i));
%             fprintf('%s: %8.6f\n',contrasts_2ndorder{i},pval);
%         end
%
%     end
% end
%
% if optInputs(varargin, 'plotcontrasts')
%
%     psc = psc_contrasts;
%     plotconds = contrasts;
%     plotorder = contrasts;
%     if optInputs(varargin, 'weightmean')
%         psc = psc_contrasts.*psc_contrasts_weights./repmat(sum(psc_contrasts_weights),size(psc_contrasts,1),1);
%     end
%
% elseif optInputs(varargin, 'plotcontrasts_2ndorder')
%
%     psc = psc_contrasts_2ndorder;
%     plotconds = contrasts_2ndorder;
%     plotorder = contrasts_2ndorder;
%
% else
%
%     plotconds = conds;
%     plotorder = read_plotorder(usubs(1), test_runtype);
%     if optInputs(varargin,'sigav_timecourse');
%         plotorder = read_plotorder(usubs(1), test_runtype, '2groups');
%     end
%
% end


%% Plot (OLD)
%
% if optInputs(varargin,'plotrescontrasts');
%
%     mybar(res_contrasts,contrasts,contrasts,varargin{:});
%
% elseif optInputs(varargin, 'psc_vs_res')
%
%     cmap = [...
%         0.1190    0.1386    0.3517
%         0.4984    0.1493    0.8308
%         0.9597    0.2575    0.5853
%         0.3404    0.8407    0.5497
%         0.5853    0.2543    0.9172
%         0.2238    0.8143    0.2858
%         0.7513    0.2435    0.7572
%         0.2551    0.9293    0.7537
%         0.5060    0.3500    0.3804
%         0.6991    0.1966    0.5678
%         0.8909    0.2511    0.0759
%         0.9593    0.6160    0.0540
%         0.5472    0.4733    0.5308];
%     inds = setdiff(1:13,5);
%     load('peaks.mat')
%     xdat = p.statmean{4}(inds);
%     ydat = mean(psc(:,inds));
%     psc_sem = stderr_withsub(psc(:,inds));
%     hold on;
%     for i = 1:length(xdat)
%         errorbar(xdat(i),ydat(i),psc_sem(i),'o','Color',cmap(i,:),'MarkerSize',10,'LineWidth',2);
%     end
%     legend(conds(inds));
%     xbounds = [min(xdat)-1, max(xdat)+1];
%     ybounds = [min(ydat)-0.25, max(ydat)+0.25];
%     xlim(xbounds);ylim(ybounds);
%     c = polyfit(xdat(1:4),ydat(1:4),1);
%     plot(linspace(xbounds(1),xbounds(2),20),linspace(xbounds(1),xbounds(2),20)*c(1)+c(2),'k--');
%     title(titlestring);
%     ylabel('% Signal Change');
%     xlabel('Resolvability');
%     hold off;
%
% elseif optInputs(varargin, 'psc_vs_res_contrasts')
%
%     xdat = res_contrasts;
%     ydat = mean(psc_contrasts,1);
%     plot(xdat,ydat,'o')
%     psc_sem = squeeze(std(psc_contrasts,[],1)/sqrt(size(psc_contrasts,1)-1));
%     hold on;
%     errorbar(xdat, ydat, psc_sem,'o');
%     xbounds = [min(xdat)-1, max(xdat)+1];
%     ybounds = [min(ydat)-0.25, max(ydat)+0.25];
%     xlim(xbounds);ylim(ybounds);
%     c = polyfit(xdat(1:2),ydat(1:2),1);
%     plot(linspace(xbounds(1),xbounds(2),20),linspace(xbounds(1),xbounds(2),20)*c(1)+c(2),'k--');
%     hold off;
%
% elseif optInputs(varargin, 'plotcontrasts') || optInputs(varargin, 'plotcontrasts_resolvability')
%
%     if optInputs(varargin, 'plotcontrasts_resolvability')
%         baselines = varargin{optInputs(varargin, 'plotcontrasts_resolvability')+1};
%         stat_baselines = nan(size(baselines));
%         psc_baselines = nan(size(baselines));
%         for j = 1:length(baselines)
%             x = strcmp(baselines{j}, contrasts);
%             psc_baselines(j) = mean(psc_contrasts(:,x));
%             stat_baselines(j) = mean(res_contrasts(:,x));
%         end
%
%         c = polyfit(stat_baselines,psc_baselines,1);
%         psc_predicted = c(2) + c(1)*res_contrasts;
%         psc_contrasts = psc_contrasts - repmat(psc_predicted, size(psc_contrasts,1), 1);
% %         figure;
% %         plot(res_contrasts,mean(psc_contrasts),'o');
% %         hold on;
% %         xbounds = [min(res_contrasts)-1, max(res_contrasts)+1];
% %         ybounds = [min(mean(psc_contrasts))-0., max(mean(psc_contrasts))+0.25];
% %         plot(linspace(xbounds(1),xbounds(2),20),linspace(xbounds(1),xbounds(2),20)*c(1)+c(2),'k--');
% %         xlim(xbounds); ylim(ybounds);
% %         hold off;
% %
% %         keyboard;
%
%     end
%
%     psc_mean = squeeze(mean(psc_contrasts));
%     %     if optInputs(varargin, 'sem_nonorm')
%     psc_sem = squeeze(std(psc_contrasts,[],1)/sqrt(size(psc_contrasts,1)-1));
%     %     else
%     %         [psc_sem psc_norm] = stderr_withsub(psc_contrasts);
%     %     end
%
%     if optInputs(varargin,'noplot')
%         return;
%     end
%
%     if optInputs(varargin, 'scatter')
%         if optInputs(varargin, 'sem_nonorm')
%             myscatter(psc_contrasts,contrasts,contrasts);
%         else
%             myscatter(mean(psc_contrasts(:)) + psc_norm, contrasts, plotorder)
%         end
%     else
%         mybar(psc_mean,contrasts,contrasts,'errorbar',psc_sem,varargin{:});
%     end
%
%     title(titlestring);
% else
%
%     psc_mean = squeeze(mean(psc));
%     if optInputs(varargin, 'sem_nonorm')
%         psc_sem = squeeze(std(psc,[],1)/sqrt(size(psc,1)-1));
%     else
%         [psc_sem psc_norm] = stderr_withsub(psc);
%     end
%
%     if optInputs(varargin,'noplot')
%         return;
%     end
%
%     plotorder = read_plotorder(usubs(1), test_runtype);
%     if optInputs(varargin,'sigav_timecourse')
%         if strcmp(read_expname,'DP')
%             x = 6;
%         else
%             x = 5;
%         end
%         myplot(win',psc_mean(:,1:x),conds(1:x),conds(1:x),'errorbar',psc_sem);
%         myplot(win',psc_mean(:,x+1:end),conds(x+1:end),conds(x+1:end),'errorbar',psc_sem);
%     elseif optInputs(varargin, 'scatter')
%         if optInputs(varargin, 'sem_nonorm')
%             myscatter(psc,conds,plotorder);
%         else
%             myscatter(mean(psc(:)) + psc_norm, conds, plotorder)
%         end
%     else
%         mybar(psc_mean,conds,plotorder,'errorbar',psc_sem,varargin{:});
%     end
%
%     title(titlestring);
%
% end


%% Scraps

% % conditions
% conds = read_conditions(s, test_runtype);
% [psc runs_used mainstr] = roi_psc(s, test_runtype, testmodel, masktype, roi, localizer, loc_runtype, locmodel, varargin{:});
% nruns = length(runs_used);
%
% if nruns ~= size(psc,1);
%     error('size of runs_used does not match dimensions of psc');
% end
%
% if sigav
%
%     % collapses across runs (which is the first dimension)
%     % [nruns ntps nconds] -> [ntp nconds]
%     psc_mean_time_conds = squeeze( nanmean(psc,1) );
%     psc_sem_time_conds = squeeze( nanstd(psc,0,1)/sqrt(nruns-1) );
%
%     % [nruns ntps nconds] -> [nruns ntp]
%     % [nruns ntp] -> ntp
%     psc_mean_run_time = nanmean(psc,3);
%     psc_mean_time = squeeze( nanmean(psc_mean_run_time,1) );
%     psc_sem_time = squeeze( nanstd(psc_mean_run_time,0,1) / sqrt(nruns-1) );
%
%     % [ntps nconds] -> nconds
%     psc_mean_peak_conds = psc_mean_time_conds(peak_smp,:);
%     psc_sem_peak_conds = psc_sem_time_conds(peak_smp,:);
%
%     % [nruns ntps nconds] -> [nruns nconds]
%     % [nruns nconds] -> nconds
%     psc_mean_plat_runs_conds = squeeze( nanmean(psc(:,plat_smp,:),2) );
%     psc_mean_plat_conds = nanmean(psc_mean_plat_runs_conds, 1);
%     psc_sem_plat_conds = nanstd(psc_mean_plat_runs_conds,0,1) / sqrt(nruns-1);
%
%     varargout = {psc_mean_plat_conds, mainstr, psc_mean_peak_conds, psc_mean_time_conds};
%
% else
%
%     % [nruns nconds] -> nconds
%     psc_mean_conds = nanmean(psc,1);
%     psc_sem_conds = nanstd(psc,0,1)/sqrt(nruns-1);
%
%     varargout = {psc_mean_conds, mainstr};
%
% end
%
%
% if ~plotfig; return; end;
%
% % this ttest procedure will not work as is
% % pvals = nan(1,length(ttestpairs));
% % for t = 1:length(ttestpairs);
% %
% %     % find condition indices
% %     c1 = strmatch(ttestpairs{t}{1},conds,'exact'); c2 = strmatch(ttestpairs{t}{2},conds,'exact');
% %
% %     % find shared conditions
% %     [zz c1inds c2inds] = intersect(testruns{c1},testruns{c2});
% %     [zz pvals(t) zz zz] = ttest(pe{c1}(c1inds),   pe{c2}(c2inds));
% %     fprintf('%s vs. %s: p = %f\n',ttestpairs{t}{1},ttestpairs{t}{2},pvals(t))
% % end
%
% % title string
% if strcmp(masktype,'anat');
%     titlestring = ['Atlas: ' strrep(roi,'_',' ') ', Sub: ' num2str(s)];
% else
%     a = regexp(localizer,'_vs_','split');
%     if length(a) == 1;
%         a = [a, 'NULL'];
%     end
%
%     b = cell(1,2);
%     for i = 1:length(a)
%
%         p = read_plotformat(a{i});
%         b{i} = p.sym;
%
%     end
%
%     if negcontrast
%         b = fliplr(b);
%     end
%
%     contraststring = [b{1} ' > ' b{2}];
%     titlestring = [read_expname ' s' num2str(s) ', Atlas: ' strrep(roi,'_',' ') ', Loc: ' contraststring ];% strrep(test_runtype, '_', ' ')  ', ' strrep(masktype,'_',' ') ', '
% end
%
% plotorder = read_plotorder(s, test_runtype);
% plotorder_inds = nan(length(conds),1);
%
% fnum = figure;
% hold on;
% names = cell(length(conds),1);
% for i = 1:length(plotorder)
%
%     p = read_plotformat(plotorder{i});
%     names{i} = p.sym;
%
%     x = strmatch(plotorder{i},conds,'exact');
%     plotorder_inds(i) = x;
%
%     if sigav
%
%         figure(fnum);
%         hold on;
%         bar(i, psc_mean_plat_conds(x), 'FaceColor', p.rgbcol/255);
%
%         if optInputs(varargin, 'allfigures')
%
%             figure(fnum+1);
%             hold on;
%             bar(i, psc_mean_peak_conds(x), 'FaceColor', p.rgbcol/255);
%
%             figure(fnum+2);
%             hold on;
%             plot(win, psc_mean_time_conds(:,x), p.line, 'Color', p.rgbcol/255, 'LineWidth', 2);
%
%             figure(fnum+3);
%             hold on;
%             plot(win, psc_mean_time, p.line, 'Color', p.rgbcol/255, 'LineWidth', 2);
%
%         end
%
%     else
%
%         bar(i, psc_mean_conds(x),'FaceColor', p.rgbcol/255);
%
%     end
%
% end
%
% if sigav
%
%
%     % Mean Response for 'Plateau'
%     figure(fnum);
%     hold on;
%     errorbar(1:length(plotorder), psc_mean_plat_conds(plotorder_inds), psc_sem_plat_conds(plotorder_inds),'k.');
%     title(['Mean Resp (' num2str(plat_time(1)) ' to ' num2str(plat_time(end)) ' sec), ' titlestring]);
%     x(1) = min(psc_mean_plat_conds - psc_sem_plat_conds) - 0.5;
%     x(1) = min(x(1),0);
%     x(2) = max(psc_mean_plat_conds + psc_sem_plat_conds) + 0.5;
%     ylim(x);
%
%     xlim([0, length(plotorder)+1]);
%     set(gca,'XTick',1:length(names));
%     set(gca,'XTickLabel',names,'FontWeight',fontweight,'FontSize',fontsize,'Position',[0.15 0.25 0.7 0.65]);
%     rotateticklabel(gca,45,0.75,fontsize,fontweight);
%     ylabel('% Signal Change');
%
%     if optInputs(varargin, 'allfigures')
%
%         % Peak Response
%         figure(fnum+1);
%         hold on;
%         errorbar(1:length(plotorder), psc_mean_peak_conds(plotorder_inds), psc_sem_peak_conds(plotorder_inds),'k.');
%         title(['Peak Resp (' num2str(peak_time) ' sec), ' titlestring]);
%         x(1) = min(psc_mean_peak_conds - psc_sem_peak_conds) - 0.5;
%         x(1) = min(x(1),0);
%         x(2) = max(psc_mean_peak_conds + psc_sem_peak_conds) + 0.5;
%         ylim(x);
%
%         xlim([0, length(plotorder)+1]);
%         set(gca,'XTick',1:length(names));
%         set(gca,'XTickLabel',names,'FontWeight',fontweight,'FontSize',fontsize,'Position',[0.15 0.25 0.7 0.65]);
%         rotateticklabel(gca,45,0.75,fontsize,fontweight);
%         ylabel('% Signal Change');
%
%         % Mean Timecourse Per Condition
%         figure(fnum+2);
%         hold on;
%         x(1) = min(psc_mean_time_conds(:) - psc_sem_time_conds(:)) - 0.5;
%         x(1) = min(x(1),0);
%         x(2) = max(psc_mean_time_conds(:) + psc_sem_time_conds(:)) + 0.5;
%         ylim(x);
%
%         xlabel('Time (s)');
%         legend(names,'Location','Best');
%         ylabel('% Signal Change');
%         title(['Timecourse, ' titlestring]);
%
%         % Mean Timecourse Averaged Across all Conditions
%         figure(fnum+3);
%         hold on;
%         errorbar(win, psc_mean_time, psc_sem_time,'k.');
%         title(['Mean Timecourse (All Conditions), ' titlestring]);
%         x(1) = min(psc_mean_time(:) - psc_sem_time(:)) - 0.5;
%         x(1) = min(x(1),0);
%         x(2) = max(psc_mean_time(:) + psc_sem_time(:)) + 0.5;
%         ylim(x);
%
%     end
%
%     if optInputs(varargin,'svg');
%         figure;
%         bar(1:length(plotorder),psc_mean_plat_conds(plotorder_inds));
%         hold on;
%         errorbar(1:length(conds),psc_mean_plat_conds(plotorder_inds),psc_sem_plat_conds(plotorder_inds),'k.');
%         hold off;
%         figfile = [figuredir  mainstr '.svg'];
%         if ~exist(figuredir, 'dir'); mkdir(figuredir); end
%         saveas(gcf,figfile,'png');
%         plot2svg_beta(figfile,gcf,'svg');
%     end
%
% else
%
%     errorbar(1:length(plotorder), psc_mean_conds(plotorder_inds), psc_sem_conds(plotorder_inds),'k.');
%     title(['Mean Response (Betas), ' titlestring]);
%     x(1) = min(psc_mean_conds(:) - psc_sem_conds(:)) - 0.5;
%     x(1) = min(x(1),0);
%     x(2) = max(psc_mean_conds(:) + psc_sem_conds(:)) + 0.5;
%     ylim(x);
%
%     xlim([0, length(plotorder)+1]);
%     set(gca,'XTick',1:length(names));
%     set(gca,'XTickLabel',names,'FontWeight',fontweight,'FontSize',fontsize,'Position',[0.15 0.25 0.7 0.65]);
%     rotateticklabel(gca,45,0.75,fontsize,fontweight);
%     ylabel('% Signal Change');
%
%     if optInputs(varargin,'svg');
%         figure;
%         bar(1:length(plotorder),psc_mean_conds(plotorder_inds));
%         hold on;
%         errorbar(1:length(conds),psc_mean_conds(plotorder_inds),psc_sem_conds(plotorder_inds),'k.');
%         hold off;
%         figfile = [figuredir  regexprep(mainstr,  '_s(\d)*',  ['_s' sprintf('%d',usubs)]) '.svg'];
%         if ~exist(figuredir, 'dir'); mkdir(figuredir); end
%         saveas(gcf,figfile,'png');
%         plot2svg_beta(figfile,gcf,'svg');
%     end
%
% end
%

% figfile = [figuredir mainstr];
% if ~exist(figdir, 'dir'); mkdir(figdir); end
% saveas(gcf,figfile,'png');


%%

% if optInputs(varargin,'resolvability') || optInputs(varargin,'periodicity')
%
%     if optInputs(varargin,'resolvability')
%         baselines = varargin{optInputs(varargin,'resolvability')+1};
%         load('peaks.mat');
%         stat = p.statmean{4};
%     elseif optInputs(varargin,'periodicity')
%         baselines = varargin{optInputs(varargin,'periodicity')+1};
%         load('praatstats.mat');
%         stat = p.statmean{1};
%     end
%
%     if ~isequal(p.conds, conds);
%         error('Conditions do not match.');
%     end
%
%     psc_baselines = nan(size(baselines));
%     stat_baselines = nan(size(baselines));
%     for j = 1:length(baselines)
%
%         evn = evgroup2evname(usubs(i),baselines{j},testmodel);
%         inds = nan(size(evn));
%         for k = 1:length(evn)
%             inds(k) = find(strcmp(evn{k},conds));
%         end
%
%         psc_baselines(j) = mean(mean(psc_sub(:,inds)));
%         stat_baselines(j) = mean(stat(inds));
%
%     end
%
%     c = polyfit(stat_baselines,psc_baselines,1);
%     psc_predicted = c(2) + c(1)*stat;
%     psc_sub = psc_sub - repmat(psc_predicted, size(psc_sub,1), 1);
%
% end

