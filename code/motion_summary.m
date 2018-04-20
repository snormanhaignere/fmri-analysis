function [mean_rms, n_timepoints_above_0p5mm, n_timepoints_above_1mm, max_total_rms_bysubject] = motion_summary(exp, usubs, runtype, varargin)

addpath(genpath([pwd '/export_fig']));
figure_directory = [params('rootdir') exp '/figures/motion/'];
if ~exist(figure_directory,'dir');
    mkdir(figure_directory);
end

mean_rms_bysubject = nan(length(usubs),1);
n_timepoints_above_0p2mm_bysubject = nan(length(usubs),1);
n_timepoints_above_0p5mm_bysubject = nan(length(usubs),1);
n_timepoints_above_1mm_bysubject = nan(length(usubs),1);
rms_all = cell(length(usubs),1);

max_abs_rms_bysubject = nan(length(usubs),1);
abs_rms_all = cell(length(usubs),1);

for i = 1:length(usubs)
    
    [runnum, ~, scans] = read_runs(exp,usubs(i),runtype,varargin{:});
    %     if optInputs(varargin,'runnum')
    %         runnum = varargin{optInputs(varargin,'runnum')+1};
    %     end
    
    mean_rms_byrun = nan(length(runnum),1);
    n_timepoints_above_0p2mm_byrun = nan(length(runnum),1);
    n_timepoints_above_0p5mm_byrun = nan(length(runnum),1);
    n_timepoints_above_1mm_byrun = nan(length(runnum),1);
    rms_all{i} = [];
    
    max_abs_rms_byrun = nan(length(runnum),1);
    
    for k = 1:length(runnum)
        
        datadir = [params('rootdir') exp '/data/'];
        preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(usubs(i)) '/' runtype '_r' num2str(runnum(k)) '/'];
        
        try
            fid = fopen([preprocdir 'motcorr_rel.rms']);
            x = textscan(fid,'%f'); fclose(fid);
        catch
            keyboard
        end
        
        if optInputs(varargin, 'monkey')
            relrms = x{1}(:)/2.8571;
        else
            relrms = x{1}(:);
        end
        
        fid = fopen([preprocdir 'motcorr_abs.rms']);
        x = textscan(fid,'%f'); fclose(fid);
        
        if optInputs(varargin, 'monkey')
            absrms = x{1}(:)/2.8571;
        else
            absrms = x{1}(:);
        end
        
        mean_rms_byrun(k) = mean(relrms);
        n_timepoints_above_0p2mm_byrun(k) = mean(relrms(:)>0.2);
        n_timepoints_above_0p5mm_byrun(k) = mean(relrms(:)>0.5);
        n_timepoints_above_1mm_byrun(k) = mean(relrms(:)>1);
        rms_all{i} = [rms_all{i}; relrms];
        
        abs_rms_all{i} = [abs_rms_all{i}; absrms];
        max_abs_rms_byrun(k) = max(abs_rms_all{i});
        
    end
    max_abs_rms_bysubject(i) = max(max_abs_rms_byrun);
    mean_rms_bysubject(i) = mean(rms_all{i});
    n_timepoints_above_0p2mm_bysubject(i) = mean(rms_all{i}(:)>0.2);
    n_timepoints_above_0p5mm_bysubject(i) = mean(rms_all{i}(:)>0.5);
    n_timepoints_above_1mm_bysubject(i) = mean(rms_all{i}(:)>1);
    
end

if length(usubs)>1
    x_label = 'Subjects';
    x_axis = usubs;
    groups = ones(size(usubs));
    mean_rms = mean_rms_bysubject;
    n_timepoints_above_0p2mm = n_timepoints_above_0p2mm_bysubject;
    n_timepoints_above_0p5mm = n_timepoints_above_0p5mm_bysubject;
    n_timepoints_above_1mm = n_timepoints_above_1mm_bysubject;
else
    x_label = 'Runs';
    x_axis = runnum;
    groups = scans;
    mean_rms = mean_rms_byrun;
    n_timepoints_above_0p2mm = n_timepoints_above_0p2mm_byrun;
    n_timepoints_above_0p5mm = n_timepoints_above_0p5mm_byrun;
    n_timepoints_above_1mm = n_timepoints_above_1mm_byrun;
end

if optInputs(varargin, 'noplot')
    return;
end

for i = 1:length(usubs)
    figure;
    plot(rms_all{i});
    ylim([0 2]);
    xlim([0,length(rms_all{i})+1]);
    title(sprintf('%d, %s',usubs(i),strrep(runtype,'_',' ')));
    xlabel('scans');
    ylabel('RMS motion (mm)');
    export_fig([figure_directory 'motion_us-' num2str(usubs(i)) '_' runtype '_all_tps.pdf'],'-pdf','-nocrop');
end

figure;
set(gcf,'Position',[0 0 max(70*length(x_axis),600) 600]);
subplot(1,2,1);
set(gca,'Position',[0.05 0.2 0.4 0.6]);
colors = colormap('Lines(10)');
unique_groups = unique(groups);
for i = 1:length(unique_groups)
    xi = groups == unique_groups(i);
    plot(find(xi), mean_rms(xi),'o-','LineWidth',2, 'Color', colors(i,:));
    if i == 1
        hold on;
    end
end
ylim([0 0.5]);
xlim([0,length(mean_rms)+1]);
set(gca, 'XTick', 1:length(x_axis), 'XTickLabel', x_axis, 'FontSize',7);
ylabel('Mean RMS');
xlabel(x_label);

subplot(1,2,2);
set(gca,'Position',[0.55 0.2 0.4 0.6]);
% nTR = read_scanparams(exp,usubs(1),runtype,'nTR','run',1,varargin{:});
unique_groups = unique(groups);
for i = 1:length(unique_groups)
    xi = groups == unique_groups(i);
    %     plot(find(xi), n_timepoints_above_0p2mm(xi),'o:','LineWidth',2, 'Color', colors(i,:));
    if i == 1
        hold on;
    end
    plot(find(xi), n_timepoints_above_0p5mm(xi),'o--','LineWidth',2, 'Color', colors(i,:));
    plot(find(xi), n_timepoints_above_1mm(xi),'o-','LineWidth',2, 'Color', colors(i,:));
end
if length(usubs) > 1
    ylim([0 0.2]);
else
    ylim([0 0.5]);
end

xlim([0,length(mean_rms)+1]);
set(gca, 'XTick', 1:length(x_axis), 'XTickLabel', x_axis, 'FontSize',7);
ylabel('Proportion TPs > 0.5 or 1 mm');
xlabel(x_label);
legend({'TPs > 0.5 mm', 'TPs > 1 mm'});
export_fig([figure_directory 'motion_us' sprintf('-%d',usubs) '_' runtype '.pdf'],'-pdf','-nocrop');

if optInputs(varargin, 'plot_abs_motion')
    for i = 1:length(usubs)
        figure;
        plot(abs_rms_all{i});
        ylim([0 10]);
        xlim([0,length(abs_rms_all{i})+1]);
        title(sprintf('%d, %s',usubs(i),strrep(runtype,'_',' ')));
        xlabel('scans');
        ylabel('Abs RMS motion (mm)');
        export_fig([figure_directory 'motion_us-' num2str(usubs(i)) '_' runtype '_abs_rms_all_tps.pdf'],'-pdf','-nocrop');
    end
    
    if length(usubs) > 1
        figure;
        plot(max_abs_rms_bysubject,'o-','LineWidth',2);
        set(gca, 'XTick', 1:length(x_axis), 'XTickLabel', usubs, 'FontSize',10);
        ylim([0 10]);
        xlim([0,length(usubs)+1]);
        xlabel('Subjects');
        ylabel('Max Abs RMS motion (mm)');
        export_fig([figure_directory 'motion_us' sprintf('-%d',usubs) '_' runtype '_max_abs_rms.pdf'],'-pdf','-nocrop');
    end
end
% figure;
% plot(relrms);
% ylim([0 2]);
% title(sprintf('%s, run %d',strrep(runtype,'_',' '),r));

% fprintf('\n\n%s, run %d\n',runtype,r);
% fprintf('Mean: %.3f\n',mean(relrms(:)));
% fprintf('Time points > 0.5 mm: %d\n',);
% fprintf('Time points > 1 mm: %d\n',);