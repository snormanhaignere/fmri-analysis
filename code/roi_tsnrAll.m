function [tsnr] = roi_tsnrAll(exp, usubs, roi, runtype, type, varargin)

addpath(genpath([pwd '/export_fig']));
figure_directory = [params('rootdir') exp '/figures/tsnr/'];
if ~exist(figure_directory,'dir');
    mkdir(figure_directory);
end

tsnr_subjects = nan(size(usubs));
for i = 1:length(usubs)
    [tsnr_runs,runs,scans] = roi_tsnr_v2(exp, usubs(i), roi, runtype, type, varargin{:});
    tsnr_subjects(i) = mean(tsnr_runs);
end

if length(usubs)>1
    x_label = 'Subjects';
    x_axis = usubs;
    groups = ones(size(usubs));
    tsnr = tsnr_subjects;
else
    x_label = 'Runs';
    x_axis = runs;
    groups = scans;
    tsnr = tsnr_runs;
end

figure;
colors = colormap('Lines(10)');
set(gcf,'Position',[0 0 max(30*length(x_axis),600) 600]);
unique_groups = unique(groups);
for i = 1:length(unique_groups)
    xi = groups == unique_groups(i);
    try
    plot(find(xi), tsnr(xi),'o-','LineWidth',2, 'Color', colors(i,:));
    catch
        keyboard;
    end
    if i == 1
        hold on;
    end
end
if median(tsnr) < 30
    ylim([0 30]);
else
    ylim([0 150]);
end
xlim([0,length(x_axis)+1]);
set(gca, 'XTick', 1:length(x_axis), 'XTickLabel', x_axis);
ylabel('tSNR');
xlabel(x_label);
export_fig([figure_directory 'tsnr_us' sprintf('-%d',usubs) '_' roi '_' runtype '_' type '.pdf'],'-pdf','-nocrop');
hold off;