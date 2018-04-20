function mybar(data,conds,plotorder,varargin)

fontsize = 14;
fontweight = 'Demi';

names = cell(length(plotorder),1);
inds = nan(length(plotorder),1);

if optInputs(varargin, 'nofigure')
    fh = gcf;
else
    fh = figure;
end
hold on;
for i = 1:length(plotorder)

    p = read_plotformat(plotorder{i});
    
    x = find(strcmp(plotorder{i},conds));
    if length(x) > 1
        warning('MATLAB:CUSTOM','Multiple condition matches in calling my bar.');
    end
    inds(i) = x(1);
    bar(i, data(inds(i)),'FaceColor', p.rgbcol/255);
    names{i} = p.sym; 
    
end

if optInputs(varargin,'errorbar')
    sem = varargin{optInputs(varargin,'errorbar')+1};
    errorbar(1:length(plotorder),data(inds),sem(inds),'k.');
end

if optInputs(varargin,'errorbar2')
    lsem = varargin{optInputs(varargin,'errorbar2')+1};
    usem = varargin{optInputs(varargin,'errorbar2')+2};
    errorbar(1:length(plotorder),data(inds),data(inds)-lsem(inds),usem(inds)-data(inds),'k.');
end

if optInputs(varargin, 'simplot')
    clf(fh);
    bar(1:length(plotorder),data(inds));
    hold on;
    if optInputs(varargin,'errorbar')
        sem = varargin{optInputs(varargin,'errorbar')+1};
        errorbar(1:length(plotorder),data(inds),sem(inds),'k.');
    end
    if optInputs(varargin,'errorbar2')
        lsem = varargin{optInputs(varargin,'errorbar2')+1};
        usem = varargin{optInputs(varargin,'errorbar2')+2};
        errorbar(1:length(plotorder),data(inds),lsem(inds),usem(inds),'k.');
    end
    return;
end

if optInputs(varargin,'ylabel')
    ylabel(varargin{optInputs(varargin,'ylabel')+1});
end

if optInputs(varargin,'title')
    title(varargin{optInputs(varargin,'title')+1});
end

if optInputs(varargin,'ylim')
    ylim(varargin{optInputs(varargin,'ylim')+1});
end

xlim([0, length(plotorder)+1]);
set(gca,'XTick',1:length(plotorder));
set(gca,'XTickLabel',names,'FontWeight',fontweight,'FontSize',fontsize,'Position',[0.15 0.35 0.7 0.4]);

tilt = 45;
if optInputs(varargin,'tilt')
    tilt = varargin{optInputs(varargin,'tilt')+1};
end
rotateticklabel(gca,tilt,0.75,fontsize,fontweight);
hold off;