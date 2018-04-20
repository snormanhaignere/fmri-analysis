function [y, t] = hrf_fsfast_gamma(TR, varargin)



%% HRF
d = [2.25, 0, 0, 0];
tau = [1.25, 8, 8, 1.25];
alpha = [2, 0.3, 0.01, 1];
hrf_sign = [1, -1, -1, 1];
names = {'BOLD','MION','MION_CUSTOM1','BOLD_CUSTOM1'};

t = 0:TR:500;
y = zeros(length(t), length(d));
for i = 1:length(d)
    t_shift = (t-d(i));
    y(t_shift>0,i) = ((t_shift(t_shift>0)/tau(i)).^alpha(i) .* exp(-t_shift(t_shift>0)/tau(i)));
    y(:,i) = hrf_sign(i)*y(:,i)/max(y(:,i));
end

for i = 1:length(names)
    if optInputs(varargin, names{i})
        xi = strcmp(names{i}, names);
        y = y(:,xi);
    end
end

if optInputs(varargin,'noplot');
    return;
end

%% Plot HRF

clear colors;
colors{1} = [0 24 200]/255;
% colors{2} = [0 70 210]/255;
colors{2} = [0 130 244]/255;
% colors{4} = [35 200 255]/255;
% colors{3} = [60 255 255]/255;
% % colors{6} = [120 244 244]/255;
% colors{4} = [180 255 255]/255;
% % colors{8} = [255 255 255]/255;
colors = fliplr(colors);
figure;
for i = 1:length(d);
    plot(t,y(:,i)','Color',colors{mod(i-1,length(colors))+1}, 'LineWidth',2);
    if i == 1
        hold on;
    end
end
legend(names);
ylabel('Response');
xlabel('Time (s)');
export_fig('hrf_gamma.pdf','-pdf','-nocrop');

%% Plot block convolved with hrf

figure;
for j = 1:length(size(y))
    sr = 1/TR;
    blockdur = 50;
    % [onsets seq durs,zz,names] = textread(parafile, '%f%d%f%f%s'); %#ok<ASGLU>
    onsets = [0,blockdur];
    durs = [blockdur,blockdur];
    seq = [1,0];
    h = y(:,j);
    
    % upsample sequence vector
    totaltime = onsets(end)+durs(end)+10;
    nsmps = totaltime*sr + 1;
    seq_upsampled = zeros(nsmps,1);
    for i = 1:length(seq)
        if seq(i)~=0;
            inds = onsets(i)*sr + (1:(durs(i)*sr));
            seq_upsampled(inds) = seq(i);
        end
    end
    
    t = (0:nsmps-1)*TR;
    model = zeros(nsmps,max(seq));
    for i = 1:size(model,2)
        boxcar = zeros(nsmps,1);
        boxcar(seq_upsampled==i) = 1;
        
        x = conv(boxcar,h);
        model(:,i) = x(1:nsmps)/sum(h);
    end
    
    plot(t,model,'Color',colors{mod(j-1,length(colors))+1}, 'LineWidth',2);
    if j == 1
        hold on;
    end
end

legend(names);
ylabel('Response');
xlabel('Time (s)');
export_fig('hrf_gamma_block.pdf','-pdf','-nocrop');


%% FFT

t = 0:TR:1e3;
y = zeros(length(t), length(d));
for i = 1:length(d)
    t_shift = (t-d(i));
    y(t_shift>0,i) = ((t_shift(t_shift>0)/tau(i)).^alpha(i) .* exp(-t_shift(t_shift>0)/tau(i)));
    y(:,i) = y(:,i)/max(y(:,i));
end

figure;
for i = 1:size(y,2)
    [px,f] = fftplot2(y(:,i), 1/TR, 'noplot');
    semilogx(f,10*log10(px),'Color',colors{mod(i-1,length(colors))+1},'LineWidth',2);
    ylim([-100 -20]);
    if i == 1
        hold on;
    end
end
legend(names);
ylabel('Power');
xlabel('Frequency (Hz)');
export_fig('hrf_spectra.pdf','-pdf','-nocrop');




