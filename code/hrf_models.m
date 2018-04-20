function irf = hrf_models(TR,varargin)
% function irf = hrf_models(TR,varargin)
% 
% plots hrf functions of different lags

%%
TR = 1
u = linspace(1.2,0.7,5);
a1 = 3;
b1 = 1.75;
d1 = 1.5;
a2 = 3;
b2 = 1.5;
d2 = 6.5;
s2 = 0.2;

colors{1} = [0 24 200]/255;
colors{2} = [0 70 210]/255;
colors{3} = [0 130 244]/255;
colors{4} = [35 200 255]/255;
colors{5} = [60 255 255]/255;
colors{6} = [120 244 244]/255;
colors{7} = [180 255 255]/255;
colors{8} = [255 255 255]/255;
colors = fliplr(colors);

warning('OFF');
tt = 0:TR:20;
irf = zeros(length(tt),length(u));
for uind = 1:length(u);
    temp1 = hrf_double_gamma(tt*u(uind),a1,b1,d1,a2,b2,d2,s2);
    irf(:,uind) = (1/max(temp1)).*temp1;
end
warning('ON');

if optInputs(varargin,'noplot'); return; end

figure;
hold on;
for uind = 1:length(u);
    plot(tt,squeeze(irf(:,uind))','Color',colors{mod(uind+1,length(colors))+1});
end
hold off