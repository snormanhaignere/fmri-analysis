load('canonical_pitch_block_response.mat');



[hg, tg] = hrf_fsfast_gamma(0.2,'BOLD_CUSTOM1','noplot');
xi = tg<max(t);
tg = tg(xi);
hg = hg(xi);
yu = interp1(t,y,tg,'cubic')';
yu = yu/max(yu(:));

boxcar = zeros(length(tg),1);
boxcar(tg<17) = 1;
yg = conv(boxcar, hg);
yg = yg(1:length(tg));
yg = yg/max(yg(:));
plot(tg',[yu,hg,yg])
