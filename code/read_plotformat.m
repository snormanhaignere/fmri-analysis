function p = read_plotformat(cond)

switch cond
%     case 'orchestra'
%         p.col = 'turquoise7'; p.line = '-'; p.sym = 'OC'; p.width = 2; p.mark = 'v'; p.rgbcol = [2 85 85];
%     case 'solo'
%         p.col = 'turquoise5'; p.line = '-'; p.sym = 'SL'; p.width = 2; p.mark = 'o'; p.rgbcol = [63 136 136];
%     case 'noise-ll'
%         p.col = 'turquoise3'; p.line = '-'; p.sym = 'N'; p.width = 2; p.mark = 'x'; p.rgbcol = [226 94 6];
%     case 'noise-ll-jumps'
%         p.col = 'turquoise1'; p.line = '-'; p.sym = 'jN'; p.width = 2; p.mark = '+'; p.rgbcol = [255 173 119];
%     case 'midi-intact'
%         p.col = 'darkgreen'; p.line = '-'; p.sym = 'IM'; p.width = 2; p.mark = 'o';
%     case 'midi-scram';
%         p.col = 'springgreen'; p.line = '-'; p.sym = 'SM'; p.width = 2; p.mark = 'o';
%     case 'harm-hh'
%         p.col = 'green'; p.line = '-'; p.sym = 'HH'; p.width = 2; p.mark = '*'; p.rgbcol = [18 54 36];
%     case 'harm-hl';
%         p.col = 'springgreen'; p.line = '-'; p.sym = 'HL'; p.width = 2; p.mark = 's'; p.rgbcol = [42 98 70];
%     case 'harm-ll';
%         p.col = 'teal-green'; p.line = '-'; p.sym = 'LL'; p.width = 2; p.mark = 'd'; p.rgbcol = [59 135 97];
%     case 'harm-ll-jitter';
%         p.col = 'darkgreen'; p.line = '-'; p.sym = 'jLL'; p.width = 2; p.mark = 'v'; p.rgbcol = [175 238 175];
%     case 'harm-ll-fast'
%         p.col = 'blue5'; p.line = ':'; p.sym = 'fLL'; p.width = 2; p.mark = '^'; p.rgbcol = [30 70 191];
%     case 'harm-ll-slow'
%         p.col = 'blue3'; p.line = ':'; p.sym = 'sLL'; p.width = 2; p.mark = '<'; p.rgbcol = [98 160 254];
        
    case 'NULL'
        p.col = 'green'; p.line = '--'; p.sym = 'Null'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];

    case 'all'
        p.col = 'green'; p.line = '--'; p.sym = 'All'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];

    case 'harm-res-unres'
        p.col = 'green'; p.line = '-'; p.sym = 'Pi Res+Unres'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36]
    case 'harm-res-unres-lowmask'
        p.col = 'green'; p.line = '-'; p.sym = 'Pi Res+Unres Mask'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
        
    case 'harm-boost'
        p.col = 'green'; p.line = '--'; p.sym = 'Res - Noise'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-unres-boost'
        p.col = 'green'; p.line = '--'; p.sym = 'Unres - Noise'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-lowmask-boost'
        p.col = 'green'; p.line = '--'; p.sym = 'Res - Noise (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-unres-lowmask-boost'
        p.col = 'green'; p.line = '--'; p.sym = 'Unres - Noise (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-boost_vs_harm-unres-boost'
        p.col = 'green'; p.line = '--'; p.sym = 'Low Res - Noise > Unres - Noise'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-lowfreq_vs_noise-lowfreq_vs2_harm-unres_vs_noise-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Low Res - Noise > Unres - Noise'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-highfreq_vs_noise-highfreq_vs2_harm-unres_vs_noise-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Hi Res - Noise > Unres - Noise'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-lowfreq-lowmask_vs_noise-lowfreq-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Low Res - Noise > Unres - Noise (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-highfreq-lowmask_vs_noise-highfreq-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Hi Res - Noise > Unres - Noise (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm_vs_noise_vs2_harm-unres_vs_noise-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Res - Noise > Unres - Noise'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-lowmask_vs_noise-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Res - Noise > Unres - Noise (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];

    case 'harm-highfreq_vs_noise-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi HiFr > No HiFr'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-highfreq-schroeder-negN-1_vs_noise-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Schr Pi HiFr > No HiFr'; p.width = 2; p.mark = '*'; p.rgbcol = max([111 28 28]-40,0);%[18 54 36];
    case 'harm-lowfreq_vs_noise-lowfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi LoFr > No LoFr'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-unres_vs_noise-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Unres > No HiFr'; p.width = 2; p.mark = '*'; p.rgbcol = [241 198 198];%[18 54 36];
    case 'harm-unres-schroeder-neg_vs_noise-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Schr Unres > No HiFr'; p.width = 2; p.mark = '*'; p.rgbcol = [241 198 198]-40;%[18 54 36];
    case 'harm-unres-lowmask_vs_harm-unres-schroeder-neg-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Unres > Unres Schr Neg (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];
    case 'harm-unres_vs_harm-unres-schroeder-neg'
        p.col = 'green'; p.line = '--'; p.sym = 'Unres > Unres Schr Neg'; p.width = 2; p.mark = '*'; p.rgbcol = [150 150 150];%[18 54 36];

    case 'harm-highfreq-lowmask_vs_noise-highfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi HiFr > No HiFr (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-lowfreq-lowmask_vs_noise-lowfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi LoFr > No LoFr (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-unres-lowmask_vs_noise-highfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Unres > No HiFr (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [241 198 198];%[18 54 36];
    case 'harm-unres-schroeder-neg-lowmask_vs_noise-highfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Schr Unres > No HiFr (Mask)'; p.width = 2; p.mark = '*'; p.rgbcol = [241 198 198]-40;%[18 54 36];
        
    case 'harm-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi HiFr+LoFr Mask'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'noise-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'No HiFr+LoFr Mask'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];

    case 'harm'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi HiFr+LoFr'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-wawomask'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi All'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-highfreq-wawomask'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi HiFr WaWo Mask'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'noise'
        p.col = 'green'; p.line = '--'; p.sym = 'Noise'; p.width = 2; p.mark = '*'; p.rgbcol = [61 76 158];%[18 54 36];
    case 'noise-wawomask'
        p.col = 'green'; p.line = '--'; p.sym = 'No All'; p.width = 2; p.mark = '*'; p.rgbcol = [61 76 158];%[18 54 36];
    case 'noise-highfreq-wawomask'
        p.col = 'green'; p.line = '--'; p.sym = 'No HiFr WaWo Mask'; p.width = 2; p.mark = '*'; p.rgbcol = [61 76 158];%[18 54 36];
        
    case 'highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'High Freq'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'lowfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Low Freq'; p.width = 2; p.mark = '*'; p.rgbcol = [61 76 158];%[18 54 36];
        
    case 'harm-highfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi HiFr'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-highfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi HiFr Low Mask'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-unres'
        p.col = 'blue5'; p.line = '--'; p.sym = 'Pi Unres'; p.width = 2; p.mark = '^'; p.rgbcol = [241 198 198];%[30 70 191];
    case 'harm-unres-lowmask'
        p.col = 'blue5'; p.line = '--'; p.sym = 'Pi Unres Low Mask'; p.width = 2; p.mark = '^'; p.rgbcol = [241 198 198];%[30 70 191];
    case 'noise-highfreq'
        p.col = 'turquoise3'; p.line = '--'; p.sym = 'Noise HiFr'; p.width = 2; p.mark = 'x'; p.rgbcol = [61 76 158];%[191 69 29];%[226 94 6];
    case 'noise-highfreq-lowmask'
        p.col = 'turquoise3'; p.line = '--'; p.sym = 'Noise HiFr Low Mask'; p.width = 2; p.mark = 'x'; p.rgbcol = [61 76 158];%[191 69 29];%[226 94 6];
    case 'harm-lowfreq'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi LoFr'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'harm-lowfreq-lowmask'
        p.col = 'green'; p.line = '--'; p.sym = 'Pi LoFr Low Mask'; p.width = 2; p.mark = '*'; p.rgbcol = [111 28 28];%[18 54 36];
    case 'noise-lowfreq'
        p.col = 'turquoise3'; p.line = '--'; p.sym = 'Noise LoFr'; p.width = 2; p.mark = 'x'; p.rgbcol = [61 76 158];%[191 69 29];%[226 94 6];
    case 'noise-lowfreq-lowmask'
        p.col = 'turquoise3'; p.line = '--'; p.sym = 'Noise LoFr Low Mask'; p.width = 2; p.mark = 'x'; p.rgbcol = [61 76 158];%[191 69 29];%[226 94 6];
    case 'harm-unres-schroeder-neg'
        p.col = 'blue5'; p.line = '--'; p.sym = 'Pi Unres Schr'; p.width = 2; p.mark = '^'; p.rgbcol = [241 198 198]-40;%[30 70 191];
    case 'harm-unres-schroeder-neg-lowmask'
        p.col = 'blue5'; p.line = '--'; p.sym = 'Pi Unres Schr Low Mask'; p.width = 2; p.mark = '^'; p.rgbcol = [241 198 198]-40;%[30 70 191];
    case 'harm-highfreq-schroeder-negN-1'
        p.col = 'blue5'; p.line = '--'; p.sym = 'Pi HiFr Schr'; p.width = 2; p.mark = '^'; p.rgbcol = max([111 28 28]-40,0);%[30 70 191];
        
    case 'speech-voiced';
        p.col = 'springgreen'; p.line = '.-'; p.sym = 'Voiced Speech'; p.width = 2; p.mark = 's'; p.rgbcol = [111 28 28];%[199 73 73];%;[154 49 49];%[42 98 70];
    case 'speech-normal';
        p.col = 'teal-green'; p.line = '.-'; p.sym = 'Normal Speech'; p.width = 2; p.mark = 'd'; p.rgbcol = [241 198 198];%[224 104 104];%[59 135 97];
    case 'speech-whisper'
        p.col = 'black'; p.line = '.-'; p.sym = 'Whisper Speech'; p.width = 2; p.mark = '*'; p.rgbcol = [61 76 158];%[207 139 71];%;[116 72 135];%[66 66 66];
    case 'speech-whisper-frags'
        p.col = 'black'; p.line = '.-'; p.sym = 'Whisper Fragments'; p.width = 2; p.mark = '*'; p.rgbcol = [61 76 158];%[207 139 71];%;[116 72 135];%[66 66 66];
        
    case 'sentences';
        p.col = 'darkgreen'; p.line = '-'; p.sym = 'Sentences'; p.width = 2; p.mark = 'v'; p.rgbcol = [0 0 0];%[175 238 175];
        
    case 'nonwords';
        p.col = 'darkgreen'; p.line = '-'; p.sym = 'Nonwords'; p.width = 2; p.mark = 'v'; p.rgbcol = [0 0 0];%[175 238 175];
        
    case 'solo';
        p.col = 'darkgreen'; p.line = '-'; p.sym = 'Solo Music'; p.width = 2; p.mark = 'v'; p.rgbcol = [111 28 28];%[175 238 175];
    case 'orchestra';
        p.col = 'darkgreen'; p.line = '-'; p.sym = 'Orchestra Music'; p.width = 2; p.mark = 'v'; p.rgbcol = [111 28 28];%[175 238 175];
    case 'drums'
        p.col = 'turquoise1'; p.line = '-'; p.sym = 'Drums'; p.width = 2; p.mark = '+'; p.rgbcol = [61 76 158];%[255 173 119];
        
    case 'env-pitch'
        p.col = 'crimson'; p.line = ':'; p.sym = 'Env Sounds'; p.width = 2; p.mark = '>'; p.rgbcol = [111 28 28];
    case 'animals-pitch'
        p.col = 'crimson'; p.line = ':'; p.sym = 'Ani Calls'; p.width = 2; p.mark = '>'; p.rgbcol = [111 28 28];
    case 'env-nopitch'
        p.col = 'red'; p.line = ':'; p.sym = 'NonPi Env'; p.width = 2; p.mark = 'p'; p.rgbcol = [61 76 158];
        
    case 'harm-f0same-freqsame'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm F0Sa FrSa'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
    case 'harm-f0same-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm F0Sa FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
    case 'harm-f0diff-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm F0Di FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
        
    case 'harm-low-f0same-freqsame'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm Low F0Sa FrSa'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
    case 'harm-low-f0same-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm Low F0Sa FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
    case 'harm-low-f0diff-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm Low F0Di FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
        
    case 'harm-high-f0same-freqsame'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm High F0Sa FrSa'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
    case 'harm-high-f0same-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm High F0Sa FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
    case 'harm-high-f0diff-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Harm High F0Di FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [111 28 28];
        
    case 'inharm-f0same-freqsame'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm F0Sa FrSa'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
    case 'inharm-f0same-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm F0Sa FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
    case 'inharm-f0diff-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm F0Di FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
        
    case 'inharm-low-f0same-freqsame'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm Low F0Sa FrSa'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
    case 'inharm-low-f0same-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm Low F0Sa FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
    case 'inharm-low-f0diff-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm Low F0Di FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
        
    case 'inharm-high-f0same-freqsame'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm High F0Sa FrSa'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
    case 'inharm-high-f0same-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm High F0Sa FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
    case 'inharm-high-f0diff-freqdiff'
        p.col = 'red'; p.line = ':'; p.sym = 'Inharm High F0Di FrDi'; p.width = 2; p.mark = 'p'; p.rgbcol = [175 238 175];
        
    otherwise
        
        p.col = 'red'; p.line = ':'; p.sym = strrep(cond,'_',' '); p.width = 2; p.mark = 'p'; p.rgbcol = [125 125 125];
        
end

% if findstr(cond,'harm-')
%     p.rgbcol = [255 120 120];
% end
% if findstr(cond,'inharm-')
%     p.rgbcol = [241 198 198];
% end
% 
% if findstr(cond,'noise-')
%     p.rgbcol = [61 76 158];
% end
% 
% if findstr(cond,'1200-2400')
%     p.rgbcol = [120 255 120];
% end
% 
% if findstr(cond,'2400-4800')
%     p.rgbcol = [120 120 255];
% end

if strfind(cond,'intact')
    p.rgbcol = [61 76 158];
end

if strfind(cond,'scram')
    p.rgbcol = [120 120 255];
end

if strfind(cond,'quilt')
    p.rgbcol = [120 120 255];
end

if strfind(cond,'german-intact')
    p.rgbcol = [111 28 28];
end

if strfind(cond,'german-quilt')
    p.rgbcol = [111 70 70];
end

if strfind(cond,'env-pitch')
    p.rgbcol = [125 125 125];
end


        %     case 'animals'
%         p.col = 'violetred4'; p.line = '--';p.sym = 'AN'; p.width = 2; p.mark = 'o';
%     case 'events'
%         p.col = 'violetred3'; p.line = '--'; p.sym = 'EV'; p.width = 2; p.mark = 'o';
%     case 'textures'
%         p.col = 'violetred2'; p.line = '--'; p.sym = 'TX'; p.width = 2; p.mark = 'o';
%     case 'irn-ll'
%         p.col = 'black'; p.line = '--'; p.sym = 'IRN'; p.width = 2; p.mark = '*'; p.rgbcol = [66 66 66];
%     case 'nonwords'
%         p.col = 'gray'; p.line = '--'; p.sym = 'NW'; p.width = 2; p.mark = 's'; p.rgbcol = [220 220 220];
%     case 'pitched'
%         p.col = 'red'; p.line = '-'; p.sym = 'PI'; p.width = 2; p.mark = 'o';
%     case 'nonpitched'
%         p.col = 'red'; p.line = '-'; p.sym = 'nPI'; p.width = 2; p.mark = 'o';
%     case 'sentences_pitchfrags';
%         p.col = 'black'; p.line = '-'; p.sym = 'SF'; p.width = 2; p.mark = 'o';
%     case 'nonwords_pitchfrags'
%         p.col = 'gray'; p.line = '-'; p.sym = 'NF'; p.width = 2; p.mark = 'o';
%     case 'speakers_ED'
%         p.col = 'orange'; p.line = '-'; p.sym = 'ED'; p.width = 2; p.mark = 'o';
%     case 'speakers_SNH'
%         p.col = 'orange'; p.line = '-'; p.sym = 'SNH'; p.width = 2; p.mark = 'o';
%     case 'speakers_BD'
%         p.col = 'orange'; p.line = '-'; p.sym = 'BD'; p.width = 2; p.mark = 'o';
%     case 'speakers_DD'
%         p.col = 'orange'; p.line = '-'; p.sym = 'DD'; p.width = 2; p.mark = 'o';

%     case 'animals-pitch'
%         p.col = 'violetred4'; p.line = '--';p.sym = 'pAN'; p.width = 2; p.mark = 'h'; p.rgbcol = [122 3 3];
