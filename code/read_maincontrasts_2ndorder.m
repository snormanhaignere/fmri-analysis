function contrasts = read_maincontrasts_2ndorder(exp,runtype)

% contrasts = {'harm-lowfreq_vs_noise-lowfreq_vs2_harm-unres_vs_noise-highfreq', 'harm-highfreq_vs_noise-highfreq_vs2_harm-unres_vs_noise-highfreq',...
%     'harm-lowfreq-lowmask_vs_noise-lowfreq-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask', 'harm-highfreq-lowmask_vs_noise-highfreq-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'};

% x1 = {'harm-lowfreq_vs_noise-lowfreq_vs2_harm-unres_vs_noise-highfreq', 'harm-lowfreq-lowmask_vs_noise-lowfreq-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'};
%
% x2 = {'harm-highfreq_vs_noise-highfreq_vs2_harm-unres_vs_noise-highfreq','harm-highfreq-lowmask_vs_noise-highfreq-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'};
%
% x3 = {'harm_vs_noise_vs2_harm-unres_vs_noise-highfreq','harm-lowmask_vs_noise-lowmask_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'};
%
% x4 = {'harm-unres_vs_harm-unres-schroeder-neg_vs2_harm-unres-lowmask_vs_harm-unres-schroeder-neg-lowmask', 'harm-unres_vs_noise-highfreq_vs2_harm-unres-lowmask_vs_noise-highfreq-lowmask'};

% contrasts= [x1,x2,x3,x4];

switch exp
    case 'pitch_resthr_v4'
        switch runtype
            case 'main'
                contrasts = {...
                    'freq1200-2400_harm-sine6-12_vs_freq1200-2400_noise_vs2_freq2400-4800_harm-sine12-24_vs_freq2400-4800_noise',...
                    };
            otherwise
                contrasts = {};
        end
    otherwise
        contrasts = {};
end