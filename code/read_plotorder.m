function plotorder = read_plotorder(exp, us, runtype, varargin)

% plotorder = read_plotorder(s, runtype)
%
% returns a list of conditions in the order to be plotted

switch exp
    
    case 'mspec'
        plotorder = {'sentences','nonwords','orchestra','solo','animals-pitch','env-pitch'};
        
    case 'pitch_resthr_v3'
        plotorder = {...
            'freq1200-2400_harm-sine3-6',     'freq1200-2400_harm-sine4-8',       'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12', ...
            'freq1200-2400_harm-sine8-16',    'freq1200-2400_harm-sine10-20',     'freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30',...
            'freq1200-2400_noise',...
            ...
            'freq2400-4800_harm-sine3-6',     'freq2400-4800_harm-sine4-8',       'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12',...
            'freq2400-4800_harm-sine8-16',    'freq2400-4800_harm-sine10-20',     'freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30',...
            'freq2400-4800_noise',...
            ...
            'freq1200-2400_harm-sine3-6','freq1200-2400_harm-schr-1ERB3-6','freq1200-2400_harm-sine15-30','freq1200-2400_harm-schr-1ERB15-30',...
            'freq2400-4800_harm-sine3-6','freq2400-4800_harm-schr-1ERB3-6','freq2400-4800_harm-sine15-30','freq2400-4800_harm-schr-1ERB15-30',...
            };
        
    case 'stax'
        plotorder = {'harm-f0125','noise-f0125','harm-f0600','noise-f0600'};
    case 'pitch_dp_v2'
        if any(us == [1 2 3 62 14 17])
            plotorder = {...
                'freq2000-4000_harm-sine-3-6','freq2000-4000_harm-schr-1ERB-3-6',...
                'freq2000-4000_harm-sine-15-30','freq2000-4000_harm-schr-1ERB-15-30',...
                'freq2000-4000_noise',...
                'freq2000-4000_harm-sine-masknoise-3-6',...
                'freq2000-4000_harm-sine-masknoise-15-30','freq2000-4000_harm-schr-1ERB-masknoise-15-30',...
                'freq2000-4000_noise-masknoise',...
                };
        else
            plotorder = {...
                'freq2000-4000_harm-sine-3-6',...
                'freq2000-4000_harm-sine-15-30',...
                'freq2000-4000_noise',...
                'freq2000-4000_harm-sine-masknoise-3-6',...
                'freq2000-4000_harm-sine-masknoise-15-30',...
                'freq2000-4000_noise-masknoise',...
                };
        end
    case 'amusia'
        plotorder = {'harm-f0167','noise-f0167','harm-f0667','noise-f0667'};
    case 'music_scram_familiar'
        plotorder = {...
            'familiar-music-intact-10sec', 'familiar-music-quilt-30ms-10sec',...
            'unfamiliar-music-intact-10sec','unfamiliar-music-quilt-30ms-10sec',...
            'bigband-intact-10sec','bigband-quilt-30ms-10sec',...
            'midi-intact-10sec','midi-pitch-scram-jitter-10sec','midi-onsetintervals-scram-jitter-10sec','midi-quilt-30ms-10sec',...
            'german-intact-10sec','german-quilt-30ms-10sec',...
            'env-pitch-10sec',...
            };
    case 'pitch_overlap'
        switch runtype
            case 'overlap_v2'
                plotorder{1} = 'freq1200-2400_harm-sine4-8';
                plotorder{2} = 'freq1200-2400_harm-sine15-30';
                plotorder{3} = 'freq1200-2400_noise';
                
                plotorder{4} = 'freq2400-4800_harm-sine4-8';
                plotorder{5} = 'freq2400-4800_harm-sine15-30';
                plotorder{6} = 'freq2400-4800_noise';
                
                plotorder{7} = '2note-lowhigh_harm4-8';
                plotorder{8} = '2note-highlow_harm4-8';
                plotorder{9} = '2note-lowhigh_harm15-30';
                plotorder{10} = '2note-highlow_harm15-30';
                plotorder{11} = '2note-lowhigh_noise';
                plotorder{12} = '2note-highlow_noise';
                
                plotorder{13} = '1note-lowfreq_harm4-8';
                plotorder{14} = '1note-lowfreq_harm15-30';
                plotorder{15} = '1note-lowfreq_noise';
                
                plotorder{16} = '1note-highfreq_harm4-8';
                plotorder{17} = '1note-highfreq_harm15-30';
                plotorder{18} = '1note-highfreq_noise';
                
                
            otherwise
                plotorder = read_conditions(exp, us, runtype, varargin{:});
        end
        
    case 'pitch_overlap_v3'
        
        switch runtype
            case 'overlap_v3'
                
                plotorder{1} = '1note-low_harm3-6';
                plotorder{2} = '1note-low_harm15-30';
                plotorder{3} = '1note-low_noise';
                
                plotorder{4} = '2note-octave-lowhigh_harm3-6';
                plotorder{5} = '2note-octave-lowhigh_harm15-30';
                plotorder{6} = '2note-octave-lowhigh_noise';
                
                plotorder{7} = '2note-overlap-low_harm3-6';
                plotorder{8} = '2note-overlap-low_harm15-30';
                plotorder{9} = '2note-overlap-low_noise';
                
                plotorder{10} = '1note-high_harm3-6';
                plotorder{11} = '1note-high_harm15-30';
                plotorder{12} = '1note-high_noise';
                
                plotorder{13} = '2note-octave-highlow_harm3-6';
                plotorder{14} = '2note-octave-highlow_harm15-30';
                plotorder{15} = '2note-octave-highlow_noise';
                
                plotorder{16} = '2note-overlap-high_harm3-6';
                plotorder{17} = '2note-overlap-high_harm15-30';
                plotorder{18} = '2note-overlap-high_noise';
                
                
            otherwise
                plotorder = read_conditions(exp, us, runtype, varargin{:});
        end
        
    case 'pitch_adapt_params'
        
        switch runtype
            case {'adapt_params_fixed_stim','adapt_params_variable_stim'}
                
                plotorder{1} = 'freq400_freqdiff0_rate4';
                plotorder{2} = 'freq400_freqdiff275_rate4';
                plotorder{3} = 'freq400_freqdiff550_rate4';
                plotorder{4} = 'freq400_freqdiff1100_rate4';
                
                plotorder{5} = 'freq400_freqdiff0_rate8';
                plotorder{6} = 'freq400_freqdiff275_rate8';
                plotorder{7} = 'freq400_freqdiff550_rate8';
                plotorder{8} = 'freq400_freqdiff1100_rate8';
                
                plotorder{9} = 'freq3200_freqdiff0_rate4';
                plotorder{10} = 'freq3200_freqdiff275_rate4';
                plotorder{11} = 'freq3200_freqdiff550_rate4';
                plotorder{12} = 'freq3200_freqdiff1100_rate4';
                
                plotorder{13} = 'freq3200_freqdiff0_rate8';
                plotorder{14} = 'freq3200_freqdiff275_rate8';
                plotorder{15} = 'freq3200_freqdiff550_rate8';
                plotorder{16} = 'freq3200_freqdiff1100_rate8';
                
                
                
            otherwise
                plotorder = read_conditions(exp, us, runtype, varargin{:});
        end
        
    case 'pitch_adapt_params_v2'
        
        switch runtype
            case {'adapt_params_variable_stim'}
                
                plotorder{1} = 'freq400_freqdiff0_sonelevel3';
                plotorder{2} = 'freq400_freqdiff275_sonelevel3';
                plotorder{3} = 'freq400_freqdiff550_sonelevel3';
                plotorder{4} = 'freq400_freqdiff1100_sonelevel3';
                
                plotorder{5} = 'freq400_freqdiff0_sonelevel8';
                plotorder{6} = 'freq400_freqdiff275_sonelevel8';
                plotorder{7} = 'freq400_freqdiff550_sonelevel8';
                plotorder{8} = 'freq400_freqdiff1100_sonelevel8';
                
                plotorder{9} = 'freq3200_freqdiff0_sonelevel3';
                plotorder{10} = 'freq3200_freqdiff275_sonelevel3';
                plotorder{11} = 'freq3200_freqdiff550_sonelevel3';
                plotorder{12} = 'freq3200_freqdiff1100_sonelevel3';
                
                plotorder{13} = 'freq3200_freqdiff0_sonelevel8';
                plotorder{14} = 'freq3200_freqdiff275_sonelevel8';
                plotorder{15} = 'freq3200_freqdiff550_sonelevel8';
                plotorder{16} = 'freq3200_freqdiff1100_sonelevel8';
                
                plotorder{17} = 'freq1131_freqdiff0_sonelevel3';
                plotorder{18} = 'freq1131_freqdiff3800_sonelevel3';
                plotorder{19} = 'freq1131_freqdiff0_sonelevel8';
                plotorder{20} = 'freq1131_freqdiff3800_sonelevel8';
                
                
                
            otherwise
                plotorder = read_conditions(exp, us, runtype, varargin{:});
        end
        
    otherwise
        plotorder = read_conditions(exp, us, runtype, varargin{:});
end
