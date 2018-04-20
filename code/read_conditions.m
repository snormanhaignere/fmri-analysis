function conds = read_conditions(exp,us,runtype,varargin)

% function conds = read_conditions(s,runtype)

switch exp
    
    case 'pitch_event'
        switch runtype
            case 'main'
                if us == 32
                    conds = {'harm-highfreq','harm-lowfreq','harm-unresolved','noise-highfreq','noise-lowfreq','targets'};
                else
                    conds = {'harm-highfreq','noise-highfreq','harm-unresolved','harm-spaced','harm-spaced-jittered','noise-spaced','targets'};
                end
                
            case 'localizer'
                
                if us == 32
                    conds = {'solo','nonwords'};
                else
                    conds = {'harm-highfreq','noise-highfreq','harm-lowfreq','harm-unresolved'};
                end
        end
        
    case 'pitch_resthr'
        
        conds = {'harm-1-2', 'harm-2-4', 'harm-3-6', 'harm-4-8', 'harm-5-10', 'harm-7-14', 'harm-9-18', 'harm-12-24', 'harm-15-30', 'harm-20-40', 'noise'};
        
    case 'pitch_f0'
        
        switch runtype
            case 'adapt6note'
                conds = {...
                    'f0same-freqsame-AA', 'f0same-freqsame-BB', 'f0same-freqsame-CC', 'f0same-freqsame-DD',...
                    'f0same-freqdiff-AB', 'f0same-freqdiff-BA', 'f0same-freqdiff-CD', 'f0same-freqdiff-DC',...
                    'f0diff-freqdiff-AC', 'f0diff-freqdiff-CA', 'f0diff-freqdiff-BD', 'f0diff-freqdiff-DB'};
            case 'localizer'
                conds = {'harm-f0125','noise-f0125','harm-f0275','noise-f0275','harm-f0600','noise-f0600'};
        end
        
    case 'pitch_f0inharm'
        
        conds = {...
            'harm-f0same-freqsame-AA', 'harm-f0same-freqsame-BB', 'harm-f0same-freqsame-CC', 'harm-f0same-freqsame-DD',...
            'harm-f0same-freqdiff-AB', 'harm-f0same-freqdiff-BA', 'harm-f0same-freqdiff-CD', 'harm-f0same-freqdiff-DC',...
            'harm-f0diff-freqdiff-AC', 'harm-f0diff-freqdiff-CA', 'harm-f0diff-freqdiff-BD', 'harm-f0diff-freqdiff-DB',...
            ...
            'inharm-f0same-freqsame-AA', 'inharm-f0same-freqsame-BB', 'inharm-f0same-freqsame-CC', 'inharm-f0same-freqsame-DD',...
            'inharm-f0same-freqdiff-AB', 'inharm-f0same-freqdiff-BA', 'inharm-f0same-freqdiff-CD', 'inharm-f0same-freqdiff-DC',...
            'inharm-f0diff-freqdiff-AC', 'inharm-f0diff-freqdiff-CA', 'inharm-f0diff-freqdiff-BD', 'inharm-f0diff-freqdiff-DB'...
            };
        
    case 'pitch_f0smallfreq'
        
        conds = {...
            'harm-low-f0same-freqsame','harm-low-f0same-freqdiff', 'harm-low-f0diff-freqdiff',...
            'inharm-low-f0same-freqsame','inharm-low-f0same-freqdiff', 'inharm-low-f0diff-freqdiff',...
            'harm-high-f0same-freqsame','harm-high-f0same-freqdiff', 'harm-high-f0diff-freqdiff',...
            'inharm-high-f0same-freqsame','inharm-high-f0same-freqdiff', 'inharm-high-f0diff-freqdiff',...
            };
        
    case 'pitch_dp'
        
        conds = {'harm-lowfreq','noise-lowfreq','harm-highfreq','noise-highfreq','harm-unres','harm-unres-schroeder-neg',...
            'harm-lowfreq-lowmask','noise-lowfreq-lowmask','harm-highfreq-lowmask','noise-highfreq-lowmask','harm-unres-lowmask','harm-unres-schroeder-neg-lowmask'...
            'harm-highfreq-schroeder-negN-1'};
        
    case 'pitch_natural'
        
        conds = {'harm-lowfreq','noise-lowfreq','harm-highfreq','noise-highfreq','harm-unres','animals-pitch','env-pitch','env-nopitch','solo','drums','songs','speech-normal','speech-whisper'};
        if optInputs(varargin,'naturalconds')
            conds = conds(6:13);
        end
        
    case 'pitch_contour'
        
        switch runtype
            case 'main'
                conds = {'harm-freqsame-melsame', 'harm-freqdiff-melsame', 'harm-freqsame-meldiff', 'harm-freqdiff-meldiff'};
            case 'localizer'
                conds = {'harm-f0125','noise-f0125','harm-f0275','noise-f0275','harm-f0600','noise-f0600'};
        end
        
    case 'pitch_f02by2'
        
        conds = {'harm-f0ref-f0same-freqsame','harm-f0ref-f0same-freqdiff','harm-f0ref-f0diff-freqsame','harm-f0ref-f0diff-freqdiff'};
        
    case 'mspec'
        
        if any(us == [33 34])
            conds = {'orchestra','solo','sentences','nonwords','animals','events','textures','songs','lyrics','midi-intact','midi-scram','drums-intact','drums-scram'};
        else
            conds = {'orchestra','solo','sentences','nonwords','animals-pitch','env-pitch','env-nopitch','songs','lyrics','drums-intact','drums-scram','mel-intact','mel-scram','mel-interp','mel-noise'};
        end
        
    case 'tonotopy_rate'
        
        conds = {...
            'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16',...
            'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16',...
            'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16',...
            'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16',...
            };
        
        
    case 'sepi'
        
        conds = {'pitch','noise','highfreq','lowfreq','sentences','nonwords','solo'};
        
    case 'pitch_resthr_v2'
        
        switch runtype
            
            case 'main'
                
                conds = {...
                    'freq1200-2400_harm-sine3-6','freq1200-2400_harm-sine4-8',  'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12', 'freq1200-2400_harm-sine8-16',...
                    'freq1200-2400_harm-sine10-20','freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30','freq1200-2400_noise',...
                    'freq2400-4800_harm-sine3-6','freq2400-4800_harm-sine4-8',  'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12', 'freq2400-4800_harm-sine8-16',...
                    'freq2400-4800_harm-sine10-20','freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30','freq2400-4800_noise',...
                    };
                
            case 'localizer'
                
                conds = {'freq200','freq400','freq800','freq1600','freq3200','freq6400','harm-f0253','noise-f0253'};
                
        end
        
    case 'amusia'
        
        if any(us == [43 44])
            conds = {'harm','highfreq'};
        else
            conds = {'harm-f0167','harm-f0667','noise-f0167','noise-f0667'};
        end
        
    case 'pitch_resthr_v3'
        
        switch runtype
            
            case 'main'
                
                conds = {...
                    'freq1200-2400_harm-sine3-6',     'freq1200-2400_harm-sine4-8',       'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12', ...
                    'freq1200-2400_harm-sine8-16',    'freq1200-2400_harm-sine10-20',     'freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30',...
                    'freq1200-2400_harm-schr-1ERB3-6','freq1200-2400_harm-schr-1ERB15-30','freq1200-2400_noise',...
                    'freq2400-4800_harm-sine3-6',     'freq2400-4800_harm-sine4-8',       'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12',...
                    'freq2400-4800_harm-sine8-16',    'freq2400-4800_harm-sine10-20',     'freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30',...
                    'freq2400-4800_harm-schr-1ERB3-6','freq2400-4800_harm-schr-1ERB15-30','freq2400-4800_noise',...
                    };
                
        end
        
    case 'pitch_resthr_v4'
        
        switch runtype
            
            case 'main'
                
                conds = {...
                    'freq1200-2400_harm-sine3-6',     'freq1200-2400_harm-sine4-8',       'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12', ...
                    'freq1200-2400_harm-sine8-16',    'freq1200-2400_harm-sine10-20',     'freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30',...
                    'freq1200-2400_noise',...
                    'freq2400-4800_harm-sine3-6',     'freq2400-4800_harm-sine4-8',       'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12',...
                    'freq2400-4800_harm-sine8-16',    'freq2400-4800_harm-sine10-20',     'freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30',...
                    'freq2400-4800_noise',...
                    };
                
            case 'localizer'
                
                conds = {'freq200','freq400','freq800','freq1600','freq3200','freq6400','harm-f0333','noise-f0333'};
                
        end
        
    case 'pitch_misc'
        
        switch runtype
            
            case 'schroeder'
                
                conds = {...
                    'freq1200-2400_harm-schr-1ERB3-6','freq1200-2400_harm-sine15-30','freq1200-2400_harm-schr-1ERB15-30','freq1200-2400_noise',...
                    'freq2400-4800_harm-schr-1ERB3-6','freq2400-4800_harm-sine15-30','freq2400-4800_harm-schr-1ERB15-30','freq2400-4800_noise',...
                    };
                
        end
        
    case 'music_12channel'
        conds = {'solo-long','solo-short','solo-long-highfreq','env-pitch','env-nopitch','midi-intact','midi-pitch-scram-jitter','midi-onsetintervals-scram-jitter'};
    case 'stax'
        conds = {'harm-f0125','harm-f0600','noise-f0125','noise-f0600'};
        
    case 'multiband'
        conds = {'harm-f0125','harm-f0600','noise-f0125','noise-f0600'};
        
    case 'pitch_dp_v2'
        if any(us == [61 106])
            conds = {...
                'freq2000-4000_harm-sine-3-6','freq2000-4000_harm-sine-15-30','freq2000-4000_noise',...
                'freq2000-4000_harm-sine-masknoise-3-6','freq2000-4000_harm-sine-masknoise-15-30','freq2000-4000_noise-masknoise',...
                };
        else
            conds = {...
                'freq2000-4000_harm-sine-3-6','freq2000-4000_harm-sine-15-30','freq2000-4000_noise',...
                'freq2000-4000_harm-sine-masknoise-3-6','freq2000-4000_harm-sine-masknoise-15-30','freq2000-4000_noise-masknoise',...
                'freq2000-4000_harm-schr-1ERB-3-6','freq2000-4000_harm-schr-1ERB-15-30','freq2000-4000_harm-schr-1ERB-masknoise-15-30',...
                };
        end
        
    case 'music_scram'
        
        conds = {...
            'midi-intact-10sec','midi-pitch-scram-jitter-10sec','midi-onsetintervals-scram-jitter-10sec','midi-quilt-30ms-10sec',...
            'german-intact-10sec','german-quilt-30ms-10sec',...
            'bigband-intact-10sec','bigband-quilt-30ms-10sec',...
            'env-pitch-10sec','env-nopitch-10sec',...
            };
        
    case 'pitch_f0adapt'
        
        switch runtype
            
            case 'f0adapt'
                
                conds = {...
                    'harm-0-f0same-freqdiff','harm-0-f0diff-freqdiff',...
                    'harm-10-f0same-freqdiff','harm-10-f0diff-freqdiff',...
                    'harm-inf-f0same-freqdiff','harm-inf-f0diff-freqdiff',...
                    'inharm-0-f0same-freqdiff','inharm-0-f0diff-freqdiff',...
                    'inharm-10-f0same-freqdiff','inharm-10-f0diff-freqdiff',...
                    'inharm-inf-f0same-freqdiff','inharm-inf-f0diff-freqdiff',...
                    };
                
            case 'localizer'
                
                conds = {'freq200', 'freq400', 'freq800', 'freq1600', 'freq3200', 'freq6400', 'harm-f0333', 'noise-f0333'};
                
        end
        
    case 'pitch_f0adapt_v3'
        
        switch runtype
            
            case 'f0adapt'
                
                clear conds
                conds{1} = 'freq200_test_harm_f0same-freqsame';
                conds{2} = 'freq200_test_harm_f0same-freqdiff';
                conds{3} = 'freq200_test_harm_f0diff-freqdiff';
                
                conds{4} = 'freq200_test_inharm_f0same-freqsame';
                conds{5} = 'freq200_test_inharm_f0same-freqdiff';
                conds{6} = 'freq200_test_inharm_f0diff-freqdiff';
                
                conds{7} = 'freq800_test_harm_f0same-freqsame';
                conds{8} = 'freq800_test_harm_f0same-freqdiff';
                conds{9} = 'freq800_test_harm_f0diff-freqdiff';
                
                conds{10} = 'freq800_test_inharm_f0same-freqsame';
                conds{11} = 'freq800_test_inharm_f0same-freqdiff';
                conds{12} = 'freq800_test_inharm_f0diff-freqdiff';
                
            case 'localizer'
                
                conds = {'freq200', 'freq400', 'freq800', 'freq1600', 'freq3200', 'freq6400', 'harm-f0333', 'noise-f0333'};
                
        end
        
    case 'music_scram_familiar'
        
        switch runtype
            case 'main'
                clear conds;
                conds{1} = 'midi-intact-10sec';
                conds{2} = 'midi-pitch-scram-jitter-10sec';
                conds{3} = 'midi-onsetintervals-scram-jitter-10sec';
                conds{4} = 'midi-quilt-30ms-10sec';
                conds{5} = 'german-intact-10sec';
                conds{6} = 'german-quilt-30ms-10sec';
                conds{7} = 'familiar-music-intact-10sec';
                conds{8} = 'familiar-music-quilt-30ms-10sec';
                conds{9} = 'unfamiliar-music-intact-10sec';
                conds{10} = 'unfamiliar-music-quilt-30ms-10sec';
                conds{11} = 'env-pitch-10sec';
                conds{12} = 'bigband-intact-10sec';
                conds{13} = 'bigband-quilt-30ms-10sec';
            otherwise
                fprintf('Error in read_conditions.m: No valid runtype\n');
                drawnow; keyboard;
        end
        
    case 'pitch_overlap'
        
        switch runtype
            
            case 'overlap'
                conds{1} = 'freq1200-2400_harm-sine3-6';
                conds{2} = 'freq2400-4800_harm-sine3-6';
                conds{3} = 'freq1200-2400_harm-sine15-30';
                conds{4} = 'freq2400-4800_harm-sine15-30';
                conds{5} = 'freq1200-2400_noise';
                conds{6} = 'freq2400-4800_noise';
                
                conds{7} = '2note-lowhigh_harm3-6';
                conds{8} = '2note-highlow_harm3-6';
                conds{9} = '2note-lowhigh_harm15-30';
                conds{10} = '2note-highlow_harm15-30';
                conds{11} = '2note-lowhigh_noise';
                conds{12} = '2note-highlow_noise';
                
                conds{13} = '1note-lowfreq_harm3-6';
                conds{14} = '1note-highfreq_harm3-6';
                conds{15} = '1note-lowfreq_harm15-30';
                conds{16} = '1note-highfreq_harm15-30';
                conds{17} = '1note-lowfreq_noise';
                conds{18} = '1note-highfreq_noise';
                
            case 'overlap_v2'
                conds{1} = 'freq1200-2400_harm-sine4-8';
                conds{2} = 'freq2400-4800_harm-sine4-8';
                conds{3} = 'freq1200-2400_harm-sine15-30';
                conds{4} = 'freq2400-4800_harm-sine15-30';
                conds{5} = 'freq1200-2400_noise';
                conds{6} = 'freq2400-4800_noise';
                
                conds{7} = '2note-lowhigh_harm4-8';
                conds{8} = '2note-highlow_harm4-8';
                conds{9} = '2note-lowhigh_harm15-30';
                conds{10} = '2note-highlow_harm15-30';
                conds{11} = '2note-lowhigh_noise';
                conds{12} = '2note-highlow_noise';
                
                conds{13} = '1note-lowfreq_harm4-8';
                conds{14} = '1note-highfreq_harm4-8';
                conds{15} = '1note-lowfreq_harm15-30';
                conds{16} = '1note-highfreq_harm15-30';
                conds{17} = '1note-lowfreq_noise';
                conds{18} = '1note-highfreq_noise';
                
            case 'tonotopy'
                conds = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
        end
        
    case 'pitch_overlap_v2'
        
        switch runtype
            
            case 'overlap_v2'
                conds{1} = 'freq1200-2400_harm-sine4-8';
                conds{2} = 'freq1200-2400_harm-sine15-30';
                conds{3} = 'freq1200-2400_noise';
                
                conds{4} = '2note-octave-lowhigh_harm4-8';
                conds{5} = '2note-octave-lowhigh_harm15-30';
                conds{6} = '2note-octave-lowhigh_noise';
                
                conds{7} = '2note-overlap-low_harm4-8';
                conds{8} = '2note-overlap-low_harm15-30';
                conds{9} = '2note-overlap-low_noise';
                
                conds{10} = '1note-low_harm4-8';
                conds{11} = '1note-low_harm15-30';
                conds{12} = '1note-low_noise';
                
                conds{13} = 'freq2400-4800_harm-sine4-8';
                conds{14} = 'freq2400-4800_harm-sine15-30';
                conds{15} = 'freq2400-4800_noise';
                
                conds{16} = '2note-octave-highlow_harm4-8';
                conds{17} = '2note-octave-highlow_harm15-30';
                conds{18} = '2note-octave-highlow_noise';
                
                conds{19} = '2note-overlap-high_harm4-8';
                conds{20} = '2note-overlap-high_harm15-30';
                conds{21} = '2note-overlap-high_noise';
                
                conds{22} = '1note-high_harm4-8';
                conds{23} = '1note-high_harm15-30';
                conds{24} = '1note-high_noise';
        end
        
    case 'pitch_overlap_v3'
        
        switch runtype
            
            case 'localizer'
                
                conds{1} = 'freq200';
                conds{2} = 'freq400';
                conds{3} = 'freq800';
                conds{4} = 'freq1600';
                conds{5} = 'freq3200';
                conds{6} = 'freq6400';
                
                conds{7} = 'freq1200-2400_harm-sine3-6';
                conds{8} = 'freq1200-2400_harm-sine15-30';
                conds{9} = 'freq1200-2400_noise';
                conds{10} = 'freq2400-4800_harm-sine3-6';
                conds{11} = 'freq2400-4800_harm-sine15-30';
                conds{12} = 'freq2400-4800_noise';
                
            case 'overlap_v3'
                
                conds{1} = '2note-octave-lowhigh_harm3-6';
                conds{2} = '2note-octave-lowhigh_harm15-30';
                conds{3} = '2note-octave-lowhigh_noise';
                
                conds{4} = '2note-overlap-low_harm3-6';
                conds{5} = '2note-overlap-low_harm15-30';
                conds{6} = '2note-overlap-low_noise';
                
                conds{7} = '1note-low_harm3-6';
                conds{8} = '1note-low_harm15-30';
                conds{9} = '1note-low_noise';
                
                conds{10} = '2note-octave-highlow_harm3-6';
                conds{11} = '2note-octave-highlow_harm15-30';
                conds{12} = '2note-octave-highlow_noise';
                
                conds{13} = '2note-overlap-high_harm3-6';
                conds{14} = '2note-overlap-high_harm15-30';
                conds{15} = '2note-overlap-high_noise';
                
                conds{16} = '1note-high_harm3-6';
                conds{17} = '1note-high_harm15-30';
                conds{18} = '1note-high_noise';
                
                
            case 'overlap_v3_lowconds'
                
                conds{1} = '2note-octave-lowhigh_harm3-6';
                conds{2} = '2note-octave-lowhigh_harm15-30';
                conds{3} = '2note-octave-lowhigh_noise';
                
                conds{4} = '2note-overlap-low_harm3-6';
                conds{5} = '2note-overlap-low_harm15-30';
                conds{6} = '2note-overlap-low_noise';
                
                conds{7} = '1note-low_harm3-6';
                conds{8} = '1note-low_harm15-30';
                conds{9} = '1note-low_noise';
                
            case 'overlap_v3_highconds'
                
                conds{1} = '2note-octave-highlow_harm3-6';
                conds{2} = '2note-octave-highlow_harm15-30';
                conds{3} = '2note-octave-highlow_noise';
                
                conds{4} = '2note-overlap-high_harm3-6';
                conds{5} = '2note-overlap-high_harm15-30';
                conds{6} = '2note-overlap-high_noise';
                
                conds{7} = '1note-high_harm3-6';
                conds{8} = '1note-high_harm15-30';
                conds{9} = '1note-high_noise';
                
        end
        
    case 'pitch_f0adapt_v2'
        
        conds = {...
            'f0200-harm-inf-f0same-freqsame','f0200-harm-inf-f0same-freqdiff','f0200-harm-inf-f0diff-freqdiff',...
            'f0200-inharm-inf-f0same-freqsame','f0200-inharm-inf-f0same-freqdiff','f0200-inharm-inf-f0diff-freqdiff',...
            ...
            'f0400-harm-inf-f0same-freqsame','f0400-harm-inf-f0same-freqdiff','f0400-harm-inf-f0diff-freqdiff',...
            'f0400-inharm-inf-f0same-freqsame','f0400-inharm-inf-f0same-freqdiff','f0400-inharm-inf-f0diff-freqdiff',...
            ...
            'f0800-harm-inf-f0same-freqsame','f0800-harm-inf-f0same-freqdiff','f0800-harm-inf-f0diff-freqdiff',...
            'f0800-inharm-inf-f0same-freqsame','f0800-inharm-inf-f0same-freqdiff','f0800-inharm-inf-f0diff-freqdiff',...
            };
        
    case 'pitch_adapt_params'
        % conditions
        conds{1} = 'freq400_freqdiff0_rate4';
        conds{2} = 'freq400_freqdiff0_rate8';
        
        conds{3} = 'freq3200_freqdiff0_rate4';
        conds{4} = 'freq3200_freqdiff0_rate8';
        
        conds{5} = 'freq400_freqdiff275_rate4';
        conds{6} = 'freq400_freqdiff275_rate8';
        
        conds{7} = 'freq3200_freqdiff275_rate4';
        conds{8} = 'freq3200_freqdiff275_rate8';
        
        conds{9} = 'freq400_freqdiff550_rate4';
        conds{10} = 'freq400_freqdiff550_rate8';
        
        conds{11} = 'freq3200_freqdiff550_rate4';
        conds{12} = 'freq3200_freqdiff550_rate8';
        
        conds{13} = 'freq400_freqdiff1100_rate4';
        conds{14} = 'freq400_freqdiff1100_rate8';
        
        conds{15} = 'freq3200_freqdiff1100_rate4';
        conds{16} = 'freq3200_freqdiff1100_rate8';
        
    case 'naturalsound'
        
        switch runtype
            case {'quilting_v1_continuous'}
                clear conds;
                conds{1} = 'music-intact';
                conds{2} = 'music-quilt-31ms';
                conds{3} = 'music-quilt-63ms';
                conds{4} = 'music-quilt-125ms';
                conds{5} = 'music-quilt-250ms';
                conds{6} = 'music-quilt-500ms';
                conds{7} = 'music-quilt-1000ms';
                conds{8} = 'music-quilt-2000ms';
                conds{9} = 'german-intact';
                conds{10} = 'german-quilt-31ms';
                conds{11} = 'german-quilt-63ms';
                conds{12} = 'german-quilt-125ms';
                conds{13} = 'german-quilt-250ms';
                conds{14} = 'german-quilt-500ms';
                conds{15} = 'german-quilt-1000ms';
                conds{16} = 'german-quilt-2000ms';
            case {'main_v2_combined','main_combined'}
                conds = strrep(mydir([params('rootdir') exp '/stims/final_stims60/']),'.wav','')';
            case {'main_v3_combined','main_v3_combined_raw'}
                conds = strrep(mydir([params('rootdir') exp '/stims/final_stims165/']),'.wav','')';
                x = regexp(conds,'(\d)*','match');
                [~, xi] = sort(str2double(cat(1,x{:})),'ascend');
                conds = conds(xi);
            case {'main_v3'}
                r = varargin{optInputs(varargin, 'run')+1};
                b = read_timings(exp,us,runtype,r,'nofix');
                conds = b.conds';
            case {'texture_combined'}
                x = load([params('rootdir') exp '/scripts_exp/selected_stimuli_36stims.mat']);
                stimulus_types = {'natural-1','natural-2','marginals','modpower','fullmodel'};
                conds = cell(length(x.stim_names)*length(stimulus_types),1);
                for i = 1:length(stimulus_types)
                    for j = 1:length(x.stim_names)
                        conds{j + (i-1)*36} = [stimulus_types{i} '-' x.stim_names{j}];
                    end
                end
            case {'natsoundloc'}
                conds = {'english','german-intact','german-quilt','bigband-intact','bigband-quilt','env-pitch'};
            case 'localizer'
                conds = {'freq200','freq400','freq800','freq1600','freq3200','freq6400','harm-f0333','noise-f0333'};
            case {'pitch_localizer'}
                conds = {'harm-f0333','noise-f0333'};
            case {'speech_localizer'}
                conds = {'german-intact','german-quilt'};
            case {'scrambling'}
                conds = {'german-v2-intact','german-v2-quilt','bigband-v2-intact','bigband-v2-quilt','midi-intact','midi-pitch-scram','midi-pitch-rhythm-scram','midi-quilt'};
            case {'scrambling_russian'}
                conds = {'russian-intact','russian-quilt','bigband-v2-intact','bigband-v2-quilt','midi-intact','midi-pitch-scram','midi-pitch-rhythm-scram','midi-quilt'};
            case {'texture'}
                conds = {'natural-1','natural-2','marginals','modpower','fullmodel'};
            case {'spectrotemporal'}
                conds = {'originals-1','originals-2','cochlear-marginals','tempmod-marginals','specmod-marginals','spectempmod'};
            case {'spectrotemporal_v2'}
                conds = {'originals-1','originals-2','cochlear-marginals','tempmod','spectempmod'};
            case {'spectrotemporal_v2_combined','spectrotemporal_v2_combined_raw'}
                x = load([params('rootdir') 'naturalsound/scripts_exp/spectrotemporal_stims_v4.mat']);
                stimulus_types = {'originals-1','originals-2','cochlear-marginals','tempmod','spectempmod'};
                conds = cell(length(x.stims)*length(stimulus_types),1);
                for i = 1:length(stimulus_types)
                    for j = 1:length(x.stims)
                        conds{j + (i-1)*length(x.stims)} = [stimulus_types{i} '-si' num2str(j) '-' x.stims{j}];
                    end
                end
            case {'spectrotemporal_combined','spectrotemporal_combined_raw'}
                x = load([params('rootdir') 'naturalsound/scripts_exp/spectrotemporal_stims.mat']);
                stimulus_types = {'originals-1','originals-2','cochlear-marginals','tempmod-marginals','specmod-marginals','spectempmod'};
                conds = cell(length(x.stims)*length(stimulus_types),1);
                for i = 1:length(stimulus_types)
                    for j = 1:length(x.stims)
                        conds{j + (i-1)*length(x.stims)} = [stimulus_types{i} '-si' num2str(j) '-' x.stims{j}];
                    end
                end
        end
        
    case 'pitch_adapt_params_v2'
        
        conds{1} = 'freq400_freqdiff0_sonelevel3';
        conds{2} = 'freq400_freqdiff0_sonelevel8';
        
        conds{3} = 'freq3200_freqdiff0_sonelevel3';
        conds{4} = 'freq3200_freqdiff0_sonelevel8';
        
        conds{5} = 'freq1131_freqdiff0_sonelevel3';
        conds{6} = 'freq1131_freqdiff0_sonelevel8';
        
        conds{7} = 'freq400_freqdiff275_sonelevel3';
        conds{8} = 'freq400_freqdiff275_sonelevel8';
        
        conds{9} = 'freq3200_freqdiff275_sonelevel3';
        conds{10} = 'freq3200_freqdiff275_sonelevel8';
        
        conds{11} = 'freq400_freqdiff550_sonelevel3';
        conds{12} = 'freq400_freqdiff550_sonelevel8';
        
        conds{13} = 'freq3200_freqdiff550_sonelevel3';
        conds{14} = 'freq3200_freqdiff550_sonelevel8';
        
        conds{15} = 'freq400_freqdiff1100_sonelevel3';
        conds{16} = 'freq400_freqdiff1100_sonelevel8';
        
        conds{17} = 'freq3200_freqdiff1100_sonelevel3';
        conds{18} = 'freq3200_freqdiff1100_sonelevel8';
        
        conds{19} = 'freq1131_freqdiff3800_sonelevel3';
        conds{20} = 'freq1131_freqdiff3800_sonelevel8';
        
    case 'pitch_adapt_params_v3'
        
        % conditions
        conds{1} = 'freq200_test_freqdiff0';
        conds{2} = 'freq200_test_freqdiff275';
        conds{3} = 'freq200_test_freqdiff550';
        conds{4} = 'freq200_test_freqdiff1100';
        conds{5} = 'freq3200_test_freqdiff0';
        conds{6} = 'freq3200_test_freqdiff275';
        conds{7} = 'freq3200_test_freqdiff550';
        conds{8} = 'freq3200_test_freqdiff1100';
        conds{9}  = 'freq200_adapt';
        conds{10} = 'freq3200_adapt';
        conds{11} = 'freq200_topup';
        conds{12} = 'freq3200_topup';
        
    case 'pitch_f0adapt_v4'
        
        switch runtype
            
            case 'f0adapt'
                
                clear conds
                conds{1} = 'freq297_test_harm_f0same-freqsame';
                conds{2} = 'freq297_test_harm_f0same-freqdiff';
                conds{3} = 'freq297_test_harm_f0diff-freqdiff';
                
                conds{4} = 'freq297_test_inharm_f0same-freqsame';
                conds{5} = 'freq297_test_inharm_f0same-freqdiff';
                conds{6} = 'freq297_test_inharm_f0diff-freqdiff';
                
                conds{7} = 'freq440_test_harm_f0same-freqsame';
                conds{8} = 'freq440_test_harm_f0same-freqdiff';
                conds{9} = 'freq440_test_harm_f0diff-freqdiff';
                
                conds{10} = 'freq440_test_inharm_f0same-freqsame';
                conds{11} = 'freq440_test_inharm_f0same-freqdiff';
                conds{12} = 'freq440_test_inharm_f0diff-freqdiff';
                
            case 'localizer'
                
                conds = {'freq200', 'freq400', 'freq800', 'freq1600', 'freq3200', 'freq6400', 'harm-f0333', 'noise-f0333'};
                
        end
        
    case 'pitch_f0adapt_v5'
        
        switch runtype
            
            case 'f0adapt'
                
                clear conds
                conds{1} = 'freq297_test_harm_f0same-freqsame';
                conds{2} = 'freq297_test_harm_f0same-freqdiff';
                conds{3} = 'freq297_test_harm_f0diff-freqdiff';
                
                conds{4} = 'freq297_test_inharm_f0same-freqsame';
                conds{5} = 'freq297_test_inharm_f0same-freqdiff';
                conds{6} = 'freq297_test_inharm_f0diff-freqdiff';
                
                conds{7} = 'freq200_test_harm_f0same-freqsame';
                conds{8} = 'freq200_test_harm_f0same-freqdiff';
                conds{9} = 'freq200_test_harm_f0diff-freqdiff';
                
                conds{10} = 'freq200_test_inharm_f0same-freqsame';
                conds{11} = 'freq200_test_inharm_f0same-freqdiff';
                conds{12} = 'freq200_test_inharm_f0diff-freqdiff';
                
            case 'localizer'
                
                conds = {'freq200', 'freq400', 'freq800', 'freq1600', 'freq3200', 'freq6400', 'harm-f0333', 'noise-f0333'};
                
        end
        
    case 'pitch_f0adapt_v6'
        
        switch runtype
            
            case 'freqf0adapt2'
                
                clear conds
                conds{1} = 'freq300_puretone_freqdiff-0';
                conds{2} = 'freq300_puretone_freqdiff-2p75';
                conds{3} = 'freq300_puretone_freqdiff-5p5';
                conds{4} = 'freq300_puretone_freqdiff-11';
                
                conds{5} = 'freq3200_puretone_freqdiff-0';
                conds{6} = 'freq3200_puretone_freqdiff-2p75';
                conds{7} = 'freq3200_puretone_freqdiff-5p5';
                conds{8} = 'freq3200_puretone_freqdiff-11';
                
                
            case 'localizer'
                
                conds = {'freq200', 'freq400', 'freq800', 'freq1600', 'freq3200', 'freq6400', 'harm-f0333', 'noise-f0333'};
                
        end
        
    case 'pitch_f0adapt_v7'
        
        switch runtype
            
            case 'freqf0adapt'
                
                clear conds
                conds{1} = 'freq300_puretone_freqdiff-0';
                conds{2} = 'freq300_puretone_freqdiff-2p75';
                conds{3} = 'freq300_puretone_freqdiff-5p5';
                conds{4} = 'freq300_puretone_freqdiff-11';
                
                conds{5} = 'freq300_harm_freqdiff-0';
                conds{6} = 'freq300_harm_freqdiff-5p5';
                
                conds{7} = 'freq300_inharm_freqdiff-0';
                conds{8} = 'freq300_inharm_freqdiff-5p5';
                
            case 'localizer'
                
                conds = {'freq200', 'freq400', 'freq800', 'freq1600', 'freq3200', 'freq6400', 'harm-f0333', 'noise-f0333'};
                
        end
        
    case 'tonotopy_monkey'
        switch runtype
            case 'localizer_1000msTA'
                conds = {'freq200', 'freq400', 'freq800', 'freq1600', 'freq3200', 'freq6400', 'freq12800'};
        end
        
    case 'pitch_localizer_monkey'
        switch runtype
            case {'pitchloc2','pitchloc2_combined_raw','pitchloc2_combined_split'}
                conds = {'harm-f0100', 'noise-f0100','harm-f0200', 'noise-f0200','harm-f0400', 'noise-f0400','harm-f0800', 'noise-f0800','harm-f01600', 'noise-f01600'};
            case 'pitchloc3'
                conds{1} = 'harm-res-freq1200-2400';
                conds{2} = 'harm-unres-freq1200-2400';
                conds{3} = 'noise-freq1200-2400';
                conds{4} = 'harm-res-freq2400-4800';
                conds{5} = 'harm-unres-freq2400-4800';
                conds{6} = 'noise-freq2400-4800';
                conds{7} = 'harm-res-freq4800-9600';
                conds{8} = 'harm-unres-freq4800-9600';
                conds{9} = 'noise-freq4800-9600';
            case 'pitchloc4'
                conds{1} = 'harm-res-freq1200-2400';
                conds{2} = 'harm-unres-freq1200-2400';
                conds{3} = 'noise-freq1200-2400';
                conds{4} = 'noise-mod-freq1200-2400';
                conds{5} = 'harm-res-freq2400-4800';
                conds{6} = 'harm-unres-freq2400-4800';
                conds{7} = 'noise-freq2400-4800';
                conds{8} = 'noise-mod-freq2400-4800';
            case {'pitchloc6','pitchloc6_combined_raw','pitchloc6_combined_split'}
                conds = {'harm-f0200-80dB', 'noise-f0200-75dB', 'harm-f0200-75dB', 'noise-f0200-70dB', 'harm-f0200-70dB', 'noise-f0200-65dB'};
                
            case {'pitchloc7_70dB', 'pitchloc7_75dB', 'pitchloc7_80dB', ...
                    'pitchloc7_70-75-80dB', 'pitchloc7_70-75-80dB_combined_raw',...
                    'pitchloc7_70-75-80dB_combined_split'}
                
                switch runtype
                    case 'pitchloc7_70dB'
                        levels = 70;
                    case 'pitchloc7_75dB'
                        levels = 75;
                    case 'pitchloc7_80dB'
                        levels = 80;
                    case {'pitchloc7_70-75-80dB',...
                            'pitchloc7_70-75-80dB_combined_raw',...
                            'pitchloc7_70-75-80dB_combined_split'}
                        levels = [70 75 80];
                    otherwise
                        error('No matching runtype for %s\n', runtype);
                end
                
                clear conds;
                count = 0;
                stimtype = {'harm','noise'};
                for k = 1:length(levels)
                    for i = 1:length(stimtype)
                        
                        f0s = {'100','200','400','800','1600'};
                        for j = 1:length(f0s)
                            
                            count = count+1;
                            conds{count} = [stimtype{i} '-f0' f0s{j} '-' num2str(levels(k)) 'dB']; %#ok<AGROW>
                        end
                    end
                end
                
            otherwise
                error('No matching runtype for %s\n', runtype);
        end
        
    case 'voice_localizer_monkey'
        switch runtype
            case {'mvocs_pitch','mvocs_pitch_combined_raw','mvocs_pitch_combined_split'}
                conds = {'pitched', 'unpitched-matched','pitched-noisevoc', 'unpitched-matched-noisevoc'};
            case {'mvocs_pitch_v2','mvocs_pitch_v2_combined_raw','mvocs_pitch_v2_combined_split'}
                conds = {'pitched', 'unpitched-matched','pitched-straight-harmvoc', 'pitched-straight-noisevoc','pitched-straight-noisevoc-5dB'};
            case {'mvocs_pitch_v3','mvocs_pitch_v3_combined_raw','mvocs_pitch_v3_combined_split'}
                conds = {'pitched-straight-harmvoc', 'pitched-straight-noisevoc-5dB'};
            case {'mvocs_pitch_v4','mvocs_pitch_v4_combined_raw','mvocs_pitch_v4_combined_split'}
                conds = {'pitched-straight-harmvoc-80dB', 'pitched-straight-noisevoc-75dB',...
                    'pitched-straight-harmvoc-75dB', 'pitched-straight-noisevoc-70dB','pitched-straight-harmvoc-70dB', 'pitched-straight-noisevoc-65dB'};
            case {'mvocs_pitch_v5', 'mvocs_pitch_v5_combined_raw', 'mvocs_pitch_v5_combined_split'}
                conds = {...
                    'pitched-straight-harmvoc-80dB', 'pitched-straight-noisevoc-80dB', ...
                    'pitched-straight-harmvoc-75dB', 'pitched-straight-noisevoc-75dB', ...
                    'pitched-straight-harmvoc-70dB', 'pitched-straight-noisevoc-70dB', ...
                    'pitched-straight-harmvoc-65dB', 'pitched-straight-noisevoc-65dB'};
        end
        
    case 'voice_localizer_human'
        switch runtype
            case {'mvocs_pitch_v2','mvocs_pitch_v2_combined_raw','mvocs_pitch_v2_combined_split'}
                conds = {'pitched', 'unpitched-matched','pitched-straight-harmvoc', 'pitched-straight-noisevoc','pitched-straight-noisevoc-5dB'};
            case {'mvocs_pitch_v4','mvocs_pitch_v4_combined_raw','mvocs_pitch_v4_combined_split'}
                conds = {'pitched-straight-harmvoc-80dB', 'pitched-straight-noisevoc-80dB','pitched-straight-harmvoc-75dB', 'pitched-straight-noisevoc-75dB','pitched-straight-harmvoc-70dB', 'pitched-straight-noisevoc-70dB','pitched-straight-harmvoc-65dB', 'pitched-straight-noisevoc-65dB'};
        end
        
    case 'resting_monkey'
        switch runtype
            case {'rest','rest_combined_raw','rest_combined_split'}
                conds = {'harm-f0100', 'noise-f0100','harm-f0200', 'noise-f0200','harm-f0400', 'noise-f0400','harm-f0800', 'noise-f0800','harm-f01600', 'noise-f01600'};
        end
        
    case 'pitch_localizer_human'
        switch runtype
            case {'pitchloc2','pitchloc2_combined_raw','pitchloc2_combined_split'}
                conds = {...
                    'harm-f0100',  'noise-dichotic-f0100', 'noise-diotic-f0100', ...
                    'harm-f0200',  'noise-dichotic-f0200', 'noise-diotic-f0200', ...
                    'harm-f0400',  'noise-dichotic-f0400', 'noise-diotic-f0400', ...
                    'harm-f0800',  'noise-dichotic-f0800', 'noise-diotic-f0800', ...
                    'harm-f01600', 'noise-dichotic-f01600', 'noise-diotic-f01600', ...
                    };
            case {'pitchloc6','pitchloc6_combined_raw','pitchloc6_combined_split'}
                conds = {'harm-f0200-80dB', 'noise-f0200-80dB', 'harm-f0200-75dB', 'noise-f0200-75dB', 'harm-f0200-70dB', 'noise-f0200-70dB', 'harm-f0200-65dB', 'noise-f0200-65dB'};
            otherwise
                error('Runtype is not right.');
        end
        
    case 'naturalsound_monkey'
        switch runtype
            case {'main'}
                r = varargin{optInputs(varargin, 'run')+1};
                b = read_timings(exp,us,runtype,r,'nofix');
                conds = b.conds';
            case {'main_combined','main_combined_raw'}
                conds = strrep(mydir([params('rootdir') exp '/stims/final_stims195/']),'.wav','')';
                x = regexp(conds,'(\d)*','match');
                ids = nan(size(x));
                for i = 1:length(x);
                    ids(i) = str2double(x{i}{1});
                end
                [~, xi] = sort(ids,'ascend');
                conds = conds(xi);
            case {'voice_localizer'}
                conds = {'macaques','animals','german'};
            case {'petkov_localizer','petkov_localizer_combined_split','petkov_localizer_combined_raw'}
                clear conds;
                conds{1} = 'MVocs';
                conds{2} = 'AVocs';
                conds{3} = 'EnvSounds';
                conds{4} = 'MVocsPhaseScram';
        end
        
    case 'naturalsound-nmf-localizer'
        
        switch runtype
            case 'main_v1'
                load([params('rootdir') 'naturalsound-nmf-localizer/data/experiment/original_files/data/nmf_all_angier_runNo1_cb1.mat'],'p','e');
                conds = p.conditions(2:end);
            otherwise
                error('Runtype is not right.');
        end
        
    case 'color_monkey'
        conds = {'dummy'};
        
    case 'tono-pitch-localizer'
        
        switch runtype
            
            case {'tono_pitch_localizer', 'tono_pitch_localizer_combined_split', 'tono_pitch_localizer_combined_raw'}
                
                conds = {'freq200','freq400','freq800','freq1600','freq3200','freq6400','harm-f0333','noise-diotic-f0333'};
                
        end
        
    case 'tono-localizer'
        
        switch runtype
            
            case {'tono_localizer', 'tono_localizer_combined_split', 'tono_localizer_combined_raw'}
                
                conds = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
                
        end
        
    otherwise
        error('No valid experiment');
        
end