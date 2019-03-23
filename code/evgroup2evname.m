function evname = evgroup2evname(exp,us,evgroup,model,varargin)

% function evname = evgroup2evname(exp,us,evgroup,model)
% pitched-straight-harmvoc-noisevoc-70-75dB
% used for creating groups of evs to be used in a contrast
% currently there are no groups

switch exp
    
    case 'pitch_resthr'
        
        switch evgroup
            case 'harm'
                evname = {'harm-1-2', 'harm-2-4', 'harm-3-6', 'harm-4-8', 'harm-5-10', 'harm-7-14', 'harm-9-18', 'harm-12-24', 'harm-15-30', 'harm-20-40'};
                
            case 'all'
                evname = read_conditions(us,'main');
                
            otherwise
                evname = {evgroup};
        end
        
    case {'pitch_resthr_v2','pitch_resthr_v3','pitch_resthr_v4','tono-pitch-localizer', 'tono-localizer'}
        
        switch evgroup
            
            case 'all_main'
                evname = read_conditions(exp,us,'main');
                
            case 'all_localizer'
                evname = read_conditions(exp,us,'localizer');
                
            case {'harm-all', 'harm'}
                evname = {...
                    'freq1200-2400_harm-sine3-6',  'freq1200-2400_harm-sine4-8',  'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12', 'freq1200-2400_harm-sine8-16',...
                    'freq1200-2400_harm-sine10-20','freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30',...
                    'freq2400-4800_harm-sine3-6',  'freq2400-4800_harm-sine4-8',  'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12', 'freq2400-4800_harm-sine8-16',...
                    'freq2400-4800_harm-sine10-20','freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30',...
                    };
                
            case 'harm34'
                
                evname = {...
                    'freq1200-2400_harm-sine3-6','freq1200-2400_harm-sine4-8',...
                    'freq2400-4800_harm-sine3-6','freq2400-4800_harm-sine4-8',...
                    };
                
            case 'harm1215'
                evname = {...
                    'freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30',...
                    'freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30',...
                    };
                
            case 'harm3456'
                evname = {...
                    'freq1200-2400_harm-sine3-6','freq1200-2400_harm-sine4-8',  'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12',...
                    'freq2400-4800_harm-sine3-6','freq2400-4800_harm-sine4-8',  'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12',...
                    };
                
            case 'freq1200-2400_harm3456'
                evname = {...
                    'freq1200-2400_harm-sine3-6','freq1200-2400_harm-sine4-8',  'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12',...
                    };
                
            case 'freq2400-4800_harm3456'
                evname = {...
                    'freq2400-4800_harm-sine3-6','freq2400-4800_harm-sine4-8',  'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12',...
                    };
                
            case 'harm101215'
                evname = {...
                    'freq1200-2400_harm-sine10-20','freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30',...
                    'freq2400-4800_harm-sine10-20','freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30',...
                    };
                
            case 'harm1200-2400'
                evname = {...
                    'freq1200-2400_harm-sine3-6','freq1200-2400_harm-sine4-8',  'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12', 'freq1200-2400_harm-sine8-16',...
                    'freq1200-2400_harm-sine10-20','freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30',...
                    };
                
            case 'harm2400-4800'
                evname = {...
                    'freq2400-4800_harm-sine3-6','freq2400-4800_harm-sine4-8',  'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12', 'freq2400-4800_harm-sine8-16',...
                    'freq2400-4800_harm-sine10-20','freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30',...
                    };
                
            case 'noise1200-2400'
                evname = {'freq1200-2400_noise'};
                
            case 'noise2400-4800'
                evname = {'freq2400-4800_noise'};
                
            case {'noise-all','noise'}
                evname = {'freq1200-2400_noise','freq2400-4800_noise'};
                
            case 'freq1200-2400'
                evname = {...
                    'freq1200-2400_harm-sine3-6','freq1200-2400_harm-sine4-8',  'freq1200-2400_harm-sine5-10', 'freq1200-2400_harm-sine6-12', 'freq1200-2400_harm-sine8-16',...
                    'freq1200-2400_harm-sine10-20','freq1200-2400_harm-sine12-24','freq1200-2400_harm-sine15-30','freq1200-2400_noise',...
                    };
                
            case 'freq2400-4800'
                evname = {...
                    'freq2400-4800_harm-sine3-6','freq2400-4800_harm-sine4-8',  'freq2400-4800_harm-sine5-10', 'freq2400-4800_harm-sine6-12', 'freq2400-4800_harm-sine8-16',...
                    'freq2400-4800_harm-sine10-20','freq2400-4800_harm-sine12-24','freq2400-4800_harm-sine15-30','freq2400-4800_noise',...
                    };
                
            case 'freq200-400'
                evname = {'freq200','freq800'};
                
            case 'freq200-6400'
                evname = {'freq200','freq6400'};
                
            case 'freq800-1600'
                evname = {'freq800','freq1600'};
                
            case 'freq3200-6400'
                evname = {'freq3200','freq6400'};
                
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq1600-3200-6400'
                evname = {'freq1600','freq3200','freq6400'};
                
            case 'freq800-1600-3200-6400'
                evname = {'freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-3200-6400'
                evname = {'freq200','freq400','freq3200','freq6400'};
                
            case 'freq200-400-800-1600'
                evname = {'freq200','freq400','freq800','freq1600'};
                
            case 'freq400-800-1600-3200-6400'
                evname = {'freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-800-1600-3200-6400'
                evname = {'freq200','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-1600-3200-6400'
                evname = {'freq200','freq400','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-800-3200-6400'
                evname = {'freq200','freq400','freq800','freq3200','freq6400'};
                
            case 'freq200-400-800-1600-6400'
                evname = {'freq200','freq400','freq800','freq1600','freq6400'};
                
            case 'freq200-400-800-1600-3200'
                evname = {'freq200','freq400','freq800','freq1600','freq3200'};
                
            case 'allfreq'
                evname = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'harm-sine3-6'
                evname = {'freq1200-2400_harm-sine3-6','freq2400-4800_harm-sine3-6'};
                
            case 'harm-sine4-8'
                evname = {'freq1200-2400_harm-sine4-8','freq2400-4800_harm-sine4-8'};
                
            case 'harm-sine5-10'
                evname = {'freq1200-2400_harm-sine5-10','freq2400-4800_harm-sine5-10'};
                
            case 'harm-sine6-12'
                evname = {'freq1200-2400_harm-sine6-12','freq2400-4800_harm-sine6-12'};
                
            case 'harm-sine8-16'
                evname = {'freq1200-2400_harm-sine8-16','freq2400-4800_harm-sine8-16'};
                
            case 'harm-sine10-20'
                evname = {'freq1200-2400_harm-sine10-20','freq2400-4800_harm-sine10-20'};
                
            case 'harm-sine12-24'
                evname = {'freq1200-2400_harm-sine12-24','freq2400-4800_harm-sine12-24'};
                
            case 'harm-sine15-30'
                evname = {'freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30'};
                
            case 'harm-schr-1ERB15-30'
                evname = {'freq1200-2400_harm-schr-1ERB15-30','freq2400-4800_harm-schr-1ERB15-30'};
                
            case 'harm-schr-1ERB3-6'
                evname = {'freq1200-2400_harm-schr-1ERB3-6','freq2400-4800_harm-schr-1ERB3-6'};
                
            otherwise
                evname = {evgroup};
                
        end
        
    case 'pitch_misc'
        
        switch evgroup
            
            case 'all_schroeder'
                evname = read_conditions(exp,us,'schroeder');
                
            case 'harm-schr-1ERB3-6'
                evname = {'freq1200-2400_harm-schr-1ERB3-6','freq2400-4800_harm-schr-1ERB3-6'};
                
            case 'harm-sine15-30'
                evname = {'freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30'};
                
            case 'harm-schr-1ERB15-30'
                evname = {'freq1200-2400_harm-schr-1ERB15-30','freq2400-4800_harm-schr-1ERB15-30'};
                
            case 'noise'
                evname = {'freq1200-2400_noise','freq2400-4800_noise'};
                
            otherwise
                evname = {evgroup};
                
        end
        
    case 'pitch_f02by2'
        
        switch evgroup
            
            case 'all'
                evname = read_conditions(exp,us,'f0adapt');
                
            case 'harm-f0ref-f0diff'
                evname = {'harm-f0ref-f0diff-freqsame','harm-f0ref-f0diff-freqdiff'};
                
            otherwise
                evname = {evgroup};
                
        end
        
    case 'pitch_f0smallfreq'
        
        switch evgroup
            
            case 'low'
                evname = {...
                    'harm-low-f0same-freqsame','harm-low-f0same-freqdiff', 'harm-low-f0diff-freqdiff'...
                    'inharm-low-f0same-freqsame','inharm-low-f0same-freqdiff', 'inharm-low-f0diff-freqdiff'...
                    };
                
            case 'high'
                evname = {...
                    'harm-high-f0same-freqsame','harm-high-f0same-freqdiff', 'harm-high-f0diff-freqdiff'...
                    'inharm-high-f0same-freqsame','inharm-high-f0same-freqdiff', 'inharm-high-f0diff-freqdiff'...
                    };
                
            case 'harm-f0same-freqsame'
                evname = {'harm-low-f0same-freqsame','harm-high-f0same-freqsame'};
                
            case 'harm-f0same-freqdiff'
                evname = {'harm-low-f0same-freqdiff','harm-high-f0same-freqdiff'};
                
            case 'harm-f0diff-freqdiff'
                evname = {'harm-low-f0diff-freqdiff','harm-high-f0diff-freqdiff'};
                
            case 'inharm-f0same-freqsame'
                evname = {'inharm-low-f0same-freqsame','inharm-high-f0same-freqsame'};
                
            case 'inharm-f0same-freqdiff'
                evname = {'inharm-low-f0same-freqdiff','inharm-high-f0same-freqdiff'};
                
            case 'inharm-f0diff-freqdiff'
                evname = {'inharm-low-f0diff-freqdiff','inharm-high-f0diff-freqdiff'};
                
            case 'harm-low'
                evname = {'harm-low-f0same-freqsame','harm-low-f0same-freqdiff', 'harm-low-f0diff-freqdiff'};
                
            case 'harm-high'
                evname = {'harm-high-f0same-freqsame','harm-high-f0same-freqdiff', 'harm-high-f0diff-freqdiff'};
                
            case 'inharm-low'
                evname = {'inharm-low-f0same-freqsame','inharm-low-f0same-freqdiff', 'inharm-low-f0diff-freqdiff'};
                
            case 'inharm-high'
                evname = {'inharm-high-f0same-freqsame','inharm-high-f0same-freqdiff', 'inharm-high-f0diff-freqdiff'};
                
            case 'harm'
                evname = {...
                    'harm-low-f0same-freqsame','harm-low-f0same-freqdiff', 'harm-low-f0diff-freqdiff',...
                    'harm-high-f0same-freqsame','harm-high-f0same-freqdiff', 'harm-high-f0diff-freqdiff'};
                
            case 'inharm'
                evname = {...
                    'inharm-low-f0same-freqsame','inharm-low-f0same-freqdiff', 'inharm-low-f0diff-freqdiff'...
                    'inharm-high-f0same-freqsame','inharm-high-f0same-freqdiff', 'inharm-high-f0diff-freqdiff'};
                
            case 'all'
                evname = read_conditions(exp,us,'f0adapt');
                
            otherwise
                evname = {evgroup};
                
        end
        
    case 'pitch_f0'
        
        switch evgroup
            case 'harm'
                evname = {'harm-f0125','harm-f0275','harm-f0600'};
            case 'noise'
                evname = {'noise-f0125','noise-f0275','noise-f0600'};
                
            case 'f0125'
                evname = {'harm-f0125','noise-f0125'};
            case 'f0275'
                evname = {'harm-f0275','noise-f0275'};
            case 'f0600'
                evname = {'harm-f0600','noise-f0600'};
                
            case 'f0125-f0275'
                evname = {'harm-f0125','noise-f0125','harm-f0275','noise-f0275'};
            case 'f0275-f0600'
                evname = {'harm-f0275','noise-f0275','harm-f0600','noise-f0600'};
            case 'f0125-f0600'
                evname = {'harm-f0125','noise-f0125','harm-f0600','noise-f0600'};
                
            case 'harm-f0125-f0275'
                evname = {'harm-f0125','harm-f0275'};
            case 'harm-f0275-f0600'
                evname = {'harm-f0275','harm-f0600'};
            case 'harm-f0125-f0600'
                evname = {'harm-f0125','harm-f0600'};
                
            case 'noise-f0125-f0275'
                evname = {'noise-f0125','noise-f0275'};
            case 'noise-f0275-f0600'
                evname = {'noise-f0275','noise-f0600'};
            case 'noise-f0125-f0600'
                evname = {'noise-f0125','noise-f0600'};
                
            case 'all_adapt6note'
                evname = read_conditions(exp,us,'adapt6note');
            case 'all_localizer'
                evname = read_conditions(exp,us,'localizer');
                
            case 'f0same-freqsame'
                evname = {'f0same-freqsame-AA', 'f0same-freqsame-BB', 'f0same-freqsame-CC', 'f0same-freqsame-DD'};
            case 'f0same-freqdiff'
                evname = {'f0same-freqdiff-AB', 'f0same-freqdiff-BA', 'f0same-freqdiff-CD', 'f0same-freqdiff-DC'};
            case 'f0diff-freqdiff'
                evname = {'f0diff-freqdiff-AC', 'f0diff-freqdiff-CA', 'f0diff-freqdiff-BD', 'f0diff-freqdiff-DB'};
                
            case 'f0same-freqsame-AABB'
                evname = {'f0same-freqsame-AA', 'f0same-freqsame-BB'};
            case 'f0same-freqsame-CCDD'
                evname = {'f0same-freqsame-CC', 'f0same-freqsame-DD'};
                
            case 'f0same-freqdiff-ABBA'
                evname = {'f0same-freqdiff-AB', 'f0same-freqdiff-BA'};
            case 'f0same-freqdiff-CDDC'
                evname = {'f0same-freqdiff-CD', 'f0same-freqdiff-DC'};
                
            case 'f0diff-freqdiff-ACCA'
                evname = {'f0diff-freqdiff-AC', 'f0diff-freqdiff-CA'};
            case 'f0diff-freqdiff-BDDB'
                evname = {'f0diff-freqdiff-BD', 'f0diff-freqdiff-DB'};
                
            case 'f0same-freqsame-AABBCCDD'
                evname = {'f0same-freqsame-AA', 'f0same-freqsame-BB','f0same-freqsame-CC', 'f0same-freqsame-DD'};
                
            case 'f0same-freqdiff-ABBACDDC'
                evname = {'f0same-freqdiff-AB', 'f0same-freqdiff-BA','f0same-freqdiff-CD', 'f0same-freqdiff-DC'};
                
            case 'f0diff-freqdiff-ACCABDDB'
                evname = {'f0diff-freqdiff-AC', 'f0diff-freqdiff-CA','f0diff-freqdiff-BD', 'f0diff-freqdiff-DB'};
                
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_f0inharm'
        
        switch evgroup
            
            case 'all'
                evname = read_conditions(exp,us,'f0adapt');
                
                
                
            case 'f0same-freqsame-AA'
                evname = {'harm-f0same-freqsame-AA', 'inharm-f0same-freqsame-AA'};
                
            case 'f0same-freqsame-DD'
                evname = {'harm-f0same-freqsame-DD', 'inharm-f0same-freqsame-DD'};
                
                
                
            case 'harm'
                evname = {...
                    'harm-f0same-freqsame-AA', 'harm-f0same-freqsame-BB', 'harm-f0same-freqsame-CC', 'harm-f0same-freqsame-DD',...
                    'harm-f0same-freqdiff-AB', 'harm-f0same-freqdiff-BA', 'harm-f0same-freqdiff-CD', 'harm-f0same-freqdiff-DC',...
                    'harm-f0diff-freqdiff-AC', 'harm-f0diff-freqdiff-CA', 'harm-f0diff-freqdiff-BD', 'harm-f0diff-freqdiff-DB'};
            case 'inharm'
                evname = {...
                    'inharm-f0same-freqsame-AA', 'inharm-f0same-freqsame-BB', 'inharm-f0same-freqsame-CC', 'inharm-f0same-freqsame-DD',...
                    'inharm-f0same-freqdiff-AB', 'inharm-f0same-freqdiff-BA', 'inharm-f0same-freqdiff-CD', 'inharm-f0same-freqdiff-DC',...
                    'inharm-f0diff-freqdiff-AC', 'inharm-f0diff-freqdiff-CA', 'inharm-f0diff-freqdiff-BD', 'inharm-f0diff-freqdiff-DB'};
                
                
                
            case 'harm-AA-CC-AC-CA'
                evname = {'harm-f0same-freqsame-AA', 'harm-f0same-freqsame-CC', 'harm-f0diff-freqdiff-AC', 'harm-f0diff-freqdiff-CA'};
                
            case 'inharm-AA-CC-AC-CA'
                evname = {'inharm-f0same-freqsame-AA', 'inharm-f0same-freqsame-CC', 'inharm-f0diff-freqdiff-AC', 'inharm-f0diff-freqdiff-CA'};
                
            case 'harm-BB-DD-BD-DB'
                evname = {'harm-f0same-freqsame-BB', 'harm-f0same-freqsame-DD', 'harm-f0diff-freqdiff-BD', 'harm-f0diff-freqdiff-DB'};
                
            case 'inharm-BB-DD-BD-DB'
                evname = {'inharm-f0same-freqsame-BB', 'inharm-f0same-freqsame-DD', 'inharm-f0diff-freqdiff-BD', 'inharm-f0diff-freqdiff-DB'};
                
                
                
            case 'harm-f0same-freqsame'
                evname = {'harm-f0same-freqsame-AA', 'harm-f0same-freqsame-BB', 'harm-f0same-freqsame-CC', 'harm-f0same-freqsame-DD'};
                
            case 'harm-f0same-freqdiff'
                evname = {'harm-f0same-freqdiff-AB', 'harm-f0same-freqdiff-BA', 'harm-f0same-freqdiff-CD', 'harm-f0same-freqdiff-DC'};
                
            case 'harm-f0diff-freqdiff'
                evname = {'harm-f0diff-freqdiff-AC', 'harm-f0diff-freqdiff-CA', 'harm-f0diff-freqdiff-BD', 'harm-f0diff-freqdiff-DB'};
                
            case 'inharm-f0same-freqsame'
                evname = {'inharm-f0same-freqsame-AA', 'inharm-f0same-freqsame-BB', 'inharm-f0same-freqsame-CC', 'inharm-f0same-freqsame-DD'};
                
            case 'inharm-f0same-freqdiff'
                evname = {'inharm-f0same-freqdiff-AB', 'inharm-f0same-freqdiff-BA', 'inharm-f0same-freqdiff-CD', 'inharm-f0same-freqdiff-DC'};
                
            case 'inharm-f0diff-freqdiff'
                evname = {'inharm-f0diff-freqdiff-AC', 'inharm-f0diff-freqdiff-CA', 'inharm-f0diff-freqdiff-BD', 'inharm-f0diff-freqdiff-DB'};
                
                
                
            case 'harm-f0same-freqsame-AABB'
                evname = {'harm-f0same-freqsame-AA', 'harm-f0same-freqsame-BB'};
            case 'harm-f0same-freqsame-CCDD'
                evname = {'harm-f0same-freqsame-CC', 'harm-f0same-freqsame-DD'};
                
            case 'harm-f0same-freqdiff-ABBA'
                evname = {'harm-f0same-freqdiff-AB', 'harm-f0same-freqdiff-BA'};
            case 'harm-f0same-freqdiff-CDDC'
                evname = {'harm-f0same-freqdiff-CD', 'harm-f0same-freqdiff-DC'};
                
            case 'harm-f0diff-freqdiff-ACCA'
                evname = {'harm-f0diff-freqdiff-AC', 'harm-f0diff-freqdiff-CA'};
            case 'harm-f0diff-freqdiff-BDDB'
                evname = {'harm-f0diff-freqdiff-BD', 'harm-f0diff-freqdiff-DB'};
                
                
                
            case 'inharm-f0same-freqsame-AABB'
                evname = {'inharm-f0same-freqsame-AA', 'inharm-f0same-freqsame-BB'};
            case 'inharm-f0same-freqsame-CCDD'
                evname = {'inharm-f0same-freqsame-CC', 'inharm-f0same-freqsame-DD'};
                
            case 'inharm-f0same-freqdiff-ABBA'
                evname = {'inharm-f0same-freqdiff-AB', 'inharm-f0same-freqdiff-BA'};
            case 'inharm-f0same-freqdiff-CDDC'
                evname = {'inharm-f0same-freqdiff-CD', 'inharm-f0same-freqdiff-DC'};
                
            case 'inharm-f0diff-freqdiff-ACCA'
                evname = {'inharm-f0diff-freqdiff-AC', 'inharm-f0diff-freqdiff-CA'};
            case 'inharm-f0diff-freqdiff-BDDB'
                evname = {'inharm-f0diff-freqdiff-BD', 'inharm-f0diff-freqdiff-DB'};
                
                
                
                
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_contour'
        
        switch evgroup
            case 'harm'
                evname = {'harm-f0125','harm-f0275','harm-f0600'};
            case 'noise'
                evname = {'noise-f0125','noise-f0275','noise-f0600'};
                
            case 'f0125'
                evname = {'harm-f0125','noise-f0125'};
            case 'f0275'
                evname = {'harm-f0275','noise-f0275'};
            case 'f0600'
                evname = {'harm-f0600','noise-f0600'};
                
            case 'f0125-f0275'
                evname = {'harm-f0125','noise-f0125','harm-f0275','noise-f0275'};
            case 'f0275-f0600'
                evname = {'harm-f0275','noise-f0275','harm-f0600','noise-f0600'};
            case 'f0125-f0600'
                evname = {'harm-f0125','noise-f0125','harm-f0600','noise-f0600'};
                
            case 'harm-f0125-f0275'
                evname = {'harm-f0125','harm-f0275'};
            case 'harm-f0275-f0600'
                evname = {'harm-f0275','harm-f0600'};
            case 'harm-f0125-f0600'
                evname = {'harm-f0125','harm-f0600'};
                
            case 'noise-f0125-f0275'
                evname = {'noise-f0125','noise-f0275'};
            case 'noise-f0275-f0600'
                evname = {'noise-f0275','noise-f0600'};
            case 'noise-f0125-f0600'
                evname = {'noise-f0125','noise-f0600'};
                
            case 'all_main'
                evname = {'harm-freqsame-melsame', 'harm-freqdiff-melsame', 'harm-freqsame-meldiff', 'harm-freqdiff-meldiff'};
            case 'all_localizer'
                evname = {'harm-f0125','noise-f0125','harm-f0275','noise-f0275','harm-f0600','noise-f0600'};
                
            case 'harm-meldiff'
                evname = {'harm-freqsame-meldiff','harm-freqdiff-meldiff'};
            case 'harm-melsame'
                evname = {'harm-freqsame-melsame','harm-freqdiff-melsame'};
            case 'harm-freqdiff'
                evname = {'harm-freqdiff-melsame','harm-freqdiff-meldiff'};
            case 'harm-freqsame'
                evname = {'harm-freqsame-melsame','harm-freqsame-meldiff'};
                
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_dp'
        
        switch evgroup
            case 'harm-wawomask'
                evname = {'harm-highfreq','harm-lowfreq','harm-highfreq-lowmask','harm-lowfreq-lowmask'};
            case 'noise-wawomask'
                evname = {'noise-highfreq','noise-lowfreq','noise-highfreq-lowmask','noise-lowfreq-lowmask'};
            case 'harm'
                evname = {'harm-highfreq','harm-lowfreq'};
            case 'noise'
                evname = {'noise-highfreq','noise-lowfreq'};
            case 'harm-lowmask'
                evname = {'harm-highfreq-lowmask','harm-lowfreq-lowmask'};
            case 'noise-lowmask'
                evname = {'noise-highfreq-lowmask','noise-lowfreq-lowmask'};
            case 'harm-highfreq-wawomask'
                evname = {'harm-highfreq','harm-highfreq-lowmask'};
            case 'noise-highfreq-wawomask'
                evname = {'noise-highfreq','noise-highfreq-lowmask'};
            case 'harm-lowfreq-wawomask'
                evname = {'harm-lowfreq','harm-lowfreq-lowmask'};
            case 'noise-lowfreq-wawomask'
                evname = {'noise-lowfreq','noise-lowfreq-lowmask'};
            case 'harm-unres-wawomask'
                evname = {'harm-unres','harm-unres-lowmask'};
            case 'harm-unres-schroeder-neg-wawomask'
                evname = {'harm-unres-schroeder-neg','harm-unres-schroeder-neg-lowmask'};
            case 'highfreq'
                evname = {'harm-highfreq','noise-highfreq'};
            case 'lowfreq'
                evname = {'harm-lowfreq','noise-lowfreq'};
            case 'highfreq-lowmask'
                evname = {'harm-highfreq-lowmask','noise-highfreq-lowmask'};
            case 'lowfreq-lowmask'
                evname = {'harm-lowfreq-lowmask','noise-lowfreq-lowmask'};
            case 'all'
                evname = {'harm-lowfreq','noise-lowfreq','harm-highfreq','noise-highfreq','harm-unres','harm-unres-schroeder-neg',...
                    'harm-lowfreq-lowmask','noise-lowfreq-lowmask','harm-highfreq-lowmask','noise-highfreq-lowmask','harm-unres-lowmask','harm-unres-schroeder-neg-lowmask'...
                    'harm-highfreq-schroeder-negN-1'};
            otherwise
                evname = {evgroup};
        end
        
    case 'mspec'
        
        if strcmp(model,'block');
            switch evgroup
                case 'pitch'
                    evname = {'env-pitch','songs','mel-scram'};
                case 'nopitch'
                    evname = {'env-nopitch','lyrics','mel-noise'};
                case 'speech'
                    evname = {'sentences','nonwords'};
                case 'nospeech'
                    evname = {'orchestra','solo','animals-pitch','env-pitch','env-nopitch','drums-intact','drums-scram','mel-intact','mel-scram','mel-interp','mel-noise'};
                case 'all'
                    evname = read_conditions(s);
                otherwise
                    evname = {evgroup};
            end
        end
        
    case 'sepi'
        
        switch evgroup
            case 'all'
                evname = {'pitch','noise','highfreq','lowfreq','sentences','nonwords','solo'};
            otherwise
                evname = {evgroup};
        end
        
    case 'tonotopy_rate'
        
        switch evgroup
            case 'freqs200'
                evname = {'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16'};
            case 'freqs550'
                evname = {'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16'};
            case 'freqs1500'
                evname = {'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16'};
            case 'freqs4150'
                evname = {'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16'};
                
                
                
            case 'freqs200-550'
                evname = {...
                    'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16',...
                    'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16',...
                    };
            case 'freqs550-1500'
                evname = {...
                    'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16',...
                    'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16'...
                    };
            case 'freqs1500-4150'
                evname = {...
                    'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16',...
                    'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16',...
                    };
            case 'freqs200-4150'
                evname = {...
                    'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16',...
                    'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16',...
                    };
                
                
                
            case 'freqs550-1500-4150'
                evname = {...
                    'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16',...
                    'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16',...
                    'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16',...
                    };
            case 'freqs200-1500-4150'
                evname = {...
                    'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16',...
                    'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16',...
                    'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16',...
                    };
            case 'freqs200-550-4150'
                evname = {...
                    'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16',...
                    'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16',...
                    'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16',...
                    };
            case 'freqs200-550-1500'
                evname = {...
                    'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16',...
                    'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16',...
                    'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16',...
                    };
                
                
                
            case 'rates1'
                evname = {'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1'};
            case 'rates2'
                evname = {'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2'};
            case 'rates4'
                evname = {'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4'};
            case 'rates8'
                evname = {'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8'};
            case 'rates16'
                evname = {'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16'};
                
                
                
            case 'rates1-2'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    };
            case 'rates2-4'
                evname = {...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    };
            case 'rates4-8'
                evname = {...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    };
            case 'rates8-16'
                evname = {...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
                
                
                
            case 'rates1-2-4'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    };
            case 'rates1-2-16'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
            case 'rates1-8-16'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
            case 'rates4-8-16'
                evname = {...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
                
                
                
            case 'rates2-4-8-16'
                evname = {...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
            case 'rates1-4-8-16'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
            case 'rates1-2-8-16'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
            case 'rates1-2-4-16'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    'freqs200_rates16','freqs550_rates16','freqs1500_rates16','freqs4150_rates16',...
                    };
            case 'rates1-2-4-8'
                evname = {...
                    'freqs200_rates1','freqs550_rates1','freqs1500_rates1','freqs4150_rates1',...
                    'freqs200_rates2','freqs550_rates2','freqs1500_rates2','freqs4150_rates2',...
                    'freqs200_rates4','freqs550_rates4','freqs1500_rates4','freqs4150_rates4',...
                    'freqs200_rates8','freqs550_rates8','freqs1500_rates8','freqs4150_rates8',...
                    };
                
            case 'all'
                evname = {...
                    'freqs200_rates1','freqs200_rates2','freqs200_rates4','freqs200_rates8','freqs200_rates16',...
                    'freqs550_rates1','freqs550_rates2','freqs550_rates4','freqs550_rates8','freqs550_rates16',...
                    'freqs1500_rates1','freqs1500_rates2','freqs1500_rates4','freqs1500_rates8','freqs1500_rates16',...
                    'freqs4150_rates1','freqs4150_rates2','freqs4150_rates4','freqs4150_rates8','freqs4150_rates16',...
                    };
                
            otherwise
                evname = {evgroup};
        end
        
    case 'amusia'
        
        switch evgroup
            
            case 'harm'
                evname = {'harm-f0167','harm-f0667'};
                
            case 'noise'
                evname = {'noise-f0167','noise-f0667'};
                
            case 'highfreq'
                evname = {'harm-f0667','noise-f0667'};
                
            case 'lowfreq'
                evname = {'harm-f0167','noise-f0167'};
                
            case 'all_localizer'
                evname = {'harm-f0167','harm-f0667','noise-f0167','noise-f0667'};
                
            otherwise
                evname = {evgroup};
                
        end
        
    case 'music_12channel'
        
        switch evgroup
            case 'solo-long-allfreq'
                evname = {'solo-long','solo-long-highfreq'};
            case 'solo'
                evname = {'solo-long','solo-long-highfreq','solo-short'};
            case 'pitch-nomusic'
                evname = {'env-pitch','midi-onsetintervals-scram-jitter'};
            otherwise
                evname = {evgroup};
                
        end
        
    case {'stax','multiband'}
        switch evgroup
            case 'harm'
                evname = {'harm-f0125','harm-f0600'};
            case 'noise'
                evname = {'noise-f0125','noise-f0600'};
            case 'f0125'
                evname = {'harm-f0125','noise-f0125'};
            case 'f0600'
                evname = {'harm-f0600','noise-f0600'};
            case 'all_localizer'
                evname = {'harm-f0125','harm-f0600','noise-f0125','noise-f0600'};
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_dp_v2'
        switch evgroup
            case 'harm-sine-3-6-wawomask'
                evname = {'freq2000-4000_harm-sine-3-6','freq2000-4000_harm-sine-masknoise-3-6'};
            case 'harm-sine-wawomask'
                evname = {'freq2000-4000_harm-sine-3-6','freq2000-4000_harm-sine-masknoise-3-6','freq2000-4000_harm-sine-15-30','freq2000-4000_harm-sine-masknoise-15-30'};
            case 'noise-wawomask'
                evname = {'freq2000-4000_noise','freq2000-4000_noise-masknoise'};
            case 'all_main'
                evname = read_conditions(exp,us,'main');
            otherwise
                evname = {evgroup};
        end
        
    case 'music_scram'
        switch evgroup
            case 'music-intact'
                evname = {'midi-intact-10sec','bigband-intact-10sec'};
            case 'music-quilt'
                evname = {'midi-quilt-30ms-10sec','bigband-quilt-30ms-10sec'};
            case 'music-scram-quilt'
                evname = {'midi-onsetintervals-scram-jitter-10sec','midi-quilt-30ms-10sec','bigband-quilt-30ms-10sec'};
            case 'all_main'
                evname = read_conditions(exp,us,'main');
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_f0adapt'
        switch evgroup
            case 'freq200-400'
                evname = {'freq200','freq400'};
                
            case 'freq3200-6400'
                evname = {'freq3200','freq6400'};
                
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq1600-3200-6400'
                evname = {'freq1600','freq3200','freq6400'};
                
            case 'allfreq'
                evname = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'harm_noise'
                evname = {'harm-f0333','noise-f0333'};
                
            case 'all_localizer'
                evname = read_conditions(exp,us,'localizer',varargin{:});
                
            case 'f0diff'
                evname = {...
                    'harm-0-f0diff-freqdiff',...
                    'harm-10-f0diff-freqdiff',...
                    'harm-inf-f0diff-freqdiff',...
                    'inharm-0-f0diff-freqdiff',...
                    'inharm-10-f0diff-freqdiff',...
                    'inharm-inf-f0diff-freqdiff',...
                    };
                
            case 'f0same'
                evname = {...
                    'harm-0-f0same-freqdiff',...
                    'harm-10-f0same-freqdiff',...
                    'harm-inf-f0same-freqdiff',...
                    'inharm-0-f0same-freqdiff',...
                    'inharm-10-f0same-freqdiff',...
                    'inharm-inf-f0same-freqdiff',...
                    };
                
            case 'f0diff-novariation'
                evname = {...
                    'harm-0-f0diff-freqdiff',...
                    'inharm-0-f0diff-freqdiff',...
                    };
                
            case 'f0same-novariation'
                evname = {...
                    'harm-0-f0same-freqdiff',...
                    'inharm-0-f0same-freqdiff',...
                    };
                
            case 'harm'
                evname = {...
                    'harm-0-f0same-freqdiff','harm-0-f0diff-freqdiff',...
                    'harm-10-f0same-freqdiff','harm-10-f0diff-freqdiff',...
                    'harm-inf-f0same-freqdiff','harm-inf-f0diff-freqdiff',...
                    };
                
            case 'inharm'
                evname = {...
                    'inharm-0-f0same-freqdiff','inharm-0-f0diff-freqdiff',...
                    'inharm-10-f0same-freqdiff','inharm-10-f0diff-freqdiff',...
                    'inharm-inf-f0same-freqdiff','inharm-inf-f0diff-freqdiff',...
                    };
                
            case 'f0diff-harm'
                evname = {...
                    'harm-0-f0diff-freqdiff',...
                    'harm-10-f0diff-freqdiff',...
                    'harm-inf-f0diff-freqdiff',...
                    };
                
            case 'f0same-harm'
                evname = {...
                    'harm-0-f0same-freqdiff',...
                    'harm-10-f0same-freqdiff',...
                    'harm-inf-f0same-freqdiff',...
                    };
                
            case 'f0diff-inharm'
                evname = {...
                    'inharm-0-f0diff-freqdiff',...
                    'inharm-10-f0diff-freqdiff',...
                    'inharm-inf-f0diff-freqdiff',...
                    };
                
            case 'f0same-inharm'
                evname = {...
                    'inharm-0-f0same-freqdiff',...
                    'inharm-10-f0same-freqdiff',...
                    'inharm-inf-f0same-freqdiff',...
                    };
                
            case 'all_f0adapt'
                evname = read_conditions(exp,us,'f0adapt',varargin{:});
                
            otherwise
                evname = {evgroup};
        end
        
    case 'music_scram_familiar'
        switch evgroup
            case 'familiar-style-intact'
                evname = {'familiar-music-intact-10sec','unfamiliar-music-intact-10sec'};
            case 'familiar-style-quilt'
                evname = {'familiar-music-quilt-30ms-10sec','unfamiliar-music-quilt-30ms-10sec'};
            case 'music-intact'
                evname = {'familiar-music-intact-10sec','unfamiliar-music-intact-10sec','midi-intact-10sec','bigband-intact-10sec'};
            case 'music-quilt'
                evname = {'familiar-music-quilt-30ms-10sec','unfamiliar-music-quilt-30ms-10sec','midi-quilt-30ms-10sec','bigband-quilt-30ms-10sec'};
            case 'music-scram-quilt'
                evname = {'familiar-music-quilt-30ms-10sec','unfamiliar-music-quilt-30ms-10sec','midi-onsetintervals-scram-jitter-10sec','midi-quilt-30ms-10sec','bigband-quilt-30ms-10sec'};
            case 'all_main'
                evname = read_conditions(exp,us,'main');
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_overlap'
        switch evgroup
            case 'manynote-harm3-6'
                evname = {'freq1200-2400_harm-sine3-6','freq2400-4800_harm-sine3-6'};
            case 'manynote-harm4-8'
                evname = {'freq1200-2400_harm-sine4-8','freq2400-4800_harm-sine4-8'};
            case 'manynote-harm15-30'
                evname = {'freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30'};
            case 'manynote-noise'
                evname = {'freq1200-2400_noise','freq2400-4800_noise'};
            case 'manynote'
                if any(us == [3 84])
                    evname = {'freq1200-2400_harm-sine3-6','freq2400-4800_harm-sine3-6','freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30','freq1200-2400_noise','freq2400-4800_noise'};
                else
                    evname = {'freq1200-2400_harm-sine4-8','freq2400-4800_harm-sine4-8','freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30','freq1200-2400_noise','freq2400-4800_noise'};
                end
            case '2note'
                if any(us == [3 84])
                    evname = {'2note-lowhigh_harm3-6','2note-highlow_harm3-6','2note-lowhigh_harm15-30','2note-highlow_harm15-30','2note-lowhigh_noise','2note-highlow_noise'};
                else
                    evname = {'2note-lowhigh_harm4-8','2note-highlow_harm4-8','2note-lowhigh_harm15-30','2note-highlow_harm15-30','2note-lowhigh_noise','2note-highlow_noise'};
                end
            case '1note'
                if any(us == [3 84])
                    evname = {'1note-lowfreq_harm3-6','1note-highfreq_harm3-6','1note-lowfreq_harm15-30','1note-highfreq_harm15-30','1note-lowfreq_noise','1note-highfreq_noise'};
                else
                    evname = {'1note-lowfreq_harm4-8','1note-highfreq_harm4-8','1note-lowfreq_harm15-30','1note-highfreq_harm15-30','1note-lowfreq_noise','1note-highfreq_noise'};
                end
            case '2note-harm3-6'
                evname = {'2note-lowhigh_harm3-6','2note-highlow_harm3-6'};
            case '2note-harm4-8'
                evname = {'2note-lowhigh_harm4-8','2note-highlow_harm4-8'};
            case '2note-harm15-30'
                evname = {'2note-lowhigh_harm15-30','2note-highlow_harm15-30'};
            case '2note-noise'
                evname = {'2note-lowhigh_noise','2note-highlow_noise'};
            case '1note-harm3-6'
                evname = {'1note-lowfreq_harm3-6','1note-highfreq_harm3-6'};
            case '1note-harm4-8'
                evname = {'1note-lowfreq_harm4-8','1note-highfreq_harm4-8'};
            case '1note-harm15-30'
                evname = {'1note-lowfreq_harm15-30','1note-highfreq_harm15-30'};
            case '1note-noise'
                evname = {'1note-lowfreq_noise','1note-highfreq_noise'};
            case 'harm3-6'
                evname = {'freq1200-2400_harm-sine3-6','freq2400-4800_harm-sine3-6','2note-lowhigh_harm3-6','2note-highlow_harm3-6','1note-lowfreq_harm3-6','1note-highfreq_harm3-6'};
            case 'harm4-8'
                evname = {'freq1200-2400_harm-sine4-8','freq2400-4800_harm-sine4-8','2note-lowhigh_harm4-8','2note-highlow_harm4-8','1note-lowfreq_harm4-8','1note-highfreq_harm4-8'};
            case 'harm15-30'
                evname = {'freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30','2note-lowhigh_harm15-30','2note-highlow_harm15-30','1note-lowfreq_harm15-30','1note-highfreq_harm15-30'};
            case 'noise'
                evname = {'freq1200-2400_noise','freq2400-4800_noise','2note-lowhigh_noise','2note-highlow_noise','1note-lowfreq_noise','1note-highfreq_noise'};
            case 'freq200-400'
                evname = {'freq200','freq400'};
                
            case 'freq800-1600'
                evname = {'freq800','freq1600'};
                
            case 'freq3200-6400'
                evname = {'freq3200','freq6400'};
                
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq1600-3200-6400'
                evname = {'freq1600','freq3200','freq6400'};
                
            case 'freq800-1600-3200-6400'
                evname = {'freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-3200-6400'
                evname = {'freq200','freq400','freq3200','freq6400'};
                
            case 'freq200-400-800-1600'
                evname = {'freq200','freq400','freq800','freq1600'};
                
            case 'freq400-800-1600-3200-6400'
                evname = {'freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-800-1600-3200-6400'
                evname = {'freq200','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-1600-3200-6400'
                evname = {'freq200','freq400','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-800-3200-6400'
                evname = {'freq200','freq400','freq800','freq3200','freq6400'};
                
            case 'freq200-400-800-1600-6400'
                evname = {'freq200','freq400','freq800','freq1600','freq6400'};
                
            case 'freq200-400-800-1600-3200'
                evname = {'freq200','freq400','freq800','freq1600','freq3200'};
                
            case 'all_tonotopy'
                evname = read_conditions(exp,us,'tonotopy');
            case 'all_overlap'
                evname = read_conditions(exp,us,'overlap');
            case 'all_overlap_v2'
                evname = read_conditions(exp,us,'overlap_v2');
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_overlap_v2'
        switch evgroup
            case 'manynote-harm4-8'
                evname = {'freq1200-2400_harm-sine4-8','freq2400-4800_harm-sine4-8'};
            case 'manynote-harm15-30'
                evname = {'freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30'};
            case 'manynote-noise'
                evname = {'freq1200-2400_noise','freq2400-4800_noise'};
                
                
            case 'manynote'
                evname = {'freq1200-2400_harm-sine4-8','freq2400-4800_harm-sine4-8','freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30','freq1200-2400_noise','freq2400-4800_noise'};
            case '2note'
                evname = {...
                    '2note-octave-lowhigh_harm4-8','2note-octave-highlow_harm4-8','2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-octave-lowhigh_noise','2note-octave-highlow_noise',...
                    '2note-overlap-low_harm4-8','2note-overlap-high_harm4-8','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30','2note-overlap-low_noise','2note-overlap-high_noise',...
                    };
            case '1note'
                evname = {'1note-low_harm4-8','1note-high_harm4-8','1note-low_harm15-30','1note-high_harm15-30','1note-low_noise','1note-high_noise'};
            case '2note-octave'
                evname = {...
                    '2note-octave-lowhigh_harm4-8','2note-octave-highlow_harm4-8','2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-octave-lowhigh_noise','2note-octave-highlow_noise',...
                    };
            case '2note-overlap'
                evname = {...
                    '2note-overlap-low_harm4-8','2note-overlap-high_harm4-8','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30','2note-overlap-low_noise','2note-overlap-high_noise',...
                    };
                
                
            case '2note-harm4-8'
                evname = {'2note-octave-lowhigh_harm4-8','2note-octave-highlow_harm4-8','2note-overlap-low_harm4-8','2note-overlap-high_harm4-8'};
            case '2note-harm15-30'
                evname = {'2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30'};
            case '2note-noise'
                evname = {'2note-octave-lowhigh_noise','2note-octave-highlow_noise','2note-overlap-low_noise','2note-overlap-high_noise'};
                
                
            case '2note-octave-harm4-8'
                evname = {'2note-octave-lowhigh_harm4-8','2note-octave-highlow_harm4-8'};
            case '2note-octave-harm15-30'
                evname = {'2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30'};
            case '2note-octave-noise'
                evname = {'2note-octave-lowhigh_noise','2note-octave-highlow_noise'};
                
                
            case '2note-overlap-harm4-8'
                evname = {'2note-overlap-low_harm4-8','2note-overlap-high_harm4-8'};
            case '2note-overlap-harm15-30'
                evname = {'2note-overlap-low_harm15-30','2note-overlap-high_harm15-30'};
            case '2note-overlap-noise'
                evname = {'2note-overlap-low_noise','2note-overlap-high_noise'};
                
                
            case '1note-harm4-8'
                evname = {'1note-low_harm4-8','1note-high_harm4-8'};
            case '1note-harm15-30'
                evname = {'1note-low_harm15-30','1note-high_harm15-30'};
            case '1note-noise'
                evname = {'1note-low_noise','1note-high_noise'};
                
            case 'harm4-8'
                evname = {'freq1200-2400_harm-sine4-8','freq2400-4800_harm-sine4-8','2note-octave-lowhigh_harm4-8','2note-octave-highlow_harm4-8','2note-overlap-low_harm4-8','2note-overlap-high_harm4-8','1note-low_harm4-8','1note-high_harm4-8'};
            case 'harm15-30'
                evname = {'freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30','2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30','1note-low_harm15-30','1note-high_harm15-30'};
            case 'noise'
                evname = {'freq1200-2400_noise','freq2400-4800_noise','2note-octave-lowhigh_noise','2note-octave-highlow_noise','2note-overlap-low_noise','2note-overlap-high_noise','1note-low_noise','1note-high_noise'};
                
                
                
            case 'freq200-400'
                evname = {'freq200','freq400'};
                
            case 'freq800-1600'
                evname = {'freq800','freq1600'};
                
            case 'freq3200-6400'
                evname = {'freq3200','freq6400'};
                
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq1600-3200-6400'
                evname = {'freq1600','freq3200','freq6400'};
                
            case 'freq800-1600-3200-6400'
                evname = {'freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-3200-6400'
                evname = {'freq200','freq400','freq3200','freq6400'};
                
            case 'freq200-400-800-1600'
                evname = {'freq200','freq400','freq800','freq1600'};
                
            case 'freq400-800-1600-3200-6400'
                evname = {'freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-800-1600-3200-6400'
                evname = {'freq200','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-1600-3200-6400'
                evname = {'freq200','freq400','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-800-3200-6400'
                evname = {'freq200','freq400','freq800','freq3200','freq6400'};
                
            case 'freq200-400-800-1600-6400'
                evname = {'freq200','freq400','freq800','freq1600','freq6400'};
                
            case 'freq200-400-800-1600-3200'
                evname = {'freq200','freq400','freq800','freq1600','freq3200'};
                
            case 'all_tonotopy'
                evname = read_conditions(exp,us,'tonotopy');
            case 'all_overlap'
                evname = read_conditions(exp,us,'overlap');
            case 'all_overlap_v2'
                evname = read_conditions(exp,us,'overlap_v2');
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_overlap_v3'
        switch evgroup
            case 'manynote-harm3-6'
                evname = {'freq1200-2400_harm-sine3-6','freq2400-4800_harm-sine3-6'};
            case 'manynote-harm15-30'
                evname = {'freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30'};
            case 'manynote-noise'
                evname = {'freq1200-2400_noise','freq2400-4800_noise'};
                
                
            case 'manynote'
                evname = {'freq1200-2400_harm-sine3-6','freq2400-4800_harm-sine3-6','freq1200-2400_harm-sine15-30','freq2400-4800_harm-sine15-30','freq1200-2400_noise','freq2400-4800_noise'};
            case '2note'
                evname = {...
                    '2note-octave-lowhigh_harm3-6','2note-octave-highlow_harm3-6','2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-octave-lowhigh_noise','2note-octave-highlow_noise',...
                    '2note-overlap-low_harm3-6','2note-overlap-high_harm3-6','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30','2note-overlap-low_noise','2note-overlap-high_noise',...
                    };
            case '1note'
                evname = {'1note-low_harm3-6','1note-high_harm3-6','1note-low_harm15-30','1note-high_harm15-30','1note-low_noise','1note-high_noise'};
            case '2note-octave'
                evname = {...
                    '2note-octave-lowhigh_harm3-6','2note-octave-highlow_harm3-6','2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-octave-lowhigh_noise','2note-octave-highlow_noise',...
                    };
            case '2note-overlap'
                evname = {...
                    '2note-overlap-low_harm3-6','2note-overlap-high_harm3-6','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30','2note-overlap-low_noise','2note-overlap-high_noise',...
                    };
                
                
            case '2note-harm3-6'
                evname = {'2note-octave-lowhigh_harm3-6','2note-octave-highlow_harm3-6','2note-overlap-low_harm3-6','2note-overlap-high_harm3-6'};
            case '2note-harm15-30'
                evname = {'2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30'};
            case '2note-noise'
                evname = {'2note-octave-lowhigh_noise','2note-octave-highlow_noise','2note-overlap-low_noise','2note-overlap-high_noise'};
                
                
            case '2note-octave-harm3-6'
                evname = {'2note-octave-lowhigh_harm3-6','2note-octave-highlow_harm3-6'};
            case '2note-octave-harm15-30'
                evname = {'2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30'};
            case '2note-octave-noise'
                evname = {'2note-octave-lowhigh_noise','2note-octave-highlow_noise'};
                
                
            case '2note-overlap-harm3-6'
                evname = {'2note-overlap-low_harm3-6','2note-overlap-high_harm3-6'};
            case '2note-overlap-harm15-30'
                evname = {'2note-overlap-low_harm15-30','2note-overlap-high_harm15-30'};
            case '2note-overlap-noise'
                evname = {'2note-overlap-low_noise','2note-overlap-high_noise'};
                
                
            case '1note-harm3-6'
                evname = {'1note-low_harm3-6','1note-high_harm3-6'};
            case '1note-harm15-30'
                evname = {'1note-low_harm15-30','1note-high_harm15-30'};
            case '1note-noise'
                evname = {'1note-low_noise','1note-high_noise'};
                
            case 'harm3-6'
                evname = {'2note-octave-lowhigh_harm3-6','2note-octave-highlow_harm3-6','2note-overlap-low_harm3-6','2note-overlap-high_harm3-6','1note-low_harm3-6','1note-high_harm3-6'};
            case 'harm15-30'
                evname = {'2note-octave-lowhigh_harm15-30','2note-octave-highlow_harm15-30','2note-overlap-low_harm15-30','2note-overlap-high_harm15-30','1note-low_harm15-30','1note-high_harm15-30'};
            case 'noise'
                evname = {'2note-octave-lowhigh_noise','2note-octave-highlow_noise','2note-overlap-low_noise','2note-overlap-high_noise','1note-low_noise','1note-high_noise'};
                
                
                
            case 'freq200-400'
                evname = {'freq200','freq400'};
                
            case 'freq800-1600'
                evname = {'freq800','freq1600'};
                
            case 'freq3200-6400'
                evname = {'freq3200','freq6400'};
                
            case 'freq1600-3200'
                evname = {'freq1600','freq3200'};
                
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq1600-3200-6400'
                evname = {'freq1600','freq3200','freq6400'};
                
            case 'freq800-1600-3200-6400'
                evname = {'freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-3200-6400'
                evname = {'freq200','freq400','freq3200','freq6400'};
                
            case 'freq200-400-800-1600'
                evname = {'freq200','freq400','freq800','freq1600'};
                
            case 'freq400-800-1600-3200-6400'
                evname = {'freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-800-1600-3200-6400'
                evname = {'freq200','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-1600-3200-6400'
                evname = {'freq200','freq400','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-800-3200-6400'
                evname = {'freq200','freq400','freq800','freq3200','freq6400'};
                
            case 'freq200-400-800-1600-6400'
                evname = {'freq200','freq400','freq800','freq1600','freq6400'};
                
            case 'freq200-400-800-1600-3200'
                evname = {'freq200','freq400','freq800','freq1600','freq3200'};
                
            case 'all_tonotopy'
                evname = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
            case 'all_overlap_v3'
                evname = read_conditions(exp,us,'overlap_v3');
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_f0adapt_v2'
        switch evgroup
            
            case 'all'
                evname = read_conditions(exp,us,'f0adapt');
                
            case 'harm-inf-f0diff-freqdiff'
                evname = {'f0200-harm-inf-f0diff-freqdiff','f0400-harm-inf-f0diff-freqdiff','f0800-harm-inf-f0diff-freqdiff'};
                
            case 'harm-inf-f0same-freqdiff'
                evname = {'f0200-harm-inf-f0same-freqdiff','f0400-harm-inf-f0same-freqdiff','f0800-harm-inf-f0same-freqdiff'};
                
            case 'harm-inf-f0same-freqsame'
                evname = {'f0200-harm-inf-f0same-freqsame','f0400-harm-inf-f0same-freqsame','f0800-harm-inf-f0same-freqsame'};
                
            case 'inharm-inf-f0diff-freqdiff'
                evname = {'f0200-inharm-inf-f0diff-freqdiff','f0400-inharm-inf-f0diff-freqdiff','f0800-inharm-inf-f0diff-freqdiff'};
                
            case 'inharm-inf-f0same-freqdiff'
                evname = {'f0200-inharm-inf-f0same-freqdiff','f0400-inharm-inf-f0same-freqdiff','f0800-inharm-inf-f0same-freqdiff'};
                
            case 'inharm-inf-f0same-freqsame'
                evname = {'f0200-inharm-inf-f0same-freqsame','f0400-inharm-inf-f0same-freqsame','f0800-inharm-inf-f0same-freqsame'};
                
            case 'inf-f0diff-freqdiff'
                evname = {...
                    'f0200-harm-inf-f0diff-freqdiff','f0400-harm-inf-f0diff-freqdiff','f0800-harm-inf-f0diff-freqdiff',...
                    'f0200-inharm-inf-f0diff-freqdiff','f0400-inharm-inf-f0diff-freqdiff','f0800-inharm-inf-f0diff-freqdiff'};
                
            case 'inf-f0same-freqdiff'
                evname = {...
                    'f0200-harm-inf-f0same-freqdiff','f0400-harm-inf-f0same-freqdiff','f0800-harm-inf-f0same-freqdiff',...
                    'f0200-inharm-inf-f0same-freqdiff','f0400-inharm-inf-f0same-freqdiff','f0800-inharm-inf-f0same-freqdiff'};
                
            case 'inf-f0same-freqsame'
                evname = {...
                    'f0200-harm-inf-f0same-freqsame','f0400-harm-inf-f0same-freqsame','f0800-harm-inf-f0same-freqsame',...
                    'f0200-inharm-inf-f0same-freqsame','f0400-inharm-inf-f0same-freqsame','f0800-inharm-inf-f0same-freqsame'};
                
            case 'f0800-inf-f0same-freqdiff'
                evname = {'f0800-harm-inf-f0same-freqdiff','f0800-inharm-inf-f0same-freqdiff'};
                
            case 'f0800-inf-f0same-freqsame'
                evname = {'f0800-harm-inf-f0same-freqsame','f0800-inharm-inf-f0same-freqsame'};
                
            case 'f0400-inf-f0same-freqdiff'
                evname = {'f0400-harm-inf-f0same-freqdiff','f0400-inharm-inf-f0same-freqdiff'};
                
            case 'f0400-inf-f0same-freqsame'
                evname = {'f0400-harm-inf-f0same-freqsame','f0400-inharm-inf-f0same-freqsame'};
                
            case 'f0200-inf-f0same-freqdiff'
                evname = {'f0200-harm-inf-f0same-freqdiff','f0200-inharm-inf-f0same-freqdiff'};
                
            case 'f0200-inf-f0same-freqsame'
                evname = {'f0200-harm-inf-f0same-freqsame','f0200-inharm-inf-f0same-freqsame'};
                
            case 'f0800'
                evname = {...
                    'f0800-harm-inf-f0same-freqsame','f0800-harm-inf-f0same-freqdiff','f0800-harm-inf-f0diff-freqdiff',...
                    'f0800-inharm-inf-f0same-freqsame','f0800-inharm-inf-f0same-freqdiff','f0800-inharm-inf-f0diff-freqdiff'};
                
            case 'f0400'
                evname = {...
                    'f0400-harm-inf-f0same-freqsame','f0400-harm-inf-f0same-freqdiff','f0400-harm-inf-f0diff-freqdiff',...
                    'f0400-inharm-inf-f0same-freqsame','f0400-inharm-inf-f0same-freqdiff','f0400-inharm-inf-f0diff-freqdiff'};
                
            case 'f0200'
                evname = {...
                    'f0200-harm-inf-f0same-freqsame','f0200-harm-inf-f0same-freqdiff','f0200-harm-inf-f0diff-freqdiff',...
                    'f0200-inharm-inf-f0same-freqsame','f0200-inharm-inf-f0same-freqdiff','f0200-inharm-inf-f0diff-freqdiff'};
                
            case 'f0200-400'
                evname = {...
                    'f0200-harm-inf-f0same-freqsame','f0200-harm-inf-f0same-freqdiff','f0200-harm-inf-f0diff-freqdiff',...
                    'f0200-inharm-inf-f0same-freqsame','f0200-inharm-inf-f0same-freqdiff','f0200-inharm-inf-f0diff-freqdiff',...
                    'f0400-harm-inf-f0same-freqsame','f0400-harm-inf-f0same-freqdiff','f0400-harm-inf-f0diff-freqdiff',...
                    'f0400-inharm-inf-f0same-freqsame','f0400-inharm-inf-f0same-freqdiff','f0400-inharm-inf-f0diff-freqdiff'};
                
            case 'f0400-800'
                evname = {...
                    'f0400-harm-inf-f0same-freqsame','f0400-harm-inf-f0same-freqdiff','f0400-harm-inf-f0diff-freqdiff',...
                    'f0400-inharm-inf-f0same-freqsame','f0400-inharm-inf-f0same-freqdiff','f0400-inharm-inf-f0diff-freqdiff',...
                    'f0800-harm-inf-f0same-freqsame','f0800-harm-inf-f0same-freqdiff','f0800-harm-inf-f0diff-freqdiff',...
                    'f0800-inharm-inf-f0same-freqsame','f0800-inharm-inf-f0same-freqdiff','f0800-inharm-inf-f0diff-freqdiff'};
                
            case 'f0200-800'
                evname = {...
                    'f0200-harm-inf-f0same-freqsame','f0200-harm-inf-f0same-freqdiff','f0200-harm-inf-f0diff-freqdiff',...
                    'f0200-inharm-inf-f0same-freqsame','f0200-inharm-inf-f0same-freqdiff','f0200-inharm-inf-f0diff-freqdiff',...
                    'f0800-harm-inf-f0same-freqsame','f0800-harm-inf-f0same-freqdiff','f0800-harm-inf-f0diff-freqdiff',...
                    'f0800-inharm-inf-f0same-freqsame','f0800-inharm-inf-f0same-freqdiff','f0800-inharm-inf-f0diff-freqdiff'};
                
            otherwise
                evname = {evgroup};
                
        end
        
    case 'pitch_adapt_params'
        switch evgroup
            case 'freq3200'
                evname = {...
                    'freq3200_freqdiff0_rate4','freq3200_freqdiff0_rate8',...
                    'freq3200_freqdiff275_rate4','freq3200_freqdiff275_rate8',...
                    'freq3200_freqdiff550_rate4','freq3200_freqdiff550_rate8',...
                    'freq3200_freqdiff1100_rate4','freq3200_freqdiff1100_rate8',...
                    };
            case 'freq400'
                evname = {...
                    'freq400_freqdiff0_rate4','freq400_freqdiff0_rate8',...
                    'freq400_freqdiff275_rate4','freq400_freqdiff275_rate8',...
                    'freq400_freqdiff550_rate4','freq400_freqdiff550_rate8',...
                    'freq400_freqdiff1100_rate4','freq400_freqdiff1100_rate8',...
                    };
            case 'rate4'
                evname = {...
                    'freq3200_freqdiff0_rate4',...
                    'freq3200_freqdiff275_rate4',...
                    'freq3200_freqdiff550_rate4',...
                    'freq3200_freqdiff1100_rate4',...
                    'freq400_freqdiff0_rate4',...
                    'freq400_freqdiff275_rate4',...
                    'freq400_freqdiff550_rate4',...
                    'freq400_freqdiff1100_rate4',...
                    };
            case 'rate8'
                evname = {...
                    'freq3200_freqdiff0_rate8',...
                    'freq3200_freqdiff275_rate8',...
                    'freq3200_freqdiff550_rate8',...
                    'freq3200_freqdiff1100_rate8',...
                    'freq400_freqdiff0_rate8',...
                    'freq400_freqdiff275_rate8',...
                    'freq400_freqdiff550_rate8',...
                    'freq400_freqdiff1100_rate8',...
                    };
            case 'freqdiff-0'
                evname = {...
                    'freq400_freqdiff0_rate4','freq400_freqdiff0_rate8',...
                    'freq3200_freqdiff0_rate4','freq3200_freqdiff0_rate8',...
                    };
            case 'freqdiff-275-550-1100'
                evname = {...
                    'freq400_freqdiff275_rate4','freq400_freqdiff275_rate8',...
                    'freq3200_freqdiff275_rate4','freq3200_freqdiff275_rate8',...
                    'freq400_freqdiff550_rate4','freq400_freqdiff550_rate8',...
                    'freq3200_freqdiff550_rate4','freq3200_freqdiff550_rate8',...
                    'freq400_freqdiff1100_rate4','freq400_freqdiff1100_rate8',...
                    'freq3200_freqdiff1100_rate4','freq3200_freqdiff1100_rate8',...
                    };
            case 'freqdiff-1100'
                evname = {...
                    'freq400_freqdiff1100_rate4','freq400_freqdiff1100_rate8',...
                    'freq3200_freqdiff1100_rate4','freq3200_freqdiff1100_rate8',...
                    };
            case 'freqdiff-275'
                evname = {...
                    'freq400_freqdiff275_rate4','freq400_freqdiff275_rate8',...
                    'freq3200_freqdiff275_rate4','freq3200_freqdiff275_rate8',...
                    };
            case 'freqdiff-550'
                evname = {...
                    'freq400_freqdiff550_rate4','freq400_freqdiff550_rate8',...
                    'freq3200_freqdiff550_rate4','freq3200_freqdiff550_rate8',...
                    };
            case 'freq400_freqdiff0'
                evname = {...
                    'freq400_freqdiff0_rate4','freq400_freqdiff0_rate8',...
                    };
            case 'freq400_freqdiff275'
                evname = {...
                    'freq400_freqdiff275_rate4','freq400_freqdiff275_rate8',...
                    };
            case 'freq400_freqdiff550'
                evname = {...
                    'freq400_freqdiff550_rate4','freq400_freqdiff550_rate8',...
                    };
            case 'freq400_freqdiff1100'
                evname = {...
                    'freq400_freqdiff1100_rate4','freq400_freqdiff1100_rate8',...
                    };
            case 'freq3200_freqdiff0'
                evname = {...
                    'freq3200_freqdiff0_rate4','freq3200_freqdiff0_rate8',...
                    };
            case 'freq3200_freqdiff275'
                evname = {...
                    'freq3200_freqdiff275_rate4','freq3200_freqdiff275_rate8',...
                    };
            case 'freq3200_freqdiff550'
                evname = {...
                    'freq3200_freqdiff550_rate4','freq3200_freqdiff550_rate8',...
                    };
            case 'freq3200_freqdiff1100'
                evname = {...
                    'freq3200_freqdiff1100_rate4','freq3200_freqdiff1100_rate8',...
                    };
                
            case 'freqdiff0_rate4'
                evname = {...
                    'freq400_freqdiff0_rate4','freq3200_freqdiff0_rate4',...
                    };
            case 'freqdiff275_rate4'
                evname = {...
                    'freq400_freqdiff275_rate4','freq3200_freqdiff275_rate4',...
                    };
            case 'freqdiff550_rate4'
                evname = {...
                    'freq400_freqdiff550_rate4','freq3200_freqdiff550_rate4',...
                    };
            case 'freqdiff1100_rate4'
                evname = {...
                    'freq400_freqdiff1100_rate4','freq3200_freqdiff1100_rate4',...
                    };
                
            case 'freqdiff0_rate8'
                evname = {...
                    'freq400_freqdiff0_rate8','freq3200_freqdiff0_rate8',...
                    };
            case 'freqdiff275_rate8'
                evname = {...
                    'freq400_freqdiff275_rate8','freq3200_freqdiff275_rate8',...
                    };
            case 'freqdiff550_rate8'
                evname = {...
                    'freq400_freqdiff550_rate8','freq3200_freqdiff550_rate8',...
                    };
            case 'freqdiff1100_rate8'
                evname = {...
                    'freq400_freqdiff1100_rate8','freq3200_freqdiff1100_rate8',...
                    };
                
            case 'adapt_params_allconds'
                evname = read_conditions(exp,us,'adapt_params_fixed_stim');
            otherwise
                evname = {evgroup};
        end
        
    case 'naturalsound'
        switch evgroup
            case 'all_main'
                evname = read_conditions(exp,us,'main_combined',varargin{:});
            case 'all_main_v3'
                evname = read_conditions(exp,us,'main_v3_combined',varargin{:});
            case 'all_texture_combined'
                evname = read_conditions(exp,us,'texture_combined',varargin{:});
            case 'all_spectrotemporal_combined'
                evname = read_conditions(exp,us,'spectrotemporal_combined',varargin{:});
            case 'all_spectrotemporal_v2_combined'
                evname = read_conditions(exp,us,'spectrotemporal_v2_combined',varargin{:});
            case 'all_main_v3_runspec'
                evname = read_conditions(exp,us,'main_v3',varargin{:});
            case 'all_pitch_localizer'
                evname = read_conditions(exp,us,'pitch_localizer',varargin{:});
            case 'all_speech_localizer'
                evname = read_conditions(exp,us,'speech_localizer',varargin{:});
            case 'freq200-400'
                evname = {'freq200','freq800'};
                
            case 'freq200-6400'
                evname = {'freq200','freq6400'};
                
            case 'freq800-1600'
                evname = {'freq800','freq1600'};
                
            case 'freq3200-6400'
                evname = {'freq3200','freq6400'};
                
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq1600-3200-6400'
                evname = {'freq1600','freq3200','freq6400'};
                
            case 'freq800-1600-3200-6400'
                evname = {'freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-3200-6400'
                evname = {'freq200','freq400','freq3200','freq6400'};
                
            case 'freq200-400-800-1600'
                evname = {'freq200','freq400','freq800','freq1600'};
                
            case 'freq400-800-1600-3200-6400'
                evname = {'freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-800-1600-3200-6400'
                evname = {'freq200','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-1600-3200-6400'
                evname = {'freq200','freq400','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-800-3200-6400'
                evname = {'freq200','freq400','freq800','freq3200','freq6400'};
                
            case 'freq200-400-800-1600-6400'
                evname = {'freq200','freq400','freq800','freq1600','freq6400'};
                
            case 'freq200-400-800-1600-3200'
                evname = {'freq200','freq400','freq800','freq1600','freq3200'};
                
            case 'midi-bigband'
                evname = {'midi-intact','bigband-v2-intact'};
                
            case 'midi-bigband-quilt'
                evname = {'midi-quilt','bigband-v2-quilt'};
                
            case 'allfreq'
                evname = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'all_localizer'
                evname = read_conditions(exp,us,'localizer',varargin{:});
                
            case 'all_scrambling'
                evname = read_conditions(exp,us,'scrambling',varargin{:});
                
            case 'all_scrambling_russian'
                evname = read_conditions(exp,us,'scrambling_russian',varargin{:});
                
            case 'all_texture'
                evname = read_conditions(exp,us,'texture',varargin{:});
                
            case 'natural'
                evname = {'natural-1','natural-2'};
                
            case 'originals'
                evname = {'originals-1','originals-2'};
                
            case 'all_spectrotemporal'
                evname = read_conditions(exp,us,'spectrotemporal',varargin{:});
                
            case 'all_spectrotemporal_v2'
                evname = read_conditions(exp,us,'spectrotemporal_v2',varargin{:});

            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_adapt_params_v2'
        
        switch evgroup
            case 'freq3200'
                evname = {...
                    'freq3200_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel8',...
                    'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
                    'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
                    'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
                    };
            case 'freq400'
                evname = {...
                    'freq400_freqdiff0_sonelevel3','freq400_freqdiff0_sonelevel8',...
                    'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
                    'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
                    'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
                    };
            case 'sonelevel3'
                evname = {...
                    'freq3200_freqdiff0_sonelevel3',...
                    'freq3200_freqdiff275_sonelevel3',...
                    'freq3200_freqdiff550_sonelevel3',...
                    'freq3200_freqdiff1100_sonelevel3',...
                    'freq400_freqdiff0_sonelevel3',...
                    'freq400_freqdiff275_sonelevel3',...
                    'freq400_freqdiff550_sonelevel3',...
                    'freq400_freqdiff1100_sonelevel3',...
                    };
            case 'sonelevel8'
                evname = {...
                    'freq3200_freqdiff0_sonelevel8',...
                    'freq3200_freqdiff275_sonelevel8',...
                    'freq3200_freqdiff550_sonelevel8',...
                    'freq3200_freqdiff1100_sonelevel8',...
                    'freq400_freqdiff0_sonelevel8',...
                    'freq400_freqdiff275_sonelevel8',...
                    'freq400_freqdiff550_sonelevel8',...
                    'freq400_freqdiff1100_sonelevel8',...
                    };
            case 'freqdiff-0'
                evname = {...
                    'freq400_freqdiff0_sonelevel3','freq400_freqdiff0_sonelevel8',...
                    'freq3200_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel8',...
                    };
            case 'freqdiff-275-550-1100'
                evname = {...
                    'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
                    'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
                    'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
                    'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
                    'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
                    'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
                    };
            case 'freqdiff-1100'
                evname = {...
                    'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
                    'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
                    };
            case 'freqdiff-275'
                evname = {...
                    'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
                    'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
                    };
            case 'freqdiff-550'
                evname = {...
                    'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
                    'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
                    };
            case 'freq400_freqdiff0'
                evname = {...
                    'freq400_freqdiff0_sonelevel3','freq400_freqdiff0_sonelevel8',...
                    };
            case 'freq400_freqdiff275'
                evname = {...
                    'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
                    };
            case 'freq400_freqdiff550'
                evname = {...
                    'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
                    };
            case 'freq400_freqdiff1100'
                evname = {...
                    'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
                    };
            case 'freq3200_freqdiff0'
                evname = {...
                    'freq3200_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel8',...
                    };
            case 'freq3200_freqdiff275'
                evname = {...
                    'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
                    };
            case 'freq3200_freqdiff550'
                evname = {...
                    'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
                    };
            case 'freq3200_freqdiff1100'
                evname = {...
                    'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
                    };
                
            case 'freqdiff0_sonelevel3'
                evname = {...
                    'freq400_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel3',...
                    };
            case 'freqdiff275_sonelevel3'
                evname = {...
                    'freq400_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel3',...
                    };
            case 'freqdiff550_sonelevel3'
                evname = {...
                    'freq400_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel3',...
                    };
            case 'freqdiff1100_sonelevel3'
                evname = {...
                    'freq400_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel3',...
                    };
                
            case 'freqdiff0_sonelevel8'
                evname = {...
                    'freq400_freqdiff0_sonelevel8','freq3200_freqdiff0_sonelevel8',...
                    };
            case 'freqdiff275_sonelevel8'
                evname = {...
                    'freq400_freqdiff275_sonelevel8','freq3200_freqdiff275_sonelevel8',...
                    };
            case 'freqdiff550_sonelevel8'
                evname = {...
                    'freq400_freqdiff550_sonelevel8','freq3200_freqdiff550_sonelevel8',...
                    };
            case 'freqdiff1100_sonelevel8'
                evname = {...
                    'freq400_freqdiff1100_sonelevel8','freq3200_freqdiff1100_sonelevel8',...
                    };
                
            case 'adapt_params_allconds'
                evname = read_conditions(exp,us,'adapt_params_fixed_stim');
            otherwise
                evname = {evgroup};
        end
        
    case 'pitch_adapt_params_v3'
        
        switch evgroup
            %             case 'freq3200'
            %                 evname = {...
            %                     'freq3200_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel8',...
            %                     'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
            %                     'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
            %                     'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
            %                     };
            %             case 'freq400'
            %                 evname = {...
            %                     'freq400_freqdiff0_sonelevel3','freq400_freqdiff0_sonelevel8',...
            %                     'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
            %                     'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
            %                     'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
            %                     };
            %             case 'sonelevel3'
            %                 evname = {...
            %                     'freq3200_freqdiff0_sonelevel3',...
            %                     'freq3200_freqdiff275_sonelevel3',...
            %                     'freq3200_freqdiff550_sonelevel3',...
            %                     'freq3200_freqdiff1100_sonelevel3',...
            %                     'freq400_freqdiff0_sonelevel3',...
            %                     'freq400_freqdiff275_sonelevel3',...
            %                     'freq400_freqdiff550_sonelevel3',...
            %                     'freq400_freqdiff1100_sonelevel3',...
            %                     };
            %             case 'sonelevel8'
            %                 evname = {...
            %                     'freq3200_freqdiff0_sonelevel8',...
            %                     'freq3200_freqdiff275_sonelevel8',...
            %                     'freq3200_freqdiff550_sonelevel8',...
            %                     'freq3200_freqdiff1100_sonelevel8',...
            %                     'freq400_freqdiff0_sonelevel8',...
            %                     'freq400_freqdiff275_sonelevel8',...
            %                     'freq400_freqdiff550_sonelevel8',...
            %                     'freq400_freqdiff1100_sonelevel8',...
            %                     };
            %             case 'freqdiff-0'
            %                 evname = {...
            %                     'freq400_freqdiff0_sonelevel3','freq400_freqdiff0_sonelevel8',...
            %                     'freq3200_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel8',...
            %                     };
            %             case 'freqdiff-275-550-1100'
            %                 evname = {...
            %                     'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
            %                     'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
            %                     'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
            %                     'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
            %                     'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
            %                     'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
            %                     };
            %             case 'freqdiff-1100'
            %                 evname = {...
            %                     'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
            %                     'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
            %                     };
            %             case 'freqdiff-275'
            %                 evname = {...
            %                     'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
            %                     'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
            %                     };
            %             case 'freqdiff-550'
            %                 evname = {...
            %                     'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
            %                     'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
            %                     };
            %             case 'freq400_freqdiff0'
            %                 evname = {...
            %                     'freq400_freqdiff0_sonelevel3','freq400_freqdiff0_sonelevel8',...
            %                     };
            %             case 'freq400_freqdiff275'
            %                 evname = {...
            %                     'freq400_freqdiff275_sonelevel3','freq400_freqdiff275_sonelevel8',...
            %                     };
            %             case 'freq400_freqdiff550'
            %                 evname = {...
            %                     'freq400_freqdiff550_sonelevel3','freq400_freqdiff550_sonelevel8',...
            %                     };
            %             case 'freq400_freqdiff1100'
            %                 evname = {...
            %                     'freq400_freqdiff1100_sonelevel3','freq400_freqdiff1100_sonelevel8',...
            %                     };
            %             case 'freq3200_freqdiff0'
            %                 evname = {...
            %                     'freq3200_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel8',...
            %                     };
            %             case 'freq3200_freqdiff275'
            %                 evname = {...
            %                     'freq3200_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel8',...
            %                     };
            %             case 'freq3200_freqdiff550'
            %                 evname = {...
            %                     'freq3200_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel8',...
            %                     };
            %             case 'freq3200_freqdiff1100'
            %                 evname = {...
            %                     'freq3200_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel8',...
            %                     };
            %
            %             case 'freqdiff0_sonelevel3'
            %                 evname = {...
            %                     'freq400_freqdiff0_sonelevel3','freq3200_freqdiff0_sonelevel3',...
            %                     };
            %             case 'freqdiff275_sonelevel3'
            %                 evname = {...
            %                     'freq400_freqdiff275_sonelevel3','freq3200_freqdiff275_sonelevel3',...
            %                     };
            %             case 'freqdiff550_sonelevel3'
            %                 evname = {...
            %                     'freq400_freqdiff550_sonelevel3','freq3200_freqdiff550_sonelevel3',...
            %                     };
            %             case 'freqdiff1100_sonelevel3'
            %                 evname = {...
            %                     'freq400_freqdiff1100_sonelevel3','freq3200_freqdiff1100_sonelevel3',...
            %                     };
            %
            %             case 'freqdiff0_sonelevel8'
            %                 evname = {...
            %                     'freq400_freqdiff0_sonelevel8','freq3200_freqdiff0_sonelevel8',...
            %                     };
            %             case 'freqdiff275_sonelevel8'
            %                 evname = {...
            %                     'freq400_freqdiff275_sonelevel8','freq3200_freqdiff275_sonelevel8',...
            %                     };
            %             case 'freqdiff550_sonelevel8'
            %                 evname = {...
            %                     'freq400_freqdiff550_sonelevel8','freq3200_freqdiff550_sonelevel8',...
            %                     };
            %             case 'freqdiff1100_sonelevel8'
            %                 evname = {...
            %                     'freq400_freqdiff1100_sonelevel8','freq3200_freqdiff1100_sonelevel8',...
            %                     };
            %
            %             case 'adapt_params_allconds'
            %                 evname = read_conditions(exp,us,'adapt_params_fixed_stim');
            otherwise
                evname = {evgroup};
        end
        
    case {'pitch_f0adapt_v3'}
        
        switch evgroup
            case 'harm_f0diff-freqdiff'
                evname = {'freq200_test_harm_f0diff-freqdiff','freq800_test_harm_f0diff-freqdiff'};
            case 'harm_f0same-freqdiff'
                evname = {'freq200_test_harm_f0same-freqdiff','freq800_test_harm_f0same-freqdiff'};
            case 'harm_f0same-freqsame'
                evname = {'freq200_test_harm_f0same-freqsame','freq800_test_harm_f0same-freqsame'};
            case 'f0same-freqdiff'
                evname = {'freq200_test_harm_f0same-freqdiff','freq800_test_harm_f0same-freqdiff','freq200_test_inharm_f0same-freqdiff','freq800_test_inharm_f0same-freqdiff'};
            case 'f0same-freqsame'
                evname = {'freq200_test_harm_f0same-freqsame','freq800_test_harm_f0same-freqsame','freq200_test_inharm_f0same-freqsame','freq800_test_inharm_f0same-freqsame'};
            case 'all_f0adapt'
                evname = read_conditions(exp,us,'f0adapt');
            otherwise
                evname = {evgroup};
        end
        
    case {'pitch_f0adapt_v4'}
        
        switch evgroup
            case 'harm_f0diff-freqdiff'
                evname = {'freq297_test_harm_f0diff-freqdiff','freq440_test_harm_f0diff-freqdiff'};
            case 'harm_f0same-freqdiff'
                evname = {'freq297_test_harm_f0same-freqdiff','freq440_test_harm_f0same-freqdiff'};
            case 'harm_f0same-freqsame'
                evname = {'freq297_test_harm_f0same-freqsame','freq440_test_harm_f0same-freqsame'};
            case 'f0same-freqdiff'
                evname = {'freq297_test_harm_f0same-freqdiff','freq440_test_harm_f0same-freqdiff','freq297_test_inharm_f0same-freqdiff','freq440_test_inharm_f0same-freqdiff'};
            case 'f0same-freqsame'
                evname = {'freq297_test_harm_f0same-freqsame','freq440_test_harm_f0same-freqsame','freq297_test_inharm_f0same-freqsame','freq440_test_inharm_f0same-freqsame'};
            case 'all_f0adapt'
                evname = read_conditions(exp,us,'f0adapt');
            otherwise
                evname = {evgroup};
        end
        
    case {'pitch_f0adapt_v5'}
        
        switch evgroup
            case 'harm_f0diff-freqdiff'
                evname = {'freq297_test_harm_f0diff-freqdiff','freq200_test_harm_f0diff-freqdiff'};
            case 'harm_f0same-freqdiff'
                evname = {'freq297_test_harm_f0same-freqdiff','freq200_test_harm_f0same-freqdiff'};
            case 'harm_f0same-freqsame'
                evname = {'freq297_test_harm_f0same-freqsame','freq200_test_harm_f0same-freqsame'};
            case 'f0same-freqdiff'
                evname = {'freq297_test_harm_f0same-freqdiff','freq200_test_harm_f0same-freqdiff','freq297_test_inharm_f0same-freqdiff','freq200_test_inharm_f0same-freqdiff'};
            case 'f0same-freqsame'
                evname = {'freq297_test_harm_f0same-freqsame','freq200_test_harm_f0same-freqsame','freq297_test_inharm_f0same-freqsame','freq200_test_inharm_f0same-freqsame'};
            case 'all_f0adapt'
                evname = read_conditions(exp,us,'f0adapt');
            otherwise
                evname = {evgroup};
        end
        
    case {'pitch_f0adapt_v6'}
        
        switch evgroup
            case 'freq300_puretone_freqdiff-2p75-5p5-11'
                evname = {'freq300_puretone_freqdiff-2p75','freq300_puretone_freqdiff-5p5','freq300_puretone_freqdiff-11'};
            case 'freq3200_puretone_freqdiff-2p75-5p5-11'
                evname = {'freq3200_puretone_freqdiff-2p75','freq3200_puretone_freqdiff-5p5','freq3200_puretone_freqdiff-11'};
            case 'freq300'
                evname = {'freq300_puretone_freqdiff-0','freq300_puretone_freqdiff-2p75','freq300_puretone_freqdiff-5p5','freq300_puretone_freqdiff-11'};
            case 'freq3200'
                evname = {'freq3200_puretone_freqdiff-0','freq3200_puretone_freqdiff-2p75','freq3200_puretone_freqdiff-5p5','freq3200_puretone_freqdiff-11'};
            case 'all_freqf0adapt2'
                evname = read_conditions(exp,us,'freqf0adapt2');
            otherwise
                evname = {evgroup};
        end
        
        
    case {'pitch_f0adapt_v7'}
        
        switch evgroup
            case 'freq300_puretone_freqdiff-2p75-5p5-11'
                evname = {'freq300_puretone_freqdiff-2p75','freq300_puretone_freqdiff-5p5','freq300_puretone_freqdiff-11'};
            case 'all_freqf0adapt'
                evname = read_conditions(exp,us,'freqf0adapt');
            case 'freq200-400'
                evname = {'freq200','freq800'};
                
            case 'freq200-6400'
                evname = {'freq200','freq6400'};
                
            case 'freq800-1600'
                evname = {'freq800','freq1600'};
                
            case 'freq3200-6400'
                evname = {'freq3200','freq6400'};
                
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq1600-3200-6400'
                evname = {'freq1600','freq3200','freq6400'};
                
            case 'freq800-1600-3200-6400'
                evname = {'freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-3200-6400'
                evname = {'freq200','freq400','freq3200','freq6400'};
                
            case 'freq200-400-800-1600'
                evname = {'freq200','freq400','freq800','freq1600'};
                
            case 'freq400-800-1600-3200-6400'
                evname = {'freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-800-1600-3200-6400'
                evname = {'freq200','freq800','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-1600-3200-6400'
                evname = {'freq200','freq400','freq1600','freq3200','freq6400'};
                
            case 'freq200-400-800-3200-6400'
                evname = {'freq200','freq400','freq800','freq3200','freq6400'};
                
            case 'freq200-400-800-1600-6400'
                evname = {'freq200','freq400','freq800','freq1600','freq6400'};
                
            case 'freq200-400-800-1600-3200'
                evname = {'freq200','freq400','freq800','freq1600','freq3200'};
                
            case 'allfreq'
                evname = {'freq200','freq400','freq800','freq1600','freq3200','freq6400'};
                
            case 'all_localizer'
                evname = read_conditions(exp,us,'localizer',varargin{:});
                
            otherwise
                evname = {evgroup};
        end
        
    case {'tonotopy_monkey'}
        
        switch evgroup
            
            case 'freq200-400-800'
                evname = {'freq200','freq400','freq800'};
                
            case 'freq3200-6400-12800'
                evname = {'freq3200','freq6400','freq12800'};
                
            case 'allfreq'
                evname = {'freq200','freq400','freq800','freq1600','freq3200','freq6400','freq12800'};
                
            otherwise
                evname = {evgroup};
        end
        
    case {'pitch_localizer_monkey'}
        
        switch evgroup
            

            case 'harm-f0100-200'
                evname = {'harm-f0100', 'harm-f0200'};
                
            case 'noise-f0100-200'
                evname = {'noise-f0100', 'noise-f0200'};
                
            case 'harm-f0800-1600'
                evname = {'harm-f0800', 'harm-f01600'};
                
            case 'noise-f0800-1600'
                evname = {'noise-f0800', 'noise-f01600'};
            
            case 'harm'
                evname = {'harm-f0100', 'harm-f0200', 'harm-f0400', 'harm-f0800', 'harm-f01600'};
                
            case 'noise'
                evname = {'noise-f0100','noise-f0200','noise-f0400','noise-f0800','noise-f01600'};
                
            case 'f0100-200'
                evname = {'harm-f0100', 'noise-f0100','harm-f0200', 'noise-f0200'};
                
            case 'f0800-1600'
                evname = {'harm-f0800', 'noise-f0800','harm-f01600', 'noise-f01600'};
                
            case 'all_pitchloc2'
                evname = read_conditions(exp,us,'pitchloc2',varargin{:});
                
            case 'harm-res'
                evname = {'harm-res-freq1200-2400','harm-res-freq2400-4800','harm-res-freq4800-9600'};
                
            case 'harm-res-ploc4'
                evname = {'harm-res-freq1200-2400','harm-res-freq2400-4800'};
                
            case 'harm-unres'
                evname = {'harm-unres-freq1200-2400','harm-unres-freq2400-4800','harm-unres-freq4800-9600'};
                
            case 'harm-unres-ploc4'
                evname = {'harm-unres-freq1200-2400','harm-unres-freq2400-4800'};
                
            case 'noise-pitchloc3'
                evname = {'noise-freq1200-2400','noise-freq2400-4800','noise-freq4800-9600'};
                
            case 'noise-ploc4'
                evname = {'noise-freq1200-2400','noise-freq2400-4800'};
                
            case 'noise-mod-ploc4'
                evname = {'noise-mod-freq1200-2400','noise-mod-freq2400-4800'};
                
            case '1200-2400'
                evname = {'harm-res-freq1200-2400', 'harm-unres-freq1200-2400','noise-freq1200-2400'};
                
            case '4800-9600'
                evname = {'harm-res-freq4800-9600', 'harm-unres-freq14800-9600','noise-freq4800-9600'};
                
            case '1200-2400-ploc4'
                evname = {'harm-res-freq1200-2400', 'harm-unres-freq1200-2400','noise-freq1200-2400','noise-mod-freq1200-2400'};
                
            case '2400-4800-ploc4'
                evname = {'harm-res-freq2400-4800', 'harm-unres-freq2400-4800','noise-freq2400-4800','noise-mod-freq2400-4800'};
                
            case 'all_pitchloc3'
                evname = read_conditions(exp,us,'pitchloc3',varargin{:});
                
            case 'all_pitchloc4'
                evname = read_conditions(exp,us,'pitchloc4',varargin{:});
                
            case 'all_pitchloc6'
                evname = read_conditions(exp,us,'pitchloc6',varargin{:});
                
                
            
            case 'harm-f0100-70-80dB'
                evname = {'harm-f0100-70dB','harm-f0100-75dB','harm-f0100-80dB'};
            
            case 'harm-f0200-70-80dB'
                evname = {'harm-f0200-70dB','harm-f0200-75dB','harm-f0200-80dB'};
                
            case 'harm-f0400-70-80dB'
                evname = {'harm-f0400-70dB','harm-f0400-75dB','harm-f0400-80dB'};
            
            case 'harm-f0800-70-80dB'
                evname = {'harm-f0800-70dB','harm-f0800-75dB','harm-f0800-80dB'};
                
            case 'harm-f01600-70-80dB'
                evname = {'harm-f01600-70dB','harm-f01600-75dB','harm-f01600-80dB'};
                
                
                
            case 'noise-f0100-70-80dB'
                evname = {'noise-f0100-70dB','noise-f0100-75dB','noise-f0100-80dB'};
                
            case 'noise-f0200-70-80dB'
                evname = {'noise-f0200-70dB','noise-f0200-75dB','noise-f0200-80dB'};
                
            case 'noise-f0400-70-80dB'
                evname = {'noise-f0400-70dB','noise-f0400-75dB','noise-f0400-80dB'};
                
            case 'noise-f0800-70-80dB'
                evname = {'noise-f0800-70dB','noise-f0800-75dB','noise-f0800-80dB'};
                
            case 'noise-f01600-70-80dB'
                evname = {'noise-f01600-70dB','noise-f01600-75dB','noise-f01600-80dB'};
                
                
                
                
            case 'noise-f0200-65-75dB'
                evname = {'noise-f0200-65dB','noise-f0200-70dB','noise-f0200-75dB'};
                
            case 'harm-f0200-75-80dB'
                evname = {'harm-f0200-75dB','harm-f0200-80dB'};
                
            case 'noise-f0200-65-70dB'
                evname = {'noise-f0200-65dB','noise-f0200-70dB'};
                
            case 'harm-f0200-70-75dB'
                evname = {'harm-f0200-70dB','harm-f0200-75dB'};
                
            case 'noise-f0200-70-75dB'
                evname = {'noise-f0200-70dB','noise-f0200-75dB'};                
                
            case 'harm-noise-f0200-70-75dB'
                evname = {'harm-f0200-70dB','harm-f0200-75dB','noise-f0200-70dB','noise-f0200-75dB'};
                
            case 'all_pitchloc7_70-75-80dB'
                evname = read_conditions(exp,us,'pitchloc7_70-75-80dB',varargin{:});
                
            case 'f0100-200_70-80dB'
                evname = {...
                    'harm-f0100-70dB', 'noise-f0100-70dB',...
                    'harm-f0100-75dB', 'noise-f0100-75dB',...
                    'harm-f0100-80dB', 'noise-f0100-80dB',...
                    'harm-f0200-70dB', 'noise-f0200-70dB',...
                    'harm-f0200-75dB', 'noise-f0200-75dB',...
                    'harm-f0200-80dB', 'noise-f0200-80dB'};    

            case 'f0800-1600_70-80dB'
                evname = {...
                    'harm-f0800-70dB',  'noise-f0800-70dB',...
                    'harm-f0800-75dB',  'noise-f0800-75dB',...
                    'harm-f0800-80dB',  'noise-f0800-80dB',...
                    'harm-f01600-70dB', 'noise-f01600-70dB',...
                    'harm-f01600-75dB', 'noise-f01600-75dB',...
                    'harm-f01600-80dB', 'noise-f01600-80dB'};  
                
            case 'harm_70-80dB'
                evname = {...
                    'harm-f0100-70dB', 'harm-f0200-70dB',...
                    'harm-f0100-75dB', 'harm-f0200-75dB',...
                    'harm-f0100-80dB', 'harm-f0200-80dB',...
                    'harm-f0400-70dB', 'harm-f0800-70dB',...
                    'harm-f0400-75dB', 'harm-f0800-75dB',...
                    'harm-f0400-80dB', 'harm-f0800-80dB',...
                    'harm-f01600-70dB', ...
                    'harm-f01600-75dB', ...
                    'harm-f01600-80dB'};
                
            case 'noise_70-80dB'
                evname = {...
                    'noise-f0100-70dB', 'noise-f0200-70dB',...
                    'noise-f0100-75dB', 'noise-f0200-75dB',...
                    'noise-f0100-80dB', 'noise-f0200-80dB',...
                    'noise-f0400-70dB', 'noise-f0800-70dB',...
                    'noise-f0400-75dB', 'noise-f0800-75dB',...
                    'noise-f0400-80dB', 'noise-f0800-80dB',...
                    'noise-f01600-70dB', ...
                    'noise-f01600-75dB', ...
                    'noise-f01600-80dB'};
                
            case 'harm_70-75dB'
                evname = {...
                    'harm-f0100-70dB', 'harm-f0200-70dB',...
                    'harm-f0100-75dB', 'harm-f0200-75dB',...
                    'harm-f0400-70dB', 'harm-f0800-70dB',...
                    'harm-f0400-75dB', 'harm-f0800-75dB',...
                    'harm-f01600-70dB', ...
                    'harm-f01600-75dB'};
                
            case 'noise_75-80dB'
                evname = {...
                    'noise-f0100-75dB', 'noise-f0200-75dB',...
                    'noise-f0100-80dB', 'noise-f0200-80dB',...
                    'noise-f0400-75dB', 'noise-f0800-75dB',...
                    'noise-f0400-80dB', 'noise-f0800-80dB',...
                    'noise-f01600-75dB', ...
                    'noise-f01600-80dB'};
                
            case '70dB'
                evname = {...
                    'harm-f0100-70dB',  'noise-f0100-70dB',...
                    'harm-f0200-70dB',  'noise-f0200-70dB',...
                    'harm-f0400-70dB',  'noise-f0400-70dB',...
                    'harm-f0800-70dB',  'noise-f0800-70dB',...
                    'harm-f01600-70dB', 'noise-f01600-70dB'};
                
            case '80dB'
                evname = {...
                    'harm-f0100-80dB',  'noise-f0100-80dB',...
                    'harm-f0200-80dB',  'noise-f0200-80dB',...
                    'harm-f0400-80dB',  'noise-f0400-80dB',...
                    'harm-f0800-80dB',  'noise-f0800-80dB',...
                    'harm-f01600-80dB', 'noise-f01600-80dB'};
                
            otherwise
                evname = {evgroup};
        end
        
    case {'voice_localizer_human'}
        
        switch evgroup
            
            case 'all_mvocs_pitch'
                evname = read_conditions(exp,us,'mvocs_pitch',varargin{:});
                
            case 'all-unpitched'
                evname = {'pitched-noisevoc','unpitched-matched','unpitched-matched-noisevoc'};
                
            case 'all_mvocs_pitch_v2'
                evname = read_conditions(exp,us,'mvocs_pitch_v2',varargin{:});
                
            case 'all_mvocs_pitch_v3'
                evname = read_conditions(exp,us,'mvocs_pitch_v3',varargin{:});
            
            case 'all_mvocs_pitch_v4'
                evname = read_conditions(exp,us,'mvocs_pitch_v4',varargin{:});
                
            case 'all_mvocs_pitch_v5'
                evname = read_conditions(exp,us,'mvocs_pitch_v5',varargin{:});
                
            case 'all_mvocs_pitch_v4_in_monkey'
                evname = read_conditions('voice_localizer_monkey',[],'mvocs_pitch_v4',varargin{:});
                
            case 'all-pitched-v2'
                evname = {'pitched','pitched-straight-harmvoc'};
                
            case 'all-unpitched-v2'
                evname = {'unpitched-matched','pitched-straight-noisevoc','pitched-straight-noisevoc-5dB'};
                
            case 'all-noisevoc-v2'
                evname = {'pitched-straight-noisevoc','pitched-straight-noisevoc-5dB'};
                
            case 'pitched-straight-harmvoc-noisevoc-5dB'
                evname = {'pitched-straight-harmvoc','pitched-straight-noisevoc-5dB'};
                
            case 'pitched-straight-harmvoc-65-80dB'
                evname = {'pitched-straight-harmvoc-80dB', 'pitched-straight-harmvoc-75dB', 'pitched-straight-harmvoc-70dB', 'pitched-straight-harmvoc-65dB'};
            
            case 'pitched-straight-noisevoc-65-80dB'
                evname = {'pitched-straight-noisevoc-80dB', 'pitched-straight-noisevoc-75dB', 'pitched-straight-noisevoc-70dB', 'pitched-straight-noisevoc-65dB'};
                
            case 'pitched-straight-harmvoc-noisevoc-65-80dB'
                evname = {...
                    'pitched-straight-harmvoc-80dB', 'pitched-straight-harmvoc-75dB', 'pitched-straight-harmvoc-70dB','pitched-straight-harmvoc-65dB',...
                    'pitched-straight-noisevoc-80dB', 'pitched-straight-noisevoc-75dB', 'pitched-straight-noisevoc-70dB', 'pitched-straight-noisevoc-65dB'};
                
            case 'pitched-straight-harmvoc-70-75dB'
                evname = {'pitched-straight-harmvoc-75dB', 'pitched-straight-harmvoc-70dB'};
                
            case 'pitched-straight-noisevoc-70-75dB'
                evname = {'pitched-straight-noisevoc-75dB', 'pitched-straight-noisevoc-70dB'};
                
            case 'pitched-straight-harmvoc-noisevoc-70-75dB'
                evname = {'pitched-straight-harmvoc-75dB', 'pitched-straight-harmvoc-70dB', 'pitched-straight-noisevoc-75dB', 'pitched-straight-noisevoc-70dB'};
            
                
            otherwise
                evname = {evgroup};
        end
        
    case {'voice_localizer_monkey'}
        
        switch evgroup
            
            case 'all_mvocs_pitch'
                evname = read_conditions(exp,us,'mvocs_pitch',varargin{:});
                
            case 'all-unpitched'
                evname = {'pitched-noisevoc','unpitched-matched','unpitched-matched-noisevoc'};
                
            case 'all_mvocs_pitch_v2'
                evname = read_conditions(exp,us,'mvocs_pitch_v2',varargin{:});
                
            case 'all_mvocs_pitch_v3'
                evname = read_conditions(exp,us,'mvocs_pitch_v3',varargin{:});
                                
            case 'all_mvocs_pitch_v4'
                evname = read_conditions(exp,us,'mvocs_pitch_v4',varargin{:});

            case 'all_mvocs_pitch_v5'
                evname = read_conditions(exp,us,'mvocs_pitch_v4',varargin{:});
                
            case 'all-pitched-v2'
                evname = {'pitched','pitched-straight-harmvoc'};
                
            case 'all-unpitched-v2'
                evname = {'unpitched-matched','pitched-straight-noisevoc','pitched-straight-noisevoc-5dB'};
                
            case 'all-noisevoc-v2'
                evname = {'pitched-straight-noisevoc','pitched-straight-noisevoc-5dB'};
                
            case 'pitched-straight-harmvoc-noisevoc-5dB'
                evname = {'pitched-straight-harmvoc','pitched-straight-noisevoc-5dB'};
                
            case 'pitched-straight-harmvoc-70-80dB'
                evname = {'pitched-straight-harmvoc-70dB','pitched-straight-harmvoc-75dB','pitched-straight-harmvoc-80dB'};
                
            case 'pitched-straight-noisevoc-65-75dB'
                evname = {'pitched-straight-noisevoc-65dB','pitched-straight-noisevoc-70dB','pitched-straight-noisevoc-75dB'};
                
            case 'pitched-straight-harmvoc-75-80dB'
                evname = {'pitched-straight-harmvoc-75dB','pitched-straight-harmvoc-80dB'};
                
            case 'pitched-straight-noisevoc-65-70dB'
                evname = {'pitched-straight-noisevoc-65dB','pitched-straight-noisevoc-70dB'};
                
            case 'pitched-straight-harmvoc-70-75dB'
                evname = {'pitched-straight-harmvoc-70dB','pitched-straight-harmvoc-75dB'};
                
            case 'pitched-straight-noisevoc-70-75dB'
                evname = {'pitched-straight-noisevoc-70dB','pitched-straight-noisevoc-75dB'};      
                
            case 'pitched-straight-harmvoc-noisevoc-70-75dB'
                evname = {'pitched-straight-harmvoc-75dB', 'pitched-straight-harmvoc-70dB',...
                    'pitched-straight-noisevoc-75dB', 'pitched-straight-noisevoc-70dB'};
                
            case 'pitched-straight-harmvoc-noisevoc-65-80dB'
                evname = {...
                    'pitched-straight-harmvoc-80dB', 'pitched-straight-harmvoc-75dB', 'pitched-straight-harmvoc-70dB','pitched-straight-harmvoc-65dB',...
                    'pitched-straight-noisevoc-80dB', 'pitched-straight-noisevoc-75dB', 'pitched-straight-noisevoc-70dB', 'pitched-straight-noisevoc-65dB'};
            
            case 'pitched-straight-harmvoc-65-80dB'
                evname = {'pitched-straight-harmvoc-80dB', 'pitched-straight-harmvoc-75dB', 'pitched-straight-harmvoc-70dB', 'pitched-straight-harmvoc-65dB'};
            
            case 'pitched-straight-noisevoc-65-80dB'
                evname = {'pitched-straight-noisevoc-80dB', 'pitched-straight-noisevoc-75dB', 'pitched-straight-noisevoc-70dB', 'pitched-straight-noisevoc-65dB'};
                
            otherwise
                evname = {evgroup};
        end
        
    case {'resting_monkey'}
        
        switch evgroup
            
            case 'harm'
                evname = {'harm-f0100', 'harm-f0200', 'harm-f0400', 'harm-f0800', 'harm-f01600'};
                
            case 'noise'
                evname = {'noise-f0100','noise-f0200','noise-f0400','noise-f0800','noise-f01600'};
                
            case 'f0100-200'
                evname = {'harm-f0100', 'noise-f0100','harm-f0200', 'noise-f0200'};
                
            case 'f0800-1600'
                evname = {'harm-f0800', 'noise-f0800','harm-f01600', 'noise-f01600'};
                
            case 'all_rest'
                evname = read_conditions(exp,us,'rest',varargin{:});
                
            otherwise
                evname = {evgroup};
        end
        
        
        
    case {'pitch_localizer_human'}
        
        switch evgroup
            
            case 'harm-f0100-200'
                evname = {'harm-f0100', 'harm-f0200'};
                
            case 'noise-diotic-f0100-200'
                evname = {'noise-diotic-f0100', 'noise-diotic-f0200'};
                
            case 'harm-f0800-1600'
                evname = {'harm-f0800', 'harm-f01600'};
                
            case 'noise-diotic-f0800-1600'
                evname = {'noise-diotic-f0800', 'noise-diotic-f01600'};
            
            case 'harm'
                evname = {'harm-f0100', 'harm-f0200', 'harm-f0400', 'harm-f0800', 'harm-f01600'};
                
            case 'noise'
                evname = {...
                    'noise-dichotic-f0100','noise-dichotic-f0200','noise-dichotic-f0400','noise-dichotic-f0800','noise-dichotic-f01600',...
                    'noise-diotic-f0100','noise-diotic-f0200','noise-diotic-f0400','noise-diotic-f0800','noise-diotic-f01600',...
                    };
                
            case 'noise-dichotic'
                evname = {...
                    'noise-dichotic-f0100','noise-dichotic-f0200','noise-dichotic-f0400','noise-dichotic-f0800','noise-dichotic-f01600',...
                    };
                
            case 'noise-diotic'
                evname = {...
                    'noise-diotic-f0100','noise-diotic-f0200','noise-diotic-f0400','noise-diotic-f0800','noise-diotic-f01600',...
                    };
                                
            case 'f0100-200'
                evname = {'harm-f0100', 'noise-dichotic-f0100', 'noise-diotic-f0100','harm-f0200', 'noise-dichotic-f0200', 'noise-diotic-f0200'};
                
            case 'f0100-200-harm-noise-dichotic'
                evname = {'harm-f0100', 'noise-dichotic-f0100', 'harm-f0200', 'noise-dichotic-f0200'};
                
            case 'f0100-200-harm-noise-diotic'
                evname = {'harm-f0100', 'noise-diotic-f0100','harm-f0200', 'noise-diotic-f0200'};
                
            case 'f0100-harm-noise-diotic'
                evname = {'harm-f0100', 'noise-diotic-f0100'};
                
            case 'harm-noise-diotic'
                evname = {'harm-f0100', 'noise-diotic-f0100', 'harm-f0200', 'noise-diotic-f0200', 'harm-f0400', 'noise-diotic-f0400', 'harm-f0800', 'noise-diotic-f0800', 'harm-f01600', 'noise-diotic-f01600'};
             
            case 'harm-noise-dichotic'
                evname = {'harm-f0100', 'noise-dichotic-f0100', 'harm-f0200', 'noise-dichotic-f0200', 'harm-f0400', 'noise-dichotic-f0400', 'harm-f0800', 'noise-dichotic-f0800', 'harm-f01600', 'noise-dichotic-f01600'};
              
            case 'f0800-1600'
                evname = {'harm-f0800', 'noise-dichotic-f0800', 'noise-diotic-f0800', 'harm-f01600', 'noise-dichotic-f01600', 'noise-diotic-f01600'};
                
            case 'f0800-1600-harm-noise-dichotic'
                evname = {'harm-f0800', 'noise-dichotic-f0800', 'harm-f01600', 'noise-dichotic-f01600'};
                
            case 'f0800-1600-harm-noise-diotic'
                evname = {'harm-f0800', 'noise-diotic-f0800','harm-f01600', 'noise-diotic-f01600'};
                
            case 'all_pitchloc2'
                evname = read_conditions(exp,us,'pitchloc2',varargin{:});
                
            case 'all_pitchloc6'
                evname = read_conditions(exp,us,'pitchloc6',varargin{:});
                
            case 'all_pitchloc6_in_monkey'
                evname = read_conditions('pitch_localizer_monkey',us,'pitchloc6',varargin{:});
                
            case 'harm-f0200-65-80dB'
                evname = {'harm-f0200-80dB', 'harm-f0200-75dB', 'harm-f0200-70dB', 'harm-f0200-65dB'};
            
            case 'noise-f0200-65-80dB'
                evname = {'noise-f0200-80dB', 'noise-f0200-75dB', 'noise-f0200-70dB', 'noise-f0200-65dB'};
                            
            case 'harm-f0200-70-75dB'
                evname = {'harm-f0200-75dB', 'harm-f0200-70dB'};
            
            case 'noise-f0200-70-75dB'
                evname = {'noise-f0200-75dB', 'noise-f0200-70dB'};
                
            case 'harm-noise-f0200-70-75dB'
                evname = {'harm-f0200-70dB','harm-f0200-75dB','noise-f0200-70dB','noise-f0200-75dB'};
                
            otherwise
                evname = {evgroup};
        end
        
        
    case {'naturalsound_monkey'}
        
        switch evgroup
            
            case 'animals-german'
                evname = {'animals','german'};
                
            case 'animals-macaques'
                evname = {'animals','macaques'};
                
            case 'all_voice_localizer'
                evname = read_conditions(exp,us,'voice_localizer',varargin{:});
                
            case 'all_main'
                evname = read_conditions(exp,us,'main',varargin{:});
                
            case 'all_main_combined'
                evname = read_conditions(exp,us,'main_combined',varargin{:});
                
            case 'MVocsControl'
                evname = {'AVocs','EnvSounds','MVocsPhaseScram'};
                
            case 'all_petkov_localizer'
                evname = read_conditions(exp,us,'petkov_localizer',varargin{:});
                
            otherwise
                evname = {evgroup};
        end
        
    case {'color_monkey'}
        evname = {evgroup};
    otherwise
        error('No valid experiment');
        
end

assert(length(evname) == length(unique(evname)));
% v1