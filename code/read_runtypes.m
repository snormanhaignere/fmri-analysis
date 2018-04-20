function runtypes = read_runtypes(exp,us,varargin)

% runtypes = read_runtypes(s,varargin)
% returns all of the runtypes for a given subject/experiment

% keyboard;
switch exp
    
    case {'pitch_event'}
        if optInputs(varargin,'func');
            runtypes = {'main','localizer'};
        else
            runtypes = {'main','localizer','struct'};
        end
        
    case {'pitch_f0inharm','pitch_f0smallfreq','pitch_f02by2'}
        
        if optInputs(varargin,'func');
            runtypes = {'f0adapt'};
        else
            runtypes = {'f0adapt','struct'};
        end
        
    case 'pitch_f0'
        
        if optInputs(varargin,'func');
            runtypes = {'adapt6note','localizer'};
        else
            runtypes = {'adapt6note','localizer','struct'};
        end
        
    case {'pitch_dp','pitch_synthetic','pitch_natural','pitch_speech','pitch_resthr','mspec','tonotopy_rate','pitch_resthr_v3','music_12channel'}
        
        if optInputs(varargin,'func');
            runtypes = {'main'};
        else
            runtypes = {'main','struct'};
        end
        
    case 'pitch_contour'
        
        if optInputs(varargin,'func');
            runtypes = {'main','localizer'};
        else
            runtypes = {'main','localizer','struct'};
        end
        
    case 'sepi'
        
        if optInputs(varargin,'func');
            runtypes = {'continuous_sepi','continuous','sparse'};
        else
            runtypes = {'continuous_sepi','continuous','sparse','struct'};
        end
        
    case {'pitch_resthr_v2','pitch_resthr_v4'}
        
        if optInputs(varargin,'func');
            runtypes = {'main','localizer'};
        else
            runtypes = {'main','localizer','struct'};
        end
        
    case 'amusia'
        
        if optInputs(varargin,'func');
            runtypes = {'localizer'};
        else
            runtypes = {'localizer','struct'};
        end
        
    case 'pitch_misc'
        
        if optInputs(varargin,'func');
            runtypes = {'schroeder'};
        else
            runtypes = {'schroeder','struct'};
        end
        
    case 'stax'
        
        if optInputs(varargin,'func');
            runtypes = {'32ch_quiet_stax'};
        else
            runtypes = {'32ch_quiet_stax'};
        end
        
    case 'multiband'
        
        if optInputs(varargin,'func');
            runtypes = {'singleband_partial','singleband_wholebrain','multiband_wholebrain'};
        else
            runtypes = {'singleband_partial','singleband_wholebrain','multiband_wholebrain','struct'};
        end
        
    case 'pitch_dp_v2'
        runtypes = {'main'};
        
    case 'music_scram'
        
        if optInputs(varargin,'func');
            runtypes = {'main'};
        else
            runtypes = {'main','struct'};
        end
        
    case 'pitch_f0adapt'
        
        if optInputs(varargin,'func');
            runtypes = {'f0adapt','localizer'};
        else
            runtypes = {'f0adapt','localizer','struct'};
        end
        
    case 'music_scram_familiar'
        
        if optInputs(varargin,'func');
            runtypes = {'main'};
        else
            runtypes = {'main','struct'};
        end
        
    case 'pitch_overlap'
        
        if any(us == [3 84])
            if optInputs(varargin,'func');
                runtypes = {'overlap'};
            else
                runtypes = {'overlap','struct'};
            end
        else
            if optInputs(varargin,'func');
                runtypes = {'overlap_v2','tonotopy'};
            else
                runtypes = {'overlap_v2','tonotopy','struct'};
            end
        end
        
    case 'pitch_overlap_v2'
        
        if optInputs(varargin,'raw')
            runtypes = {'overlap_v2_lowconds','overlap_v2_highconds','struct'};
        elseif optInputs(varargin,'func');
            runtypes = {'overlap_v2'};
        else
            runtypes = {'overlap_v2','struct'};
        end
        
    case 'pitch_overlap_v3'
        
        if optInputs(varargin,'raw')
            runtypes = {'overlap_v3_lowconds','overlap_v3_highconds','localizer','struct'};
        elseif optInputs(varargin,'func');
            runtypes = {'overlap_v3','overlap_v3_lowconds','overlap_v3_highconds','localizer'};
        else
            runtypes = {'overlap_v3','overlap_v3_lowconds','overlap_v3_highconds','localizer','struct'};
        end
        
    case 'pitch_f0adapt_v2'
        
        if optInputs(varargin,'func');
            runtypes = {'f0adapt'};
        else
            runtypes = {'f0adapt','struct'};
        end
        
    case {'pitch_f0adapt_v3','pitch_f0adapt_v4','pitch_f0adapt_v5'}
        
        if optInputs(varargin,'func');
            runtypes = {'f0adapt'};
        else
            runtypes = {'f0adapt','struct'};
        end
        
    case {'pitch_f0adapt_v6'}
        
        if optInputs(varargin,'func');
            runtypes = {'freqf0adapt2'};
        else
            runtypes = {'freqf0adapt2','struct'};
        end
        
    case {'pitch_f0adapt_v7'}
        
        if optInputs(varargin,'func');
            runtypes = {'freqf0adapt','localizer'};
        else
            runtypes = {'freqf0adapt','localizer','struct'};
        end
        
    case 'pitch_adapt_params'
        
        if optInputs(varargin,'func');
            runtypes = {'adapt_params_fixed_stim','adapt_params_variable_stim'};
        else
            runtypes = {'adapt_params_fixed_stim','adapt_params_variable_stim','struct'};
        end
        
    case 'naturalsound'
        if optInputs(varargin,'func');
            runtypes = {'quilting_v1_continuous','main_combined','main_v2_combined','main_v3_combined','natsoundloc','localizer','scrambling','scrambling_russian','texture','texture_combined','spectrotemporal','spectrotemporal_combined','spectrotemporal_v2','spectrotemporal_v2_combined'};
        elseif optInputs(varargin, 'raw_func')
            runtypes = {'quilting_v1_continuous','main_v3','main','main_v2','natsoundloc','localizer','scrambling','scrambling_russian','texture','spectrotemporal','spectrotemporal_v2'};
        else
            runtypes = {'quilting_v1_continuous','main_v3','main','main_v2','natsoundloc','localizer','scrambling','scrambling_russian','texture','spectrotemporal','spectrotemporal_v2','struct'};
        end
        
        if us ~= 1
            runtypes(ismember(runtypes,{'natsoundloc','main','main_v2','main_combined','main_v2_combined'})) = [];
        end
        
        if us == 136
            if optInputs(varargin,'func');
                runtypes = {'pitch_localizer','speech_localizer'};
            elseif optInputs(varargin, 'raw_func')
                runtypes = {'pitch_localizer','speech_localizer'};
            else
                runtypes = {'pitch_localizer','speech_localizer','struct'};
            end
        end
        
    case 'pitch_adapt_params_v2'
        
        if optInputs(varargin,'func');
            runtypes = {'adapt_params_variable_stim'};
        else
            runtypes = {'adapt_params_variable_stim','struct'};
        end
        
    case 'pitch_adapt_params_v3'
        
        if optInputs(varargin,'func');
            runtypes = {'adapt_params'};
        else
            runtypes = {'adapt_params','struct'};
        end
        
    case 'tonotopy_monkey'
        if optInputs(varargin,'func');
            runtypes = {'localizer_1000msTA'};
        else
            runtypes = {'localizer_1000msTA'};
        end
        
    case 'pitch_localizer_monkey'
        switch us
            case {157, 158, 170}
                if optInputs(varargin,'func');
                    runtypes = {'pitchloc2','pitchloc3','pitchloc4',...
                        'pitchloc2_combined_raw','pitchloc2_combined_split',...
                        'pitchloc6_combined_raw','pitchloc6_combined_split'};
                elseif optInputs(varargin,'raw_func');
                    runtypes = {'pitchloc2','pitchloc3','pitchloc4','pitchloc6'};
                else
                    runtypes = {'pitchloc2','pitchloc3','pitchloc4','pitchloc6'};
                end
            case {372, 373}
                if optInputs(varargin,'func');
                    runtypes = {...
                        'pitchloc7_70dB','pitchloc7_75dB','pitchloc7_80dB'...
                        'pitchloc7_70-75-80dB',...
                        'pitchloc7_70-75-80dB_combined_raw',...
                        'pitchloc7_70-75-80dB_combined_split',...
                        };
                elseif optInputs(varargin,'raw_func');
                    runtypes = {...
                        'pitchloc7_70dB','pitchloc7_75dB','pitchloc7_80dB'};
                else
                    runtypes = {...
                        'pitchloc7_70dB','pitchloc7_75dB','pitchloc7_80dB'};
                end
            otherwise
                error('No matching case for us %d\n', us);
        end
        
    case 'voice_localizer_monkey'
        if optInputs(varargin,'func');
            runtypes = {'mvocs_pitch','mvocs_pitch_v2','mvocs_pitch_v3','mvocs_pitch_v4','mvocs_pitch_v5','mvocs_pitch_combined_raw','mvocs_pitch_v2_combined_raw','mvocs_pitch_v3_combined_raw','mvocs_pitch_v4_combined_raw','mvocs_pitch_v5_combined_raw','mvocs_pitch_combined_split','mvocs_pitch_v2_combined_split','mvocs_pitch_v3_combined_split','mvocs_pitch_v4_combined_split','mvocs_pitch_v5_combined_split'};
        elseif optInputs(varargin,'raw_func');
            runtypes = {'mvocs_pitch','mvocs_pitch_v2','mvocs_pitch_v3','mvocs_pitch_v4','mvocs_pitch_v5'};
        else
            runtypes = {'mvocs_pitch','mvocs_pitch_v2','mvocs_pitch_v3','mvocs_pitch_v4','mvocs_pitch_v5'};
        end
        
    case 'voice_localizer_human'
        if optInputs(varargin,'func');
            runtypes = {'mvocs_pitch_v2','mvocs_pitch_v2_combined_raw','mvocs_pitch_v2_combined_split','mvocs_pitch_v4','mvocs_pitch_v4_combined_raw','mvocs_pitch_v4_combined_split'};
        elseif optInputs(varargin,'raw_func');
            runtypes = {'mvocs_pitch_v4'};
        else
            runtypes = {'mvocs_pitch_v4','struct'};
        end
        
    case 'resting_monkey'
        if optInputs(varargin,'func');
            runtypes = {'rest','rest_combined_raw','rest_combined_split'};
        elseif optInputs(varargin,'raw_func');
            runtypes = {'rest'};
        else
            runtypes = {'rest'};
        end
                
    case 'pitch_localizer_human'
        if optInputs(varargin,'func');
            runtypes = {'pitchloc2','pitchloc2_combined_raw','pitchloc2_combined_split','pitchloc6_combined_raw','pitchloc6_combined_split','pitchloc6'};
        elseif optInputs(varargin,'raw_func');
            runtypes = {'pitchloc2','pitchloc6'};
        else
            runtypes = {'pitchloc2','pitchloc6','struct'};
        end
        
    case 'tono-pitch-localizer'
        if optInputs(varargin,'func');
            runtypes = {'tono_pitch_localizer','tono_pitch_localizer_combined_raw','tono_pitch_localizer_combined_split'};
        elseif optInputs(varargin,'raw_func');
            runtypes = {'tono_pitch_localizer'};
        else
            runtypes = {'tono_pitch_localizer','struct'};
        end
        
    case 'tono-localizer'
        if optInputs(varargin,'func');
            runtypes = {'tono_localizer','tono_localizer_combined_raw',...
                'tono_localizer_combined_split'};
        elseif optInputs(varargin,'raw_func');
            runtypes = {'tono_localizer'};
        else
            runtypes = {'tono_localizer','struct'};
        end
        
    case 'naturalsound_monkey'
        %         if optInputs(varargin,'func');
        %             runtypes = {'main'};
        %         else
        %             runtypes = {'main'};
        %         end
        if us == 158;
            runtypes = {'main','voice_localizer','petkov_localizer'};
        else
            runtypes = {'main','petkov_localizer'};
        end
        
    case 'color_monkey'
        %         if optInputs(varargin,'func');
        %             runtypes = {'main'};
        %         else
        %             runtypes = {'main'};
        %         end
        runtypes = {'main'};
        
        
    case 'naturalsound-infant'
        
        runtypes = {'main_v1'};% TESTING
        
    case 'naturalsound-nmf-localizer'
        
        if optInputs(varargin,'func');
            runtypes = {'main_v1'};
        elseif optInputs(varargin,'raw_func');
            runtypes = {'main_v1'};
        else
            runtypes = {'main_v1','struct'};
        end
        
end