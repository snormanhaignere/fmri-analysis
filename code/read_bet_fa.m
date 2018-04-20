function favalue = read_bet_fa(exp, us, runtype, r, varargin)

switch exp
    %     case 'pitch_overlap_v3'
    %         switch runtype
    %             case 'localizer'
    %                 favalue = 0.15;
    %             case 'overlap_v3'
    %                 favalue = 0.15;
    %             otherwise
    %                 error('No valid runtype');
    %         end
    case 'tonotopy_monkey'
        if us == 157
            favalue = 0.3;
        elseif us == 158 && strcmp(runtype, 'localizer_1000msTA') && r == 1
            favalue = 0.35;
        elseif us == 158 && strcmp(runtype, 'localizer_1000msTA') && r == 2
            favalue = 0.5;
        elseif us == 158 && strcmp(runtype, 'localizer_1000msTA') && r == 3
            favalue = 0.25;
        elseif us == 158 && strcmp(runtype, 'localizer_1000msTA') && any(r == [5 7 8 9 10 12 13 15 16 17 18 19])
            favalue = 0.3; % check 
        else
            favalue = 0.25;
        end
        
    case 'pitch_localizer_monkey'
        if us == 157 && strcmp(runtype, 'pitchloc2') && r == 1
            favalue = 0.3;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [131 134 137])
            favalue = 0.25;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [123 125 129 132])
            favalue = 0.35;
        elseif us == 158 && strcmp(runtype, 'pitchloc6_combined_raw') && r == 9
            favalue = 0.4;
        elseif us == 158 && strcmp(runtype, 'pitchloc6_combined_raw') && r == 10
            favalue = 0.4;
        elseif us == 372 && strcmp(runtype, 'pitchloc7_70-75-80dB_combined_raw') && r == 2
            favalue = 0.4;
        elseif us == 373 && strcmp(runtype, 'pitchloc7_70-75-80dB_combined_raw') && r == 3
            favalue = 0.4;
        else
            favalue = 0.3;
        end
        
    case 'voice_localizer_monkey'
        if us == 158 && strcmp(runtype, 'mvocs_pitch_combined_raw') && any(r == 1)
            favalue = 0.4;
        elseif us == 158 && strcmp(runtype, 'mvocs_pitch_v3_combined_raw') && any(r == 4)
            favalue = 0.4;
        elseif us == 158 && strcmp(runtype, 'mvocs_pitch_v4_combined_raw') && any(r == 6)
            favalue = 0.4;
        elseif us == 157 && strcmp(runtype, 'mvocs_pitch_v4_combined_raw') && any(r == 6)
            favalue = 0.4;
        else
            favalue = 0.3;
        end
        
    case 'naturalsound_monkey'
        if us == 158 && strcmp(runtype, 'voice_localizer') && any(r == [1 4 6])
            favalue = 0.35;
        elseif us == 157 && strcmp(runtype, 'main_combined_raw') && r == 3
            favalue = 0.5;
        elseif us == 158 && strcmp(runtype, 'petkov_localizer_combined_raw')
            favalue = 0.5;
        else
            favalue = 0.3;
        end
        
    case 'resting_monkey'
        favalue = 0.3;
        
    case 'color_monkey'
        favalue = 0.3;
        
    otherwise
        favalue = 0.15;
end