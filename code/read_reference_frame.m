function reference_frame = read_reference_frame(exp, us, runtype, r, varargin)

[~, ~, ~, ~, ~, ~, ~, nTR] = read_scanparams(exp,us,runtype,'run',r,varargin{:});

switch exp
    case 'tonotopy_monkey'
        if us == 157 && strcmp(runtype, 'localizer_1000msTA') && any(r == [5 9 14 17])
            reference_frame = 2;
        elseif us == 158 && strcmp(runtype, 'localizer_1000msTA') && any(r == [2 5 11])
            reference_frame = 6;
        elseif us == 158 && strcmp(runtype, 'localizer_1000msTA') && any(r == [8])
            reference_frame = 19;
        elseif us == 158 && strcmp(runtype, 'localizer_1000msTA') && any(r == [16])
            reference_frame = 14;
        else
            reference_frame = 0;
        end
        fprintf('Using reference frame %d for %s run %d\n', reference_frame, runtype, r);
        
    case 'pitch_localizer_monkey'
        if us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [14])
            reference_frame = 21;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [22])
            reference_frame = 11;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [23 27 30])
            reference_frame = 2;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [31 43])
            reference_frame = 6;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [32])
            reference_frame = 5;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [34])
            reference_frame = 15;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [35])
            reference_frame = 9;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [46])
            reference_frame = 7;
            
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [123])
            reference_frame = 2;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [136])
            reference_frame = 32;
        elseif us == 157 && strcmp(runtype, 'pitchloc2') && any(r == [137])
            reference_frame = 7;

            
            
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [2])
            reference_frame = 14;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [12])
            reference_frame = 21;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [13])
            reference_frame = 23;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [27])
            reference_frame = 17;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [30])
            reference_frame = 1;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [32])
            reference_frame = 9;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [35])
            reference_frame = 16;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [41 49])
            reference_frame = 10;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [45])
            reference_frame = 4;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [48])
            reference_frame = 4;
        elseif us == 158 && strcmp(runtype, 'pitchloc2') && any(r == [50])
            reference_frame = 12;
            
        elseif us == 157 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 1)
            reference_frame = 1;
        elseif us == 157 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 2)
            reference_frame = 2;
        elseif us == 157 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 3)
            reference_frame = 1;
        elseif us == 157 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 4)
            reference_frame = 2275;
        elseif us == 157 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 5)
            reference_frame = 34;
            
        elseif us == 158 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 1)
            reference_frame = 1;
        elseif us == 158 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 2)
            reference_frame = 0;
        elseif us == 158 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 3)
            reference_frame = 2;
        elseif us == 158 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 4)
            reference_frame = 3;
        elseif us == 158 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 5)
            reference_frame = 1;
        elseif us == 158 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 6)
            reference_frame = 0;
            
        elseif us == 170 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 1)
            reference_frame = 0;
        elseif us == 170 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 2)
            reference_frame = 1;
        elseif us == 170 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 3)
            reference_frame = 0;
        elseif us == 170 && strcmp(runtype, 'pitchloc2_combined_raw') && any(r == 4)
            reference_frame = 0;
            
            
        elseif us == 158 && strcmp(runtype, 'pitchloc6_combined_raw') && any(r == 9)
            reference_frame = 0+161; % first image of second run
        elseif us == 158 && strcmp(runtype, 'pitchloc6_combined_raw') && any(r == 10)
            reference_frame = 0+5*161; % first image of fifth run (run "47")
          
        elseif us == 372 && strcmp(runtype, 'pitchloc7_70-75-80dB_combined_raw') && any(r == 4)
            reference_frame = 5; % 6th image of run 13 (scan 4)
        elseif us == 372 && strcmp(runtype, 'pitchloc7_70-75-80dB_combined_raw') && any(r == 5)
            reference_frame = 5; % 6th image of run 18 (scan 5)
        elseif us == 372 && strcmp(runtype, 'pitchloc7_70-75-80dB_combined_raw') && any(r == 6)
            reference_frame = 9; % 6th image of run 18 (scan 5)
            
        else
            reference_frame = 0;
        end
        
        
    case 'voice_localizer_monkey'
        if us == 157 && strcmp(runtype, 'mvocs_pitch_v3_combined_raw') && any(r == 2)
            reference_frame = 41;
        elseif us == 158 && strcmp(runtype, 'mvocs_pitch_v3_combined_raw') && any(r == 4)
            reference_frame = 161+12; % 13th image of second run
        elseif us == 157 && strcmp(runtype, 'mvocs_pitch_v4_combined_raw') && any(r == 4)
            reference_frame = 90; % 13th image of second run
        elseif us == 170 && strcmp(runtype, 'mvocs_pitch_v4_combined_raw') && any(r == 3)
            reference_frame = 5; % 13th image of second run
        elseif us == 157 && strcmp(runtype, 'mvocs_pitch_v4_combined_raw') && any(r == 5)
            reference_frame = 31; % 13th image of second run
        elseif us == 157 && strcmp(runtype, 'mvocs_pitch_v4_combined_raw') && any(r == 6)
            reference_frame = 147+161; % 148th image of second run
        elseif us == 170 && strcmp(runtype, 'mvocs_pitch_v4_combined_raw') && any(r == 4)
            reference_frame = 14+161; % 14th image of second run
        else
            reference_frame = 0;
        end
        
    case 'naturalsound_monkey'
        if us == 157 && strcmp(runtype, 'main') && any(r == [8])
            reference_frame = 1;
        elseif us == 158 && strcmp(runtype, 'main_combined_raw') && any(r == [1 3])
            reference_frame = 0;            
        elseif us == 158 && strcmp(runtype, 'main_combined_raw') && r == 2
            reference_frame = 3;
        elseif us == 157 && strcmp(runtype, 'main_combined_raw') && r == 1
            reference_frame = 6;            
        elseif us == 157 && strcmp(runtype, 'main_combined_raw') && r == 2
            reference_frame = 9;  
        elseif us == 157 && strcmp(runtype, 'main_combined_raw') && r == 3
            reference_frame = 0;  
        else
            reference_frame = 0;
        end
        
    case 'naturalsound'
        if strcmp(runtype, 'main_v3_combined_raw')
            reference_frame = 0;
        else
            reference_frame = floor(nTR/2);
        end
        
    case 'color_monkey'
        reference_frame = 0;
        
    otherwise
        reference_frame = floor(nTR/2);
end