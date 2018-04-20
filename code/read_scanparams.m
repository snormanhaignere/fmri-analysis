function varargout = read_scanparams(exp,us,runtype,varargin)

hrf_type = 'BOLD';
disdaqs = 0;
nullwin = [];
switch exp
    
    case 'pitch_event'
        
        switch runtype
            case 'main'
                if us == 32
                    nTR = 97;
                    TR = 4;
                    
                elseif us == 26 && optInputs(varargin, 'run') && varargin{optInputs(varargin, 'run')+1} == 4
                    nTR = 176;
                    TR = 2.0;
                    
                elseif us == 26
                    nTR = '187';
                    TR = '2.0';
                else
                    nTR = '195';
                    TR = '2.0';
                end
                
            case 'localizer'
                
                if s == 32
                    nTR = 174;
                    TR = 2;
                    
                elseif us == 26 && optInputs(varargin, 'run') && varargin{optInputs(varargin, 'run')+1} == 3
                    nTR = 187;
                    TR = 2;
                    
                elseif s == 26
                    nTR = 176;
                    TR = 2;
                    
                else
                    nTR = 216;
                    TR = 2;
                    
                end
        end
        
        % design parameters
        blockdur = 17;
        nulldur = 17;
        TA = 1;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur) + stim2scan;
        
    case 'pitch_resthr'
        
        % design parameters
        blockdur = 17;
        nulldur = 17;
        TR = 3.4;
        TA = 1;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur) + stim2scan;
        nTR = 142;
        
    case 'pitch_f0'
        
        switch runtype
            case 'adapt6note'
                % design parameters
                blockdur = 20.4;
                nulldur = 20.4;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur-TR) + stim2scan;
                nTR = 97;
            case 'localizer'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur-TR) + stim2scan;
                nTR = 111;
        end
        
        
    case 'pitch_f0adapt'
        
        switch runtype
            case {'f0adapt','localizer'}
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 79;
            otherwise
                error('No valid runtype')
        end
        
    case 'pitch_f0inharm'
        
        % design parameters
        blockdur = 17;
        nulldur = 17;
        TR = 3.4;
        TA = 1;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur+TR*2) + stim2scan;
        nTR = 154;
        
    case 'pitch_f0smallfreq'
        
        % design parameters
        blockdur = 13.6;
        nulldur = 13.6;
        TR = 3.4;
        TA = 1;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur+TR*2) + stim2scan;
        nTR = 124;
        
    case 'pitch_f02by2'
        
        % design parameters
        blockdur = 17;
        nulldur = 17;
        TR = 3.4;
        TA = 1;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur+TR*2) + stim2scan;
        nTR = 104;
        
    case 'pitch_natural'
        
        if any(us == [19 13 23 24 25 26 27 28])
            
            % design parameters
            blockdur = 25.6;
            nulldur = 19.2;
            TR = 3.2;
            TA = 0.8;
            stimdur = 2;
            stim2scan = 2.2;
            win = (-TR:TR:blockdur) + stim2scan;
            
        else
            
            error('Bad subject ID');
            
        end
        
    case 'pitch_contour'
        
        switch runtype
            case 'main'
                % design parameters
                blockdur = 17.6;
                nulldur = 17.6;
                TR = 4.4;
                TA = 2;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur-TR) + stim2scan;
            case 'localizer'
                % design parameters
                blockdur = 17.6;
                nulldur = 17.6;
                TR = 4.4;
                TA = 2;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur-TR) + stim2scan;
        end
        
    case 'pitch_dp'
        
        % design parameters
        blockdur = 25.6;
        nulldur = 19.2;
        TR = 3.2;
        TA = 0.8;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur+TR) + stim2scan;
        
    case 'mspec'
        
        % design parameters
        blockdur = 9;
        nulldur = NaN;
        TR = 2;
        TA = 2;
        stimdur = 2;
        stim2scan = 0;
        win = (0:TR:blockdur) + stim2scan;
        
    case 'tonotopy_rate'
        
        % design parameters
        blockdur = 17;
        nulldur = 17;
        TR = 3.4;
        TA = 1;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur+TR*2) + stim2scan;
        nTR = 129;
        
    case 'sepi'
        
        switch runtype
            
            case 'sparse'
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 94;
                
            case 'continuous_sepi'
                
                % design parameters
                blockdur = 18;
                nulldur = 18;
                TR = 3;
                TA = 3;
                stimdur = 3;
                stim2scan = 0;
                win = (0:TR:blockdur+TR*1) + stim2scan;
                nTR = 111;
                
            case 'continuous'
                
                % design parameters
                blockdur = 18;
                nulldur = 18;
                TR = 3;
                TA = 3;
                stimdur = 3;
                stim2scan = 0;
                win = (0:TR:blockdur+TR*1) + stim2scan;
                nTR = 111;
        end
        
    case {'pitch_resthr_v2','pitch_resthr_v4'}
        
        switch runtype
            
            case 'main'
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 114;
                
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    r = varargin{optInputs(varargin, 'run')+1};
                    if us == 68 && any(r == 7)
                        nTR = 113;
                    end
                end
                
            case 'localizer'
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 129;
        end
        
    case 'pitch_resthr_v3'
        
        switch runtype
            
            case 'main'
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 134;
                
        end
        
    case 'pitch_misc'
        
        switch runtype
            
            case 'schroeder'
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 104;
                
        end
        
    case 'amusia'
        
        if any(us == [43 44])
            % design parameters
            blockdur = 20.4;
            nulldur = 20.4;
            TR = 3.4;
            TA = 1;
            stimdur = 2;
            stim2scan = 2.2;
            win = (TR:TR:blockdur-TR) + stim2scan;
            nTR = 121;
        else
            blockdur = 20.4;
            nulldur = 20.4;
            TR = 3.4;
            TA = 1;
            stimdur = 2;
            stim2scan = 2.2;
            win = (-TR:TR:blockdur+TR) + stim2scan;
            nTR = 155;
        end
        
    case 'music_12channel'
        blockdur = 16;
        nulldur = 16;
        TR = 2;
        TA = 2;
        stimdur = 0;
        stim2scan = 0;
        win = (-TR:TR:blockdur+2*TR) + stim2scan;
        nTR = 201;
        
    case 'stax'
        % design parameters
        blockdur = 20.4655;
        nulldur = 20.4655;
        TR = 4.0931;
        TA = 1.6931;
        stimdur = 2;
        stim2scan = 2.2;
        win = (-TR:TR:blockdur-TR) + stim2scan;
        nTR = 79;
        disdaqs = 4;
        
    case 'pitch_dp_v2'
        % design parameters
        blockdur = 16.3724;
        nulldur = 16.3724;
        TR = 4.0931;
        TA = 1.6931;
        stimdur = 2;
        stim2scan = 2.4;
        win = (-TR:TR:blockdur+2*TR) + stim2scan;
        disdaqs = 4;
        nTR = 92;
        if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
            r = varargin{optInputs(varargin, 'run')+1};
            if us == 3 && any(r == [1 2])
                nTR = 88;
            end
        end
        
    case 'multiband'
        
        switch runtype
            
            case {'singleband_partial','multiband_wholebrain'}
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 79;
                
            case 'singleband_wholebrain'
                
                % design parameters
                blockdur = 17.6;
                nulldur = 17.6;
                TR = 4.4;
                TA = 2;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 64;
        end
        
    case 'music_scram'
        
        % design parameters
        blockdur = 10;
        nulldur = 10;
        TR = 3.5-0.0069;
        TA = 3.5-0.0069;
        stimdur = 10;
        stim2scan = 0;
        win = (-TR:TR:TR*5) + stim2scan;
        nTR = 104;
        disdaqs = 4;
        
    case 'music_scram_familiar'
        
        % design parameters
        blockdur = 10;
        nulldur = 10;
        TR = 3.5;
        TA = 3.5;
        stimdur = 10;
        stim2scan = 0;
        win = (-TR:TR:TR*5) + stim2scan;
        nTR = 132;
        disdaqs = 4;
        
    case 'pitch_overlap'
        
        switch runtype
            
            case {'overlap','overlap_v2'}
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 114;
                
            case {'tonotopy'}
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 79;
                
            otherwise
                error('No valid experiment');
                
        end
        
    case 'pitch_overlap_v2'
        
        switch runtype
            
            case {'overlap_v2_lowconds','overlap_v2_highconds'}
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 79;
                
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    r = varargin{optInputs(varargin, 'run')+1};
                    if us == 86 && strcmp(runtype, 'overlap_v2_highconds') && any(r == 4)
                        nTR = 78;
                    end
                end
                
            case {'overlap_v2'}
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*1) + stim2scan;
                nTR = 79*2;
                
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    r = varargin{optInputs(varargin, 'run')+1};
                    if us == 86 && any(r == 4)
                        nTR = 79*2-1;
                    end
                end
                
            otherwise
                error('No valid experiment');
                
        end
        
    case 'pitch_overlap_v3'
        
        switch runtype
            
            case {'overlap_v3_lowconds','overlap_v3_highconds'}
                
                % design parameters
                blockdur = 19.5;
                nulldur = 19.5;
                TR = 3.9;
                TA = 1.5;
                stimdur = 2;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 64;
                disdaqs = 4;
                
            case {'overlap_v3'}
                
                % design parameters
                blockdur = 19.5;
                nulldur = 19.5;
                TR = 3.9;
                TA = 1.5;
                stimdur = 2;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 64*2;
                
            case {'localizer'}
                
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 79;
                
            otherwise
                error('No valid runtype');
                
        end
        
    case 'pitch_f0adapt_v2'
        
        switch runtype
            
            case 'f0adapt'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 114;
                
            otherwise
                error('No valid runtype');
        end
        
    case 'pitch_adapt_params'
        switch runtype
            case 'adapt_params_fixed_stim'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 104;
                
            case 'adapt_params_variable_stim'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 104;
        end
        
    case 'naturalsound'
        switch runtype
            case 'main'
                % design parameters
                blockdur = 22;
                nulldur = 22;
                TR = 4.4;
                TA = 2;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 79;
                
            case 'main_v2'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 104;
                
            case 'main_v3'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 99;
                
            case 'main_combined'
                % design parameters
                blockdur = 22;
                nulldur = 22;
                TR = 4.4;
                TA = 2;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 79*5;
                
            case 'main_v2_combined'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 104*4;
                
            case {'main_v3_combined', 'main_v3_combined_raw'}
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 99*11;
                
            case 'natsoundloc'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 79;
            case 'localizer'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 79;
            case 'pitch_localizer'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 79;
            case 'speech_localizer'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 69;
            case 'scrambling'
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 104;
                if us == 164 && optInputs(varargin, 'run') && varargin{optInputs(varargin, 'run')+1} == 1
                    nTR = 100;
                    disdaqs = 4;
                elseif us == 164 && optInputs(varargin, 'run') && any(varargin{optInputs(varargin, 'run')+1} == [2 3])
                    disdaqs = 4;
                end
            case 'scrambling_russian'
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 104;
            case 'texture'
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 94;
            case {'texture_combined','texture_combined_raw'}
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 94*12;
            case 'spectrotemporal'
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 109;
            case 'spectrotemporal_v2'
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 94;
            case {'spectrotemporal_v2_combined','spectrotemporal_v2_combined_raw'}
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 94*12;
            case {'spectrotemporal_combined','spectrotemporal_combined_raw'}
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 109*10;
            case {'quilting_v1_continuous'}
                blockdur = 12.24;
                nulldur = 12.24;
                TR = 2.04;
                TA = 2.04;
                stimdur = 2.04; % NOT CLEAR WHAT THIS DOES HERE
                stim2scan = 0;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                disdaqs = 4;
                nTR = 245;
            otherwise
                error('No valid runtype');
        end
    case 'pitch_adapt_params_v2'
        switch runtype
            case 'adapt_params_variable_stim'
                % design parameters
                blockdur = 13.6;
                nulldur = 13.6;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                nTR = 100;
        end
    case 'pitch_adapt_params_v3'
        switch runtype
            case 'adapt_params'
                % design parameters
                blockdur = [];
                nulldur = 27;
                TR = 1.5;
                TA = 1.5;
                stimdur = [];
                stim2scan = 0;
                win = (0:TR:18);
                nTR = 226;
                disdaqs = 4;
        end
        
    case {'pitch_f0adapt_v3','pitch_f0adapt_v4','pitch_f0adapt_v5'}
        switch runtype
            case 'f0adapt'
                % design parameters
                blockdur = [];
                nulldur = 18;
                TR = 18;
                TA = 1;
                stimdur = [];
                stim2scan = 4.5;
                win = 4.5 + [-18,0];
                nullwin = 16.75;
                nTR = 17;
        end
    case {'pitch_f0adapt_v6'}
        switch runtype
            case 'freqf0adapt2'
                % design parameters
                blockdur = [];
                nulldur = 26;
                TR = 26;
                TA = 1.2;
                stimdur = [];
                stim2scan = 26-(1.2+0.25);
                win = stim2scan + [-TR,0];
                nullwin = stim2scan;
                nTR = 11;
        end
        
    case {'pitch_f0adapt_v7'}
        switch runtype
            case 'freqf0adapt'
                % design parameters
                blockdur = [];
                nulldur = 26;
                TR = 26;
                TA = 1.2;
                stimdur = [];
                stim2scan = 26-(1.2+0.3);
                win = stim2scan + [-TR,0];
                nullwin = stim2scan;
                nTR = 11;
                
            case 'localizer'
                % design parameters
                blockdur = 17;
                nulldur = 17;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 79;
        end
    case {'tonotopy_monkey'}
        hrf_type = 'MION';
        switch runtype
            case 'localizer_1000msTA'
                % design parameters
                blockdur = 23.8;
                nulldur = 23.8;
                TR = 3.4;
                TA = 1.2;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 130;
        end
    case {'pitch_localizer_monkey'}
        hrf_type = 'MION';
        switch runtype
            case {'pitchloc2'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 141;
            case {'pitchloc2_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 141;
                hrf_type = 'MION_CUSTOM1';
            case 'pitchloc3'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 131;
            case 'pitchloc4'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 121;
            case 'pitchloc2_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'pitchloc2', 'scans', scan);
                    nTR = 141*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
            case {'pitchloc6','pitchloc6_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 161;
                hrf_type = 'MION_CUSTOM1';
                
            case 'pitchloc6_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'pitchloc6', 'scans', scan);
                    nTR = 161*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc6_combined_raw');
                end
                
            case {'pitchloc7_70dB','pitchloc7_75dB','pitchloc7_80dB'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 141;
                hrf_type = 'MION_CUSTOM1';
                
            case {'pitchloc7_70-75-80dB','pitchloc7_70-75-80dB_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 141*3;
                hrf_type = 'MION_CUSTOM1';
                
            case 'pitchloc7_70-75-80dB_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run')
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(...
                        exp, us, 'pitchloc7_70-75-80dB', 'scans', scan);
                    nTR = 141*length(runs);
                else
                    error(['Need to supply run in read_scanparams.m' ...
                        'for pitchloc7_combined_raw']);
                end
        end
        
    case {'voice_localizer_monkey'}
        hrf_type = 'MION';
        switch runtype
            case {'mvocs_pitch','mvocs_pitch_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 161;
                hrf_type = 'MION_CUSTOM1';
            case {'mvocs_pitch_v2','mvocs_pitch_v2_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 191;
                hrf_type = 'MION_CUSTOM1';
            case {'mvocs_pitch_v3','mvocs_pitch_v3_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 161;
                hrf_type = 'MION_CUSTOM1';
            case {'mvocs_pitch_v4','mvocs_pitch_v4_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 161;
                hrf_type = 'MION_CUSTOM1';
                
            case {'mvocs_pitch_v5','mvocs_pitch_v5_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 111;
                hrf_type = 'MION_CUSTOM1';
                
            case 'mvocs_pitch_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'mvocs_pitch', 'scans', scan);
                    nTR = 161*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
            case 'mvocs_pitch_v2_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run')
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'mvocs_pitch_v2', 'scans', scan);
                    nTR = 191*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
            case 'mvocs_pitch_v3_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run')
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'mvocs_pitch_v3', 'scans', scan);
                    nTR = 161*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for mvocs_pitch_v3_combined_raw');
                end
            case 'mvocs_pitch_v4_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run')
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'mvocs_pitch_v4', 'scans', scan);
                    nTR = 161*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for mvocs_pitch_v3_combined_raw');
                end
                
            case 'mvocs_pitch_v5_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'MION_CUSTOM1';
                if optInputs(varargin, 'run')
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'mvocs_pitch_v5', 'scans', scan);
                    nTR = 111*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for mvocs_pitch_v3_combined_raw');
                end
        end
        
        
    case {'voice_localizer_human'}
        hrf_type = 'BOLD';
        switch runtype
            case {'mvocs_pitch_v2','mvocs_pitch_v2_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 191;
            case 'mvocs_pitch_v2_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'mvocs_pitch', 'scans', scan);
                    nTR = 191*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
            case {'mvocs_pitch_v4','mvocs_pitch_v4_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 111;
            case 'mvocs_pitch_v4_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                if optInputs(varargin, 'run')
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'mvocs_pitch', 'scans', scan);
                    nTR = 111*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
        end
        
    case {'resting_monkey'}
        hrf_type = 'MION_CUSTOM';
        switch runtype
            case {'rest','rest_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 2;
                TA = 2;
                stimdur = 0;
                stim2scan = 0;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 240;
            case 'rest_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 2;
                TA = 2;
                stimdur = 0;
                stim2scan = 0;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'rest', 'scans', scan);
                    nTR = 141*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
        end
        
    case {'pitch_localizer_human'}
        hrf_type = 'BOLD';
        switch runtype
            case {'pitchloc2','pitchloc2_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 191;
            case {'pitchloc2_combined_raw'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'pitchloc2', 'scans', scan);
                    nTR = 191*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
                
            case {'pitchloc6','pitchloc6_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 111;
                hrf_type = 'BOLD';
                
            case 'pitchloc6_combined_raw'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2.4;
                stim2scan = 2.4;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                hrf_type = 'BOLD';
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'pitchloc6', 'scans', scan);
                    nTR = 111*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc6_combined_raw');
                end
                
            otherwise
                error('Runtype does not match those in script.')
        end
        
    case {'naturalsound_monkey'}
        hrf_type = 'MION';
        switch runtype
            case 'main'
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 201;
            case {'main_combined','main_combined_raw'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 201*13;
            case {'voice_localizer','petkov_localizer','petkov_localizer_combined_split'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                nTR = 161;
            case {'petkov_localizer_combined_raw'}
                % design parameters
                blockdur = 34;
                nulldur = 34;
                TR = 3.4;
                TA = 1;
                stimdur = 2;
                stim2scan = 2.2;
                win = (-TR:TR:blockdur+TR) + stim2scan;
                if optInputs(varargin, 'run') %% subject 3 only collected 88 TRs
                    scan = varargin{optInputs(varargin, 'run')+1};
                    runs = read_runs(exp, us, 'petkov_localizer', 'scans', scan);
                    nTR = 161*length(runs);
                else
                    error('Need to supply run in read_scanparams.m for pitchloc2_combined_raw');
                end
        end
    case {'color_monkey'}
        
        hrf_type = 'MION';
        blockdur = 34;
        nulldur = 34;
        TR = 2;
        TA = 0;
        stimdur = 2;
        stim2scan = 0;
        win = (-TR:TR:blockdur+TR) + stim2scan;
        nTR = 272;
                
    case {'naturalsound-infant'}
        
        % design parameters
        blockdur = 12;
        nulldur = 12;
        TR = 3;
        TA = 3;
        stimdur = 2.4;
        stim2scan = 0;
        win = (-TR:TR:blockdur+TR*2) + stim2scan;
        disdaqs = 3;
        
        if isempty(optInputs(varargin, 'run'))
            error('Need to supply run in read_scanparams.m for naturalsound-infant');
        end
            
        run = varargin{optInputs(varargin, 'run')+1};
        switch run
            case 1
                nTR = 379 - disdaqs;
            case 2
                nTR = 615 - disdaqs;
            otherwise
                error('No matching run.');
        end
        
    case {'tono-pitch-localizer'}
        
        % design parameters
        blockdur = 2.96*3;
        nulldur = 2.96*3;
        TR = 2.96;
        TA = 0.56;
        stimdur = 2;
        stim2scan = 0.2;
        win = (-TR:TR:blockdur+TR*2) + stim2scan;
        disdaqs = 4;
        nTR = 98-disdaqs;
        
    case {'tono-localizer'}
        
        % design parameters
        switch us
            case {367}
                blockdur = 2.96*3;
                nulldur = 2.96*3;
                TR = 2.96;
                TA = 0.56;
                stimdur = 2;
                stim2scan = 0.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                disdaqs = 4;
                nTR = 98-disdaqs;
            case {375}
                blockdur = 2.97*3;
                nulldur = 2.97*3;
                TR = 2.97;
                TA = 0.57;
                stimdur = 2;
                stim2scan = 0.2;
                win = (-TR:TR:blockdur+TR*2) + stim2scan;
                disdaqs = 4;
                nTR = 98-disdaqs;
            otherwise
                error('Scan parameters are subject specific due to a slight switch in parameter values for 3T1')
        end
        
    case {'naturalsound-nmf-localizer'}
        
        % design parameters
        blockdur = 2.96*3;
        nulldur = 2.96*3;
        TR = 2.96;
        TA = 0.56;
        stimdur = 2;
        stim2scan = 0.2;
        win = (-TR:TR:blockdur+TR*2) + stim2scan;
        disdaqs = 3;
        nTR = 104-disdaqs;
        
    otherwise
        error('No valid experiment');
end

varargout = {blockdur, nulldur, TR, TA, stimdur, stim2scan, win, nTR, disdaqs, nullwin, hrf_type};

if optInputs(varargin, 'hrf_type')
    varargout = {hrf_type};
end

if optInputs(varargin, 'TR')
    varargout = {TR};
end

if optInputs(varargin, 'nTR')
    varargout = {nTR};
end

if optInputs(varargin, 'TA')
    varargout = {TA};
end