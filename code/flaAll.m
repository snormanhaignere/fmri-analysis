    function flaAll(exp,usubs,model,stage,varargin)

% stfAll(exp,usubs,'noplot',varargin{:});
for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp,usubs(i),'func');
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(runtypes)
        
        runnum = read_runs(exp,usubs(i),runtypes{j},varargin{:});
        if optInputs(varargin, 'runnum');
            runnum = varargin{optInputs(varargin,'runnum')+1};
        end
        
        for k = 1:length(runnum)
            
            fprintf('%s, subject %d, %s, run %d, %s\n',exp,usubs(i),runtypes{j},runnum(k),model);
            
            switch exp
                case {'tono-pitch-localizer', 'tono-localizer', 'voice_localizer_human','voice_localizer_monkey','amusia','color_monkey','naturalsound_monkey','pitch_localizer_monkey','pitch_localizer_human','tonotopy_monkey','pitch_f0adapt_v3','pitch_f0adapt_v4','pitch_f0adapt_v5','pitch_f0adapt_v6','pitch_f0adapt_v7','pitch_dp_v2','music_scram','music_scram_familiar','pitch_overlap_v2','pitch_overlap_v3','pitch_adapt_params','naturalsound','pitch_adapt_params_v2','pitch_adapt_params_v3'}
                    fla_stats(exp,usubs(i),runtypes{j},runnum(k),model,stage,varargin{:});
                case {}
                    fla(exp,usubs(i),runtypes{j},runnum(k),model,stage,varargin{:});
                otherwise
                    error('No valid experiment');
            end
        end
    end
end
