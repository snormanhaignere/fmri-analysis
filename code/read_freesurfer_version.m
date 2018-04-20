function version = read_freesurfer_version(exp,varargin)

switch exp
    case {'tono-localizer', 'tono-pitch-localizer','voice_localizer_human','resting_monkey','color_monkey','naturalsound_monkey','pitch_localizer_monkey','voice_localizer_monkey','pitch_localizer_human','tonotopy_monkey','pitch_f0adapt_v3','pitch_f0adapt_v4','pitch_f0adapt_v5','pitch_f0adapt_v6','pitch_f0adapt_v7','pitch_f0adapt','music_scram_familiar','pitch_overlap_v2','pitch_overlap_v3','pitch_adapt_params','pitch_adapt_params_v3','naturalsound'}
        version = '5.3.0';
    otherwise
        version = 'emf-5.1.0';
end
    