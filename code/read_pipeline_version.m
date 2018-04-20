function pipeline_version = read_pipeline_version(exp, varargin)

switch exp
    case {'tono-localizer', 'tono-pitch-localizer','voice_localizer_human','amusia','pitch_resthr_v4','color_monkey','naturalsound_monkey','pitch_localizer_monkey','voice_localizer_monkey','pitch_localizer_human','tonotopy_monkey','pitch_overlap_v3','naturalsound','pitch_f0adapt_v5','pitch_f0adapt_v6','pitch_f0adapt_v7'}
        pipeline_version = 2;
    otherwise
        pipeline_version = 1;
end