function [fwhm, actual_fwhm] = read_smooth(exp,us,runtype,varargin)

switch exp
    case {'tonotopy_monkey','pitch_localizer_monkey','voice_localizer_monkey'}
        if us == 170
            fwhm = 2;
        else
            fwhm = 2.8571; % 1mm
        end
        actual_fwhm = fwhm;
    case {'color_monkey','naturalsound_monkey'}
        fwhm = 2.8571*2; % 2mm
        actual_fwhm = fwhm;
    case {'resting_monkey'}
        fwhm = 2.8571*1; % 2mm
        actual_fwhm = fwhm;
    case {'naturalsound','pitch_resthr_v4','amusia'}
        fwhm = 0.03;
        actual_fwhm = 3;
    case {'pitch_localizer_human','voice_localizer_human','naturalsound-nmf-localizer','tono-pitch-localizer', 'tono-localizer'}
        fwhm = 3;
        actual_fwhm = 3;
    case {'naturalsound-infant'}
        fwhm = 1;
        actual_fwhm = 1;
    otherwise
        error('No valid experiment');
end
