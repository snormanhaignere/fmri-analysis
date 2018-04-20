function P = freesurfer_roi_stats(exp, us, roi, hemi, varargin)

% function P = freesurfer_roi_stats(exp, us, roi, hemi, varargin)
% 
% Wrapper function for the freesurfer command "mris_anatomical_stats".
% Calculates a variety of useful statistics for a fixed ROI.
% 
% -- Example --
% exp = 'amusia';
% us = 45;
% roi = 'stp';
% hemi = 'rh';
% P = freesurfer_roi_stats(exp, us, roi, hemi)
% 
% Last edited by Sam Norman-Haignere on 2015-06-18

freesurfer_version = read_freesurfer_version(exp);

% subject id
if isempty(us)
    subjid = exp;
else
    subjid = [exp '_us' num2str(us)];
end

label_file = [hemi '.' roi '.label'];
label_stats_file = [params('rootdir') 'freesurfer/' subjid '/label/' hemi '.' roi '.stats'];

% use mris_anatomical_stats to calculate roi stats
if ~exist(label_stats_file,'file') || optInputs(varargin, 'overwrite')
    unix_freesurfer_version(freesurfer_version,['mris_anatomical_stats -l ' label_file ' -f ' label_stats_file ' ' subjid ' ' hemi]);
end

% parse the text file
clear P;
fid = fopen(label_stats_file,'r');
while 1

    textline = fgetl(fid);
    if isnumeric(textline) && textline == -1
        break;
    end
    
    if ~strcmp(textline(1),'#')
        split_textline = regexp(textline,'\s+','split');
        P.nvert = str2double(split_textline{2});
        P.surfarea = str2double(split_textline{3});
        P.grayvol = str2double(split_textline{4});
        P.avthickness = str2double(split_textline{5});
        P.stdthickness = str2double(split_textline{6});
        P.meancurv = str2double(split_textline{7});
        P.gauscurv = str2double(split_textline{8});
        P.foldindex = str2double(split_textline{9});
        P.curvindex = str2double(split_textline{10});
    end
end
fclose all;



