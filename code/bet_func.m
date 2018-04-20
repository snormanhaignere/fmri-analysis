function bet_func(exp,us,runtype,r,varargin)

% function bet_func(exp,us,runtype,r,varargin)
% computes brain mask using FSL's bet2 algorithm

% main parameters higher fa values result in more brain being removed
favalue = read_bet_fa(exp, us, runtype, r, varargin{:});

% version of fsl to use
fsl_version = read_fsl_version(exp,varargin{:});

% preprocessing directory
preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

if optInputs(varargin, 'noslicetimecorrection')
    input_file = [preprocdir 'motcorr'];
else
    input_file = [preprocdir 'slicetimecorr'];
end
mean_func_file = [preprocdir 'mean_func'];
bet_file = [preprocdir 'bet'];
mask_file = [preprocdir 'brain_mask'];
thresh_file = [preprocdir 'brain_thresh'];
mode_file = [preprocdir 'mode.txt'];

% uses mean of functional file as input
if ~exist([bet_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version,['fslmaths '  input_file ' -Tmean ' mean_func_file]);
    unix_fsl(fsl_version,['bet2 ' mean_func_file ' ' bet_file ' -f ' num2str(favalue) ' -m']);
end

% post-processing
if ~exist([mask_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite')
    
    % create mask
    y = unix_fsl(fsl_version,['fslstats ' bet_file ' -p 2 -p 98']);
    thresh_vals = str2double(strtrim(regexp(strtrim(y),' ','split')));
    unix_fsl(fsl_version, ['fslmaths ' bet_file ' -thr ' num2str(thresh_vals(2)/10) ' -Tmin -bin ' mask_file ' -odt char']);
    
    % store mode of input data, could be used for intensity normalization
    mode = unix_fsl(fsl_version, ['fslstats ' input_file ' -k ' mask_file ' -p 50']);
    fid = fopen(mode_file,'w');
    fprintf(fid, '%f', str2double(strtrim(mode))); fclose(fid);
    
    % dilate mask
    unix_fsl(fsl_version, ['fslmaths ' mask_file ' -dilF ' mask_file]);
    
    % rethreshold slice time corrected brain
    unix_fsl(fsl_version, ['fslmaths ' input_file ' -mas ' mask_file ' ' thresh_file]);
    
    % recompute mean func file
    unix_fsl(fsl_version, ['fslmaths ' thresh_file ' -Tmean ' mean_func_file]);
    
end

% quick overwrite to check that new fa value improves brain extraction
if optInputs(varargin, 'quick_overwrite')
    % bet2
    unix_fsl(fsl_version,['bet2 ' mean_func_file ' ' bet_file ' -f ' num2str(favalue) ' -m']);
    
    % create mask
    y = unix_fsl(fsl_version,['fslstats ' bet_file ' -p 2 -p 98']);
    thresh_vals = str2double(strtrim(regexp(strtrim(y),' ','split')));
    unix_fsl(fsl_version, ['fslmaths ' bet_file ' -thr ' num2str(thresh_vals(2)/10) ' -Tmin -bin ' mask_file ' -odt char']);
    
    % dilate mask
    unix_fsl(fsl_version, ['fslmaths ' mask_file ' -dilF ' mask_file]);
end

% apply mask to example_func file
example_func_file = [preprocdir 'example_func'];
example_func_bet_file = [preprocdir 'example_func_bet'];
if ~exist([example_func_bet_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite')
    unix_fsl(fsl_version, ['fslmaths ' example_func_file ' -mas ' mask_file ' ' example_func_bet_file]);
end


%% Old Code


% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

% fa values
% switch exp
%     %     case 'pitch_overlap_v3'
%     %         switch runtype
%     %             case 'localizer'
%     %                 favalue = 0.15;
%     %             case 'overlap_v3'
%     %                 favalue = 0.15;
%     %             otherwise
%     %                 error('No valid runtype');
%     %         end
%     case 'tonotopy_monkey'
%         if us == 157
%             favalue = 0.1;
%         else
%             favalue = 0.5;
%         end
%     otherwise
%         favalue = 0.15;
% end
% if ~exist(mask_file,'file') || optInputs(varargin, 'overwrite')
%     y = unix_fsl(fsl_version,['fslstats ' bet_file ' -p 2 -p 98']);
%     thresh_vals = str2double(strtrim(regexp(strtrim(y),' ','split')));
%     unix_fsl(fsl_version, ['fslmaths ' bet_file ' -thr ' num2str(thresh_vals(2)/10) ' -Tmin -bin ' mask_file ' -odt char']);
%     mode = unix_fsl(fsl_version, ['fslstats ' slicetime_file ' -k ' mask_file ' -p 50']);
%     mode_file = [niidir runtype '_r' num2str(r) '_LAS_mode.txt'];
%     fid = fopen(mode_file,'w');
%     fprintf(fid, '%f', str2double(strtrim(mode))); fclose(fid);
%     unix_fsl(fsl_version, ['fslmaths ' mask_file ' -dilF ' mask_file]);
% end

% /usr/share/fsl/5.0/bin/fslmaths prefiltered_func_data_bet -thr 99.4399536 -Tmin -bin mask -odt char
%
% /usr/share/fsl/5.0/bin/fslstats prefiltered_func_data_st -k mask -p 50
% 559.589661
%
% /usr/share/fsl/5.0/bin/fslmaths mask -dilF mask


% if exist([lasfile '_brain.nii.gz'],'file')
%     fprintf('Deleting file: %s\n',[lasfile '_brain.nii.gz']);
%     delete([lasfile '_brain.nii.gz']);
% end
%
% if exist([lasfile '_brain_mask.nii.gz'],'file')
%     fprintf('Deleting file: %s\n',[lasfile '_brain_mask.nii.gz']);
%     delete([lasfile '_brain_mask.nii.gz']);
% end


% scriptsdir = [pwd '/'];
% datadir = strrep(scriptsdir, 'scripts/', 'data/');
%
% betcommand = '/usr/local/fsl/bin/bet';
% rawfile = [datadir 'brain/func/raw/sub' num2str(usubs) '/250000-' num2str(runorders(runnum))];
% betfile = [rawfile '_brain'];
%
% favalue = 0.5;
%
% cd('/usr/local/fsl')
% unix([betcommand ' ' rawfile ' ' betfile ' -R -f ' num2str(favalue) ' -g 0' ])
% cd(scriptsdir);