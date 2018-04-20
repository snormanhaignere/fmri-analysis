function bet_func(exp,us,runtype,r,varargin)

% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

datadir = [params('rootdir') exp '/data/'];

% fa values
switch exp
    otherwise
        fafunc = 0.08;
        fastruct = 0.15;
end


niidir = [datadir 'brain/nifti/usub' num2str(us) '/'];
slicetime_file = [niidir runtype '_r' num2str(r) '_LAS_slicetime'];
bet_file = [niidir runtype '_r' num2str(r) '_LAS_bet'];

fsl_version = read_fsl_version(exp);

if ~strcmp(runtype, 'struct')
    mean_func_file = [niidir runtype '_r' num2str(r) '_LAS_mean_func'];
    if ~exist(mean_func_file,'file') || optInputs(varargin, 'overwrite')
        unix(fsl_version,['fslmaths '  slicetime_file ' -Tmean ' mean_func_file]);
    end

if ~exist(bet_file,'file') || optInputs(varargin, 'overwrite')
    unix(fsl_version,['bet2 ' mean_func_file ' mask -f 0.3 -n -m; /usr/share/fsl/5.0/bin/immv mask_mask mask']);
end

[~, ~, ~, ~, ~, ~, ~, nTR] = read_scanparams(exp,us,runtype,varargin{:});
if ~exist([slicetime_file '.nii.gz'],'file')
    unix_fsl(fsl_version, ['slicetimer -i ' mcflirt_file ' --out=' slicetime_file ' -r ' num2str(nTR)]);
end


if strcmp(runtype,'struct');
    favalue = fastruct;
else
    favalue = fafunc;
end

fprintf(['bet ' lasfile ' ' lasfile '_brain' ' -R  -f ' num2str(favalue) ' -g 0 -n -m' '\n']);
unix(['bet ' lasfile ' ' lasfile '_brain' ' -R  -f ' num2str(favalue) ' -g 0 -n -m']);

fprintf(['fslmaths ' lasfile ' -mas ' lasfile '_brain_mask' ' ' lasfile '_brain' '\n']);
unix(['fslmaths ' lasfile ' -mas ' lasfile '_brain_mask' ' ' lasfile '_brain']);


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