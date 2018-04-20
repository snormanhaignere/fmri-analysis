function mcflirt(exp,us,runtype,r,varargin)

% function mcflirt(exp,us,runtype,r,varargin)
% motion corrects data using FSL's mcflirt algorithm

addpath(genpath('/mindhive/nklab/u/svnh/export_fig_v3'));

% fsl version for this experiment
fsl_version = read_fsl_version(exp,varargin{:});

% zero-indexed functional volume to which the run is motion-corrected
reference_frame = read_reference_frame(exp, us, runtype, r, varargin{:});

% preprocessing directory
preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r) '/'];
if ~exist(preprocdir,'dir');
    mkdir(preprocdir);
end

% input file
raw_file = [preprocdir 'raw'];

% volume used for motion-correction
example_func_file = [preprocdir 'example_func'];

% motion corrected file
mcflirt_file = [preprocdir 'motcorr'];

% create reference volume
if ~exist([example_func_file '.nii.gz'], 'file') || optInputs(varargin, 'overwrite');
    unix_fsl(fsl_version, ['fslroi ' raw_file ' ' example_func_file ' ' num2str(reference_frame) ' 1']);
end

% motion correct
if ~exist([mcflirt_file '.nii.gz'],'file') || optInputs(varargin, 'overwrite');
    unix_fsl(fsl_version, ['mcflirt -in ' raw_file ' -out ' mcflirt_file ' -reffile ' example_func_file ' -mats -plots -rmsrel -rmsabs']);
end

% plot rms motion, see motion_summary for more extensive plots
if optInputs(varargin, 'plot')
    fid = fopen([preprocdir 'motcorr_rel.rms']);
    x = textscan(fid,'%f'); fclose(fid);
    if optInputs(varargin, 'monkey')
      relrms = x{1}(:)/2.8571;
    else
      relrms = x{1}(:);
    end
    fprintf('\n\n%s, run %d\n',runtype,r);
    fprintf('Mean: %.3f\n',mean(relrms(:)));
    fprintf('Time points > 0.5 mm: %d\n',sum(relrms(:)>0.5));
    fprintf('Time points > 1 mm: %d\n',sum(relrms(:)>1));
    figure;
    plot(relrms);
    ylim([0 3]);
    xlabel('Scans'); ylabel('Time-Point-to-Time-Point Motion (mm)');
    title(sprintf('%s, run %d',strrep(runtype,'_',' '),r));
    box off;
    export_fig([preprocdir 'relrms.pdf'],'-pdf','-nocrop','-transparent');
end