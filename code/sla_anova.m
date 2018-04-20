function sla_anova(exp,us,runtype,runnum,condlevels,model,varargin)

% function sla(s,runtype,runnum,contr,model,varargin)
%
% main function for running second level analyses
%
% EDIT: in future versions, should add runtype to gfeat and design

analysisdir = [params('rootdir') exp '/analysis/'];
flasubdir = [analysisdir 'fla/usub' num2str(us) '/'];
fsftempdir = [params('rootdir') 'fmri-analysis/code/fsf-templates/sla-anova/'];

fsfsubdir = [analysisdir 'fsf/sla/usub' num2str(us) '/'];
if ~exist(fsfsubdir, 'dir');
    mkdir(fsfsubdir);
end

slasubdir = [analysisdir 'sla/usub' num2str(us) '/'];
if ~exist(slasubdir, 'dir');
    mkdir(slasubdir);
end

fwhm = read_smooth(exp, varargin{:});

x = sprintf('%s',condlevels{:});
if length(x) > 100
    gfeat = [slasubdir runtype '_anova_' DataHash(x) '_r' sprintf('%d',runnum) '_' model '_' num2str(fwhm*100,'%.0f') 'mm'];
else
    gfeat = [slasubdir runtype '_anova_' x '_r' sprintf('%d',runnum) '_' model '_' num2str(fwhm*100,'%.0f') 'mm'];
end

% fsl_version = read_fsl_version(exp,varargin);
fsl_version = '4.1.3';

if exist([gfeat '.gfeat/cope1.feat/stats/zstat1.nii.gz'],'file') && ~optInputs(varargin, 'overwrite');
    return;
end

%% parameters

voxelthresh = '2.3';
clusterthresh = '0.05';

copefiles = cell(length(condlevels),length(runnum));
for r = 1:length(runnum)
    for c = 1:length(condlevels)
        featdir = [flasubdir runtype '_r' num2str(runnum(r)) '_' model  '_' num2str(fwhm*100,'%.0f') 'mm.feat/'];
        fid = fopen([featdir 'contrastnames.txt'],'r');
        tmp = textscan(fid,'%s\n'); fclose(fid); contrastnames = tmp{1};
        ind = strmatch(condlevels{c},contrastnames,'exact');
        if isempty(ind); error('Error: no matching contrast'); end
        copefiles{c,r} = [featdir 'stats/cope' num2str(ind) '.nii.gz'];
    end
end

x = sprintf('%s',condlevels{:});
if length(x) > 100
    fsffile = [fsfsubdir 'design_' runtype '_anova_' DataHash(x) '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm*100,'%.0f') 'mm.fsf'];
else
    fsffile = [fsfsubdir 'design_' runtype '_anova_' x '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm*100,'%.0f') 'mm.fsf'];
end
fidfsf = fopen(fsffile,'w');

%% Design and Contrast Matrix

nruns = length(runnum);
nconds = length(condlevels);
runs_regressors = nan(nconds*nruns, nruns);
for i = 1:nruns
    runs_regressors(:,i) = circshift([ones(nconds,1); zeros(nconds*(nruns-1),1)], (i-1)*nconds);
end

x = [eye(nconds-1); -ones(1,nconds-1)];
condition_regressors = repmat(x, nruns, 1);

design_matrix = [runs_regressors, condition_regressors];
contrast_matrix = [zeros(length(condlevels)-1, length(runnum)), eye(length(condlevels)-1)];

evnames = cell(1,nruns + nconds-1);
for i = 1:nruns
    evnames{i} = ['run' num2str(i)];
end
for i = 1:nconds-1
    evnames{i+nruns} = [condlevels{i} '-' condlevels{end}];
end
contrastnames = evnames(nruns+1:end);

figure(1); clf(1);
subplot(2,1,1);
imagesc(design_matrix);
hold on;
title('Design Matrix');
subplot(2,1,2);
imagesc(contrast_matrix);
title('Contrast Matrix');
drawnow;

%% Misc Parameters

templatefile = [fsftempdir 'part1.txt'];
fidtemp = fopen(templatefile,'r');

clear keyname; clear value;
keyname{1}  = 'WATCHFEAT';        value{1}  = '0';
keyname{2}  = 'OUTDIRECTORY';     value{2}  = gfeat;
keyname{3}  = 'NCOPES';           value{3}  = num2str(numel(copefiles));
keyname{4}  = 'CLUSTERTHRESH';    value{4}  = clusterthresh;
keyname{5}  = 'VOXELTHRESH';      value{5}  = voxelthresh;
keyname{6}  = 'STANDARDBRAIN';    value{6}  = '/software/fsl/fsl-4.1.3/data/standard/MNI152_T1_2mm_brain';
keyname{7}  = 'NEVSREAL';         value{7}  = num2str(length(evnames));
keyname{8}  = 'NEVSORIG';         value{8}  = num2str(length(evnames));
keyname{9}  = 'NCONREAL';         value{9}  = num2str(length(contrastnames));

try
    reptextMANYFAST(fidtemp, fidfsf, keyname, value);
    fclose(fidtemp);
catch
    keyboard;
end
%% Cope files

for i = 1:length(runnum)
    for j = 1:length(condlevels)
        conindex = j + (i-1)*length(condlevels);
        fprintf(fidfsf,'# 4D AVW data or FEAT directory (%d)\n',conindex);
        fprintf(fidfsf,'set feat_files(%d) "%s"\n',conindex,copefiles{j,i});
        fprintf(fidfsf,'\n');
    end
end

%% Boiler Plate

fprintf(fidfsf,'# Add confound EVs text file\n');
fprintf(fidfsf,'set fmri(confoundevs) 0\n');
fprintf('\n\n');

%% Design Matrix

nevs = length(evnames);
for i = 1:nevs
    
    fprintf(fidfsf,'# EV %d title\n',i);
    fprintf(fidfsf,'set fmri(evtitle%d) "%s"\n',i,evnames{i});
    fprintf(fidfsf,'\n');
    fprintf(fidfsf,'# Basic waveform shape (EV %d)\n',i);
    fprintf(fidfsf,'# 0 : Square\n');
    fprintf(fidfsf,'# 1 : Sinusoid\n');
    fprintf(fidfsf,'# 2 : Custom (1 entry per volume)\n');
    fprintf(fidfsf,'# 3 : Custom (3 column format)\n');
    fprintf(fidfsf,'# 4 : Interaction\n');
    fprintf(fidfsf,'# 10 : Empty (all zeros)\n');
    fprintf(fidfsf,'set fmri(shape%d) 2\n',i);
    fprintf(fidfsf,'\n');
    fprintf(fidfsf,'# Convolution (EV %d)\n',i);
    fprintf(fidfsf,'# 0 : None\n');
    fprintf(fidfsf,'# 1 : Gaussian\n');
    fprintf(fidfsf,'# 2 : Gamma\n');
    fprintf(fidfsf,'# 3 : Double-Gamma HRF\n');
    fprintf(fidfsf,'# 4 : Gamma basis functions\n');
    fprintf(fidfsf,'# 5 : Sine basis functions\n');
    fprintf(fidfsf,'# 6 : FIR basis functions\n');
    fprintf(fidfsf,'set fmri(convolve%d) 0\n',i);
    fprintf(fidfsf,'\n');
    fprintf(fidfsf,'# Convolve phase (EV %d)\n',i);
    fprintf(fidfsf,'set fmri(convolve_phase%d) 0\n',i);
    fprintf(fidfsf,'\n');
    fprintf(fidfsf,'# Apply temporal filtering (EV %d)\n',i);
    fprintf(fidfsf,'set fmri(tempfilt_yn%d) 0\n',i);
    fprintf(fidfsf,'\n');
    fprintf(fidfsf,'# Add temporal derivative (EV %d)\n',i);
    fprintf(fidfsf,'set fmri(deriv_yn%d) 0\n',i);
    fprintf(fidfsf,'\n');
    fprintf(fidfsf,'# Custom EV file (EV %d)\n',i);
    fprintf(fidfsf,'set fmri(custom%d) "dummy"\n',i);
    fprintf(fidfsf,'\n');
    
    % no orthogonalization
    for j = 0:nevs
        fprintf(fidfsf,'# Orthogonalise EV %d wrt EV %d\n',i, j);
        fprintf(fidfsf,'set fmri(ortho%d.%d) 0\n',i, j);
        fprintf(fidfsf,'\n');
    end
    
    % design matrix value
    for j = 1:numel(copefiles);
        fprintf(fidfsf,'# Higher-level EV value for EV %d and input %d\n', i, j);
        fprintf(fidfsf,'set fmri(evg%d.%d) %d\n', j, i, design_matrix(j,i));
        fprintf(fidfsf,'\n');
    end
end

%% Boiler Plate

fprintf(fidfsf,'# Setup Orthogonalisation at higher level?\n')
fprintf(fidfsf,'set fmri(level2orth) 0\n\n');

for i = 1:numel(copefiles)
    fprintf(fidfsf,['# Group membership for input ' num2str(i) '\n']);
    fprintf(fidfsf,['set fmri(groupmem.' num2str(i) ') 1\n\n']);
end

fprintf(fidfsf,'# Contrast & F-tests mode\n');
fprintf(fidfsf,'# real : control real EVs\n');
fprintf(fidfsf,'# orig : control original EVs\n');
fprintf(fidfsf,'set fmri(con_mode_old) real\n');
fprintf(fidfsf,'set fmri(con_mode) real\n\n');

%% Contrast Matrix

for i = 1:length(contrastnames)
    fprintf(fidfsf,'# Display images for contrast_real %d\n',i);
    fprintf(fidfsf,'set fmri(conpic_real.%d) 1\n',i);
    fprintf(fidfsf,'\n');
    
    fprintf(fidfsf,'# Title for contrast_real %d\n',i);
    fprintf(fidfsf,'set fmri(conname_real.%d) %s\n',i,contrastnames{i});
    fprintf(fidfsf,'\n');
    
    for j = 1:length(evnames)
        fprintf(fidfsf,'# Real contrast_real vector %d element %d\n',i,j);
        fprintf(fidfsf,'set fmri(con_real%d.%d) %d\n',i,j,contrast_matrix(i,j));
        fprintf(fidfsf,'\n');
    end
    
    fprintf(fidfsf,'# F-test 1 element %d\n',i);
    fprintf(fidfsf,'set fmri(ftest_real1.%d) 1\n',i);
    fprintf(fidfsf,'\n');
end

fprintf(fidfsf,'# Contrast masking - use >0 instead of thresholding?\n');
fprintf(fidfsf,'set fmri(conmask_zerothresh_yn) 0\n');
fprintf(fidfsf,'\n');

%% Boiler Plate

for i = 1:length(contrastnames)+1
    for j = setdiff(1:length(contrastnames)+1,i)
        fprintf(fidfsf,'# Mask real contrast/F-test %d with real contrast/F-test %d?\n',i,j);
        fprintf(fidfsf,'set fmri(conmask%d_%d) 0\n',i,j);
        fprintf(fidfsf,'\n');
    end
end

fprintf(fidfsf,'# Do contrast masking at all?\n');
fprintf(fidfsf,'set fmri(conmask1_1) 0\n');
fprintf(fidfsf,'\n');

templatefile = [fsftempdir 'part2.txt'];
fidtemp = fopen(templatefile,'r');
reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);
fclose(fidfsf);

%%
if exist([gfeat '.gfeat'], 'dir');
    fprintf('Deleting file: %s\n',[gfeat '.gfeat']);
    unix(['rm -r ' gfeat '.gfeat']);
end

unix_fsl(fsl_version,['feat ' fsffile],varargin{:});
