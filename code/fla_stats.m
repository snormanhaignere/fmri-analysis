function fla_stats(exp,us,runtype,r,model,stage,varargin)

% function fla(s,runtype,r,model,stage,varargin)
%
% main script for running a first level analysis with FSL
% input is a subject number, a runtype, and a run number
% as well as a model type (block or event), and a stage (glm or contrast)

analysisdir = [params('rootdir') exp '/analysis/'];
datadir = [params('rootdir') exp '/data/'];
fsftempdir = [params('rootdir') 'fmri-analysis/code/fsf-templates/fla/'];
fsfsubdir = [analysisdir 'fsf/fla/usub' num2str(us) '/'];
flasubdir = [analysisdir 'fla/usub' num2str(us) '/'];
stfdir = [analysisdir 'stf/'];
fsl_version = read_fsl_version(exp,varargin{:});

if ~exist(fsfsubdir, 'dir'); mkdir(fsfsubdir); end
if ~exist(flasubdir, 'dir'); mkdir(flasubdir); end

%% general parameters

[~, ~, TR, ~, ~, ~, ~, nTR] = read_scanparams(exp,us,runtype,'run',r);
nvolumes = num2str(nTR);
trsecs = num2str(TR);
slicetiming = '1'; % 0 = none, 1 = regular up (0,1,2...), 5 = interleaved (0,2,4...1,3,5...)
hpcutoff = num2str(params('hpcutoff'));
fwhm = read_smooth(exp, varargin{:});

motionparams = '1'; 
regnonlinear = '0';
tempderiv = params('tempderiv',model);

% whether or not to automatically bring up log file
watchfeat = '0';
if optInputs(varargin,'watchfeat');
    watchfeat = '1';
end

featdir = [flasubdir runtype '_r' num2str(r) '_' model  '_' num2str(fwhm*100,'%.0f') 'mm'];
if optInputs(varargin, 'input_fname')
    funcvolume = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/' varargin{optInputs(varargin, 'input_fname')+1} '.nii.gz'];
elseif optInputs(varargin, 'detrend')
    funcvolume = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/smooth' num2str(fwhm*100,'%.0f') 'mm_intnorm_detrend1.nii.gz'];
else
    funcvolume = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/smooth' num2str(fwhm*100,'%.0f') 'mm_intnorm_hpfilt' num2str(params('hpcutoff')) '.nii.gz'];
end

if optInputs(varargin, 'functional_reference_volume')
    allruns = read_runs(exp, us, runtype, varargin{:});
    anatomicalvolume = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r' num2str(allruns(1)) '/example_func.nii.gz'];
    dof_anatomical = '12';
else
    anatomicalvolume = [analysisdir 'preprocess/usub' num2str(us) '/struct_r1/brain.nii.gz'];
    dof_anatomical = '6';
end

if any(strcmp(stage, {'glm','glm-noreg','film_gls'})) && exist([featdir '.feat/stats/pe1.nii.gz'],'file') && ~optInputs(varargin,'overwrite');
    return;
end

if any(strcmp(stage, {'reg'})) && exist([featdir '.feat/reg/example_func2highres.mat'],'file') && ~optInputs(varargin,'overwrite');
    return;
end

if any(strcmp(stage, {'contrasts'})) && exist([featdir '.feat/stats/zstat1.nii.gz'],'file') && ~optInputs(varargin,'overwrite');
    return;
end

%% setup evs
conds = read_conditions(exp,us,runtype,'run',r);
nevs = length(conds);
evname = cell(1,nevs);
evfile = cell(1,nevs);
for c = 1:nevs
    evname{c} = conds{c};
    if optInputs(varargin, 'custom_stf')
        evfile{c} = [stfdir 'usub' num2str(us) '/' runtype '_r' num2str(r) '_' conds{c} '_custom.stf'];
    else
        evfile{c} = [stfdir 'usub' num2str(us) '/' runtype '_r' num2str(r) '_' conds{c} '.stf'];
    end
end

% %% polynomial regressors
% if optInputs(varargin, 'detrend')
%     polynames = {'linear'};
%     for i = 1:length(polynames)
%         evname = [evname, polynames{i}];
%         polyfile = [stfdir 'usub' num2str(us) '/' runtype '_r'  num2str(r) '_' polynames{i} '.stf'];
%         evfile = [evfile, polyfile];
%         t = (1:nTR)' - floor(nTR/2);
%         fid = fopen(polyfile,'w');
%         fprintf(fid,repmat('%.6f\n',1,nTR),t.^i);
%         fclose(fid);
%     end
%     nevs = length(evname);
% end

%% setup contrasts

contrastnames = read_contrasts(exp,us,runtype,model,'run',r);
contrasts = zeros(length(contrastnames),nevs);

for c = 1:length(contrastnames);
    
    % an ev group name is a name for a set of conditions that can be used in a contrast
    evgroup = regexp(contrastnames{c},'_vs_','split');
    
    if length(evgroup)==1
        
        evs2find1 = evgroup2evname(exp,us,evgroup{1},model,'run',r);
        evinds1 = strmatchMANY( evs2find1,  evname );
        
        if isempty(evinds1) || length(evs2find1) ~= length(evinds1)
            fprintf(['Error in fla_stats: Failed to match contrasts to evs: ' contrastnames{c} '\n']); drawnow;
            keyboard;
        end
        
        contrasts(c,evinds1) = 1;
        
    else
        
        evs2find1 = evgroup2evname(exp,us,evgroup{1},model,'run',r);
        evs2find2 = evgroup2evname(exp,us,evgroup{2},model,'run',r);
        evinds1 = strmatchMANY( evs2find1,  evname );
        evinds2 = strmatchMANY( evs2find2,  evname );
        
        if isempty(evinds1) || isempty(evinds2) || length(evs2find1) ~= length(evinds1) || length(evs2find2) ~= length(evinds2)
            fprintf(['Error in fla_stats: Failed to match contrasts to evs: ' contrastnames{c} '\n']); drawnow;
            keyboard;
        end
        
        contrasts(c,evinds1) = length(evinds2);
        contrasts(c,evinds2) = -length(evinds1);
        
    end
end

%% complicated contrasts

[~,~,cc] = read_contrasts(exp,us,runtype,model,'run',r);
if ~isempty(cc)
    ccmat = nan(size(cc,1),nevs);
    
    for i = 1:size(cc,1)
        weightcell = cc{i,2};
        weightmat = zeros(1,nevs);
        for j = 1:size(weightcell,1)
            x = strcmp(weightcell{j,1},evname);
            weightmat(x) = weightmat(x) + weightcell{j,2};
        end
        ccmat(i,:) = weightmat;
    end
    
    contrasts = [contrasts; ccmat];
    contrastnames = [contrastnames, cc(:,1)'];
end

%% F-Tests
[~,~,~,~,ftests] = read_contrasts(exp,us,runtype,model,'run',r);
if ~isempty(ftests)
    ftmat = zeros(size(ftests,1),length(contrastnames));
    for i = 1:size(ftests,1)
        for j = 1:length(ftests{i,2})
            ftmat(i,strcmp(ftests{i,2}{j},contrastnames)) = 1;
        end
    end
end

%% open fsf file
fsffile = [fsfsubdir  'design_' runtype '_r' num2str(r) '_' model '_' num2str(fwhm*100,'%.0f') 'mm_' stage '.fsf'];
fidfsf = fopen(fsffile,'w');
if fidfsf == -1
    error('Error in fla_stats: unable to open %s\n',fsffile);
end

%% generate fsf file in parts, part 1, initial setup

fidtemp = fopen([fsftempdir 'part1.txt'],'r');
if fidtemp == -1
    error('Error in fla_stats: unable to open %s\n',[fsftempdir 'part1.txt']);
end

if strcmp(tempderiv,'1')
    nevsreal = nevs*2;
elseif strcmp(tempderiv,'0')
    nevsreal = nevs;
end

if strcmp(stage, 'glm')
    inputvolume = funcvolume;
    stagenum = '2';
    filterfunc = '0';
    mainstats = '1';
    regimages = '1';
    poststats = '0';
elseif strcmp(stage, 'contrasts')
    inputvolume = [featdir '.feat'];
    stagenum = '4';
    filterfunc = '0';
    mainstats = '0';
    regimages = '0';
    poststats = '1';
elseif strcmp(stage, 'glm-noreg')
    inputvolume = funcvolume;
    stagenum = '2';
    filterfunc = '0';
    mainstats = '1';
    regimages = '0';
    poststats = '0';
elseif strcmp(stage, 'film_gls')
    inputvolume = funcvolume;
    stagenum = '-1';
    filterfunc = '0';
    mainstats = '0';
    regimages = '0';
    poststats = '0';
elseif strcmp(stage, 'reg')
    inputvolume = [featdir '.feat'];
    stagenum = '0';
    filterfunc = '0';
    mainstats = '0';
    regimages = '1';
    poststats = '0';
else
    error('No matching stage in fla_stats.m');
end

switch fsl_version
    case '5.0'
        mni_standard = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz';
    case '4.1.3'
        mni_standard = '/software/fsl/fsl-4.1.3/data/standard/MNI152_T1_2mm_brain.nii.gz';
    otherwise
        error('No valid fsl version');
end

if optInputs(varargin, 'subject_specific_standard')
    mni_standard = anatomicalvolume;
end

if optInputs(varargin, 'functional_reference_volume')
    allruns = read_runs(exp, us, runtype, varargin{:});
    mni_standard = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(allruns(1)) '/example_func.nii.gz'];
end


if optInputs(varargin, 'MION')
    gamma_sigma = 6;
    gamma_delay = 12;
else
    gamma_sigma = 3;
    gamma_delay = 6;
end

if optInputs(varargin, '180')
    search_angle = '180';
elseif optInputs(varargin, '30')
    search_angle = '30';
else
    search_angle = '90';
end

clear keyname; clear value;
keyname{1}  = 'STAGENUM';         value{1}  = stagenum;
keyname{2}  = 'WATCHFEAT';        value{2}  = watchfeat;
keyname{3}  = 'OUTDIRECTORY';     value{3}  = featdir;
keyname{4}  = 'TRSECS';           value{4}  = trsecs;
keyname{5}  = 'NVOLUMES';         value{5}  = nvolumes;
keyname{6}  = 'FILTERFUNC';       value{6}  = filterfunc;
keyname{7}  = 'SMOOTHWIDTH';      value{7}  = num2str(fwhm);
keyname{8}  = 'MAINSTATS';        value{8}  = mainstats;
keyname{9}  = 'MOTIONPARAMS';     value{9}  = motionparams;
keyname{10} = 'NEVSORIG';         value{10} = num2str(nevs);
keyname{11} = 'NEVSREAL';         value{11} = num2str(nevsreal);
keyname{12} = 'NCONORIG';         value{12} = num2str(length(contrastnames));
keyname{13} = 'NCONREAL';         value{13} = num2str(length(contrastnames));
keyname{14} = 'POSTSTATS';        value{14} = poststats;
keyname{15} = 'REGIMAGES';        value{15} = regimages;
keyname{16} = 'DOFANATOMICAL';    value{16} = dof_anatomical;
keyname{17} = 'DOFSTANDARD';      value{17} = '12';
keyname{18} = 'STANDARDBRAIN';    value{18} = mni_standard;
keyname{19} = 'INPUTVOLUME';      value{19} = inputvolume;
keyname{20} = 'ANATOMICALVOLUME'; value{20} = anatomicalvolume;
keyname{21} = 'REGNONLINEAR';     value{21} = regnonlinear;
keyname{22} = 'SLICETIMING';      value{22} = slicetiming;
keyname{23} = 'HPCUTOFF';         value{23} = hpcutoff;


if optInputs(varargin, 'confound_regressors') ||  strcmp(exp, 'sepi') && us == 1 && any(strcmp(runtype, {'continuous', 'continuous_sepi'})) && any(r == [1 2])
    keyname{24} = 'CONFOUNDREGRESSOR';value{24} = '1';
else
    keyname{24} = 'CONFOUNDREGRESSOR';value{24} = '0';
end
keyname{25} = 'SEARCHANGLE';      value{25} = search_angle;
keyname{26} = 'NFTESTSREAL';      value{26} = num2str(length(ftests));
keyname{27} = 'NFTESTSORIG';      value{27} = num2str(length(ftests));

if strcmp(exp, 'sepi') && us == 1 && any(strcmp(runtype, {'continuous', 'continuous_sepi'})) && any(r == [1 2])
    fprintf(fidfsf, '\n\n');
    fprintf(fidfsf, '# Confound EVs text file for analysis 1\n');
    fprintf(fidfsf, 'set confoundev_files(1) "/mindhive/nklab/u/svnh/sepi/analysis/fla/usub1/tp1-2.txt"\n\n');
elseif optInputs(varargin, 'confound_regressors')
    confound_file = [stfdir 'usub' num2str(us) '/' runtype '_r' num2str(r) '_confound_regressors.stf'];
    fprintf(fidfsf, '\n\n');
    fprintf(fidfsf, '# Confound EVs text file for analysis 1\n');
    fprintf(fidfsf, ['set confoundev_files(1) "' confound_file '"\n\n']);
end

reptextMANYFAST(fidtemp, fidfsf, keyname, value);
fclose(fidtemp);

%% part 2,3, glm

fidtemp1 = fopen([fsftempdir 'part2.txt'], 'r');
fidtemp2 = fopen([fsftempdir 'part3.txt'], 'r');

clear keyname; clear value;
for e1 = 1:nevs
    keyname{1}{1} = 'EVINDEX';   value{1}{1} = num2str(e1); %#ok<*AGROW>
    keyname{1}{2} = 'EVNAME';    value{1}{2} = evname{e1};
    if optInputs(varargin, 'custom_stf')% || (optInputs(varargin, 'detrend') && any(strcmp(evname{e1}, polynames)))
        keyname{1}{3} = 'STFTYPE';   value{1}{3} = '2'; % 1-column format
        keyname{1}{4} = 'CONVTYPE';  value{1}{4} = '0'; % no convolution
    else
        keyname{1}{3} = 'STFTYPE';   value{1}{3} = '3'; % 3-column format
        keyname{1}{4} = 'CONVTYPE';  value{1}{4} = '2'; % convolve with gamma
    end
    keyname{1}{5} = 'EVFILE';     value{1}{5} = evfile{e1};
    keyname{1}{6} = 'TEMPDERIV';  value{1}{6} = tempderiv;
    keyname{1}{7} = 'GAMMASIGMA'; value{1}{7} = num2str(gamma_sigma);
    keyname{1}{8} = 'GAMMADELAY'; value{1}{8} = num2str(gamma_delay);
    
    reptextMANYFAST(fidtemp1, fidfsf, keyname{1}, value{1});
    fseek(fidtemp1,0,-1);
    
    for e2 = 0:nevs
        keyname{2}{1} = 'EVINDEX1';  value{2}{1} = num2str(e1);
        keyname{2}{2} = 'EVINDEX2';  value{2}{2} = num2str(e2);
        keyname{2}{3} = 'VALUE';     value{2}{3} = '0';
        reptextMANYFAST(fidtemp2, fidfsf, keyname{2}, value{2});
        fseek(fidtemp2,0,-1);
    end
end

fclose(fidtemp1);
fclose(fidtemp2);

%% part 4
fidtemp = fopen([fsftempdir 'part4.txt'],'r');
reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);

%% part 5,6, real contrasts
fidtemp1 = fopen([fsftempdir 'part5.txt'], 'r');
fidtemp2 = fopen([fsftempdir 'part6.txt'], 'r');

clear keyname; clear value;
for c = 1:size(contrasts,1)
    keyname{1}{1} = 'CTINDEX';   value{1}{1} = num2str(c); %#ok<*AGROW>
    keyname{1}{2} = 'CTNAME';    value{1}{2} = contrastnames{c};
    
    reptextMANYFAST(fidtemp1, fidfsf, keyname{1}, value{1});
    fseek(fidtemp1,0,-1);
    
    for e = 1:nevsreal
        keyname{2}{1} = 'CTINDEX';  value{2}{1} = num2str(c);
        keyname{2}{2} = 'EVINDEX';  value{2}{2} = num2str(e);
        
        if strcmp(tempderiv,'1') % adding temporal derivatives doubles number of regressors
            if mod(e,2)==1
                keyname{2}{3} = 'VALUE'; value{2}{3} = num2str(contrasts(c,(e+1)/2));
            else
                keyname{2}{3} = 'VALUE'; value{2}{3} = '0';
            end
        elseif strcmp(tempderiv,'0')
            keyname{2}{3} = 'VALUE'; value{2}{3} = num2str(contrasts(c,e));
        end
        
        reptextMANYFAST(fidtemp2, fidfsf, keyname{2}, value{2});
        fseek(fidtemp2,0,-1);
    end
    
    for f = 1:size(ftests,1)
        fprintf(fidfsf, ['# F-test ' num2str(f) ' element ' num2str(c) '\n']);
        fprintf(fidfsf, ['set fmri(ftest_real' num2str(f) '.' num2str(c) ') ' num2str(ftmat(f,c)) '\n\n']);
    end
end

fclose(fidtemp1);
fclose(fidtemp2);

%% part 7,8, original contrasts

fidtemp1 = fopen([fsftempdir 'part7.txt'], 'r');
fidtemp2 = fopen([fsftempdir 'part8.txt'], 'r');

clear keyname; clear value;
for c = 1:size(contrasts,1)
    keyname{1}{1} = 'CTINDEX';   value{1}{1} = num2str(c); %#ok<*AGROW>
    keyname{1}{2} = 'CTNAME';    value{1}{2} = contrastnames{c};
    
    reptextMANYFAST(fidtemp1, fidfsf, keyname{1}, value{1});
    fseek(fidtemp1,0,-1);
    
    for e = 1:nevs
        keyname{2}{1} = 'CTINDEX';  value{2}{1} = num2str(c);
        keyname{2}{2} = 'EVINDEX';  value{2}{2} = num2str(e);
        keyname{2}{3} = 'VALUE';    value{2}{3} = num2str(contrasts(c,e));
        
        reptextMANYFAST(fidtemp2, fidfsf, keyname{2}, value{2});
        fseek(fidtemp2,0,-1);
    end
    
    for f = 1:size(ftests,1)
        fprintf(fidfsf, ['# F-test ' num2str(f) ' element ' num2str(c) '\n']);
        fprintf(fidfsf, ['set fmri(ftest_orig' num2str(f) '.' num2str(c) ') ' num2str(ftmat(f,c)) '\n\n']);
    end
    
end

fclose(fidtemp1);
fclose(fidtemp2);


%% part 9

fidtemp = fopen([fsftempdir 'part9.txt'],'r');
reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);

%% part10, f-test

fidtemp = fopen([fsftempdir 'part10.txt'], 'r');

clear keyname; clear value;
for c1 = 1:size(contrasts,1)
    
    for c2 = setdiff(1:size(contrasts,1),c1)
        keyname{1} = 'CTINDEX1';  value{1} = num2str(c1);
        keyname{2} = 'CTINDEX2';  value{2} = num2str(c2);
        
        reptextMANYFAST(fidtemp, fidfsf, keyname, value);
        fseek(fidtemp,0,-1);
    end
end

fclose(fidtemp);

%% part 11

fidtemp = fopen([fsftempdir 'part11.txt'],'r');
reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);
fclose(fidfsf);

%% Just run film_gls

if strcmp(stage,'film_gls')
    if ~exist([featdir '.feat'],'dir');
        mkdir([featdir '.feat']);
    end
    
    unix(['cp '  fsffile ' ' featdir '.feat/design.fsf']);
    unix(['cp '  funcvolume ' ' featdir '.feat/filtered_func_data.nii.gz']);
    
    if optInputs(varargin, 'confound_regressors')
        confound_file = [stfdir 'usub' num2str(us) '/' runtype '_r' num2str(r) '_confound_regressors.stf'];
        unix_fsl(fsl_version, ['feat_model ' featdir './design ' confound_file]);
    else
        unix_fsl(fsl_version, ['feat_model ' featdir '.feat/design ']);
    end
        
    unix_fsl('4.1', ['film_gls -rn ' featdir '.feat/stats ' featdir '.feat/filtered_func_data.nii.gz '  featdir '.feat/design.mat 0.001']);
    delete([featdir '.feat/filtered_func_data.nii.gz']);
    return;
end


mask_file = [featdir '.feat/mask.nii.gz'];
mean_func = [featdir '.feat/mean_func.nii.gz'];
unix_fsl(fsl_version,['fslmaths ' funcvolume ' -Tmin -bin ' mask_file ' -odt char']);
unix_fsl(fsl_version,['fslmaths ' funcvolume ' -Tmean ' mean_func]);

%% Create files for registration if stats has not already been run
if strcmp(stage,'reg') && ~exist([featdir '.feat'],'dir');
    mkdir([featdir '.feat']);
end

if strcmp(stage,'reg')
    example_func_file = [featdir '.feat/example_func.nii.gz'];
    example_func_raw_file = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/example_func_bet.nii.gz'];
    if ~exist(example_func_file,'file')
        copyfile(example_func_raw_file,example_func_file,'f');
    end
end

%% call feat!

if ~optInputs(varargin,'norun')
    
    if optInputs(varargin,'overwrite') && any(strcmp(stage, {'glm','glm-noreg'})) && exist([featdir '.feat'], 'dir');
        fprintf('Deleting file: %s\n',[featdir '.feat']);
        unix(['rm -r ' featdir '.feat']);
    end
    
    fprintf(['feat ' fsffile '\n']);
    %     if strcmp(exp, 'pitch_localizer_human')
    %         unix_fsl('4.3.0', ['feat ' fsffile]);
    %     else
    
    unix_fsl('4.1', ['feat ' fsffile]);
    %     end
    
end
%% write extra files to feat directory

while true
    if exist([featdir '.feat'],'dir');
        fid = fopen([featdir '.feat/evname.txt'],'w');
        for j = 1:length(evname);
            fprintf(fid, '%s\n', evname{j});
        end
        fclose(fid);
        
        fid = fopen([featdir '.feat/contrastnames.txt'],'w');
        for j = 1:length(contrastnames);
            fprintf(fid, '%s\n', contrastnames{j});
        end
        fclose(fid);
        break;
    end
end

%% Move example_func files

example_func_file = [featdir '.feat/example_func.nii.gz'];
example_func_smooth = [featdir '.feat/example_func_smooth.nii.gz'];
example_func_raw_file = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/example_func_bet.nii.gz'];

if exist(example_func_file, 'file') && ~exist(example_func_smooth,'file')
    movefile(example_func_file,example_func_smooth);
end

if exist(example_func_smooth,'file') && exist(example_func_raw_file,'file')
    copyfile(example_func_raw_file,example_func_file,'f');
end

% example_func_file = [featdir '.feat/reg/example_func.nii.gz'];
% example_func_smooth = [featdir '.feat/reg/example_func_smooth.nii.gz'];
% example_func_raw_file = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r'   num2str(r) '/example_func.nii.gz'];
%
% if exist(example_func_file, 'file') && ~exist(example_func_smooth,'file')
%     movefile(example_func_file,example_func_smooth);
% end
%
% if exist(example_func_smooth,'file') && ~exist(example_func_file,'file') && exist(example_func_raw_file,'file')
%     copyfile(example_func_raw_file,example_func_file);
% end

%% scraps
% if s <= 2
%     error('This should be checked before being run again');
%
%     contrastnames = read
% else
%     %         contrastnames = [evname,'orchestra_vs_nonwords','nonwords_vs_orchestra','songs_vs_lyrics','lyrics_vs_songs', ...
%     %             'mel-intact_vs_mel-scram', 'mel-scram_vs_mel-intact', 'mel-intact_vs_mel-interp',  'mel-interp_vs_mel-intact', 'drums-intact_vs_drums-scram', 'drums-scram_vs_drums-intact', ...
%     %             'env-pitch_vs_env-nopitch', 'env-nopitch_vs_env-pitch', 'mel-intact_vs_mel-noise', 'mel-noise_vs_mel-intact', 'animals-pitch_vs_nonwords','nonwords_vs_animals-pitch',...
%     %             'sentences_vs_nonwords', 'nonwords_vs_sentences', ...
%     %             'orchestra_vs_animals-pitch', 'orchestra_vs_env-pitch', 'orchestra_vs_env-nopitch', ...
%     %             'drums-intact_vs_nonwords', 'drums-intact_vs_animals-pitch', 'drums-intact_vs_env-pitch', 'drums-intact_vs_env-nopitch'];
%
%     contrastnames = read_contrasts(s,model);
%
%
%     %         contrastnames = [evname, 'orchestra_vs_nonwords', 'orchestra_vs_animals-pitch', 'mel-intact_vs_mel-scram', 'mel-intact_vs_mel-interp', 'songs_vs_lyrics', ...
%     %             'env-pitch_vs_env-nopitch', 'mel-scram_vs_mel-noise',...
%     %             'env-nopitch_vs_sentences',...
%     %             'animals-pitch_vs_env-pitch', 'animals-pitch_vs_orchestra', 'nonwords_vs_animals-pitch', 'nonwords_vs_orchestra', 'sentences_vs_nonwords', ...
%     %             'drums-intact_vs_drums-scram', 'drums-intact_vs_nonwords', 'drums-intact_vs_animals-pitch', 'drums-intact_vs_env-pitch'];
%
% end




