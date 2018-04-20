function fla(exp,us,runtype,r,model,stage,varargin)

% function fla(s,runtype,r,model,stage,varargin)
%
% main script for running a first level analysis with FSL
% input is a subject number, a runtype, and a run number
% as well as a model type (block or event), and a stage (glm or contrast)

analysisdir = [params('rootdir') exp '/analysis/'];
datadir = [params('rootdir') exp '/data/'];
fsftempdir = [pwd '/fsf-templates/fla/'];
fsfsubdir = [analysisdir 'fsf/fla/usub' num2str(us) '/'];
flasubdir = [analysisdir 'fla/usub' num2str(us) '/'];
stfdir = [analysisdir 'stf/'];
fsl_version = read_fsl_version(exp);

if ~exist(fsfsubdir, 'dir'); mkdir(fsfsubdir); end
if ~exist(flasubdir, 'dir'); mkdir(flasubdir); end

%% general parameters

[~, ~, TR, ~, ~, ~, ~, nTR] = read_scanparams(exp,us,runtype,'run',r);
nvolumes = num2str(nTR);
trsecs = num2str(TR);
slicetiming = '1'; % 0 = none, 1 = regular up (0,1,2...), 5 = interleaved (0,2,4...1,3,5...)
hpcutoff = '250';

motionparams = '1';
regnonlinear = '0';
tempderiv = params('tempderiv',model);

% whether or not to automatically bring up log file
watchfeat = '0';
if optInputs(varargin,'watchfeat');
    watchfeat = '1';
end

% whether or not to overwrite feat, only relevant for glm stage
overwrite = false;
if optInputs(varargin,'overwrite');
    overwrite = true;
end

featdir = [flasubdir runtype '_r' num2str(r) '_' model  '_' params('smooth') 'mm'];
funcvolume = [datadir 'brain/nifti/usub' num2str(us) '/' runtype '_r'   num2str(r) '_LAS_brain'];
anatomicalvolume = [datadir 'brain/nifti/usub' num2str(us) '/struct_r1_LAS_brain'];

if strcmp(stage, 'glm') && exist([featdir '.feat'],'dir') && ~optInputs(varargin,'overwrite');
    return;
end

%% setup evs

conds = read_conditions(exp,us,runtype);
nevs = length(conds);
evname = cell(1,nevs);
evfile = cell(1,nevs);
for c = 1:nevs
    evname{c} = conds{c};
    evfile{c} = [stfdir 'usub' num2str(us) '/' runtype '_r' num2str(r) '_' conds{c} '.stf'];
end

%% setup contrasts

contrastnames = read_contrasts(exp,us,runtype,model);
contrasts = zeros(length(contrastnames),nevs);

for c = 1:length(contrastnames);
    
    % an ev group name is a name for a set of conditions that can be used in a contrast
    evgroup = regexp(contrastnames{c},'_vs_','split');
    
    if length(evgroup)==1
        
        evinds = strmatchMANY( evgroup2evname(exp,us,evgroup{1},model),  evname );
        if isempty(evinds); error(['Failed to match contrasts to evs: ' contrastnames{c}]); end
        
        contrasts(c,evinds) = 1;
        
    else
        
        evinds1 = strmatchMANY( evgroup2evname(exp,us,evgroup{1},model),  evname );
        evinds2 = strmatchMANY( evgroup2evname(exp,us,evgroup{2},model),  evname );
        
        if isempty(evinds1) || isempty(evinds2); error(['Failed to match contrasts to evs: ' contrastnames{c}]); end
        
        contrasts(c,evinds1) = length(evinds2);
        contrasts(c,evinds2) = -length(evinds1);
        
    end
end

%% complicated contrasts

[~,~,cc] = read_contrasts(exp,us,runtype,model);
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

%%
[~,~,~,~,ftests] = read_contrasts(exp,us,runtype,model);
if ~isempty(ftests)
    ftmat = zeros(size(ftests,1),length(contrastnames));
    for i = 1:size(ftests,1)
        for j = 1:length(ftests{i,2})
            ftmat(i,strcmp(ftests{i,2}{j},contrastnames)) = 1;
        end
    end
end

%% open fsf file
fsffile = [fsfsubdir  'design_' runtype '_r' num2str(r) '_' model '_' params('smooth') 'mm_' stage '.fsf'];
fidfsf = fopen(fsffile,'w');
if fidfsf == -1
    error('Error in fla.m: unable to open %s\n',fsffile);
end

%% generate fsf file in parts, part 1, initial setup

fidtemp = fopen([fsftempdir 'part1.txt'],'r');
if fidtemp == -1
    error('Error in fla.m: unable to open %s\n',fidtemp);
end

if strcmp(tempderiv,'1')
    nevsreal = nevs*2;
elseif strcmp(tempderiv,'0')
    nevsreal = nevs;
end

if strcmp(stage, 'glm')
    inputvolume = funcvolume;
    stagenum = '3';
    filterfunc = '1';
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
end

clear keyname; clear value;
keyname{1}  = 'STAGENUM';         value{1}  = stagenum;
keyname{2}  = 'WATCHFEAT';        value{2}  = watchfeat;
keyname{3}  = 'OUTDIRECTORY';     value{3}  = featdir;
keyname{4}  = 'TRSECS';           value{4}  = trsecs;
keyname{5}  = 'NVOLUMES';         value{5}  = nvolumes;
keyname{6}  = 'FILTERFUNC';       value{6}  = filterfunc;
keyname{7}  = 'SMOOTHWIDTH';      value{7}  = params('smooth');
keyname{8}  = 'MAINSTATS';        value{8}  = mainstats;
keyname{9}  = 'MOTIONPARAMS';     value{9}  = motionparams;
keyname{10} = 'NEVSORIG';         value{10} = num2str(nevs);
keyname{11} = 'NEVSREAL';         value{11} = num2str(nevsreal);
keyname{12} = 'NCONORIG';         value{12} = num2str(length(contrastnames));
keyname{13} = 'NCONREAL';         value{13} = num2str(length(contrastnames));
keyname{14} = 'NFTESTSREAL';      value{14} = num2str(length(ftests));
keyname{15} = 'NFTESTSORIG';      value{15} = num2str(length(ftests));
keyname{16} = 'POSTSTATS';        value{16} = poststats;
keyname{17} = 'REGIMAGES';        value{17} = regimages;
keyname{18} = 'DOFANATOMICAL';    value{18} = '6';
keyname{19} = 'DOFSTANDARD';      value{19} = '12';
keyname{20} = 'STANDARDBRAIN';    value{20} = '/software/fsl/fsl-4.1.3/data/standard/MNI152_T1_2mm_brain';
keyname{21} = 'INPUTVOLUME';      value{21} = inputvolume;
keyname{22} = 'ANATOMICALVOLUME'; value{22} = anatomicalvolume;
keyname{23} = 'REGNONLINEAR';     value{23} = regnonlinear;
keyname{24} = 'SLICETIMING';      value{24} = slicetiming;
keyname{25} = 'HPCUTOFF';         value{25} = hpcutoff;

if strcmp(exp, 'sepi') && us == 1 && any(strcmp(runtype, {'continuous', 'continuous_sepi'})) && any(r == [1 2])
    keyname{26} = 'CONFOUNDREGRESSOR'; value{26} = '1';
else
    keyname{26} = 'CONFOUNDREGRESSOR'; value{26} = '0';
end

if strcmp(exp, 'sepi') && us == 1 && any(strcmp(runtype, {'continuous', 'continuous_sepi'})) && any(r == [1 2])
    fprintf(fidfsf, '\n\n');
    fprintf(fidfsf, '# Confound EVs text file for analysis 1\n');
    fprintf(fidfsf, 'set confoundev_files(1) "/mindhive/nklab/u/svnh/sepi/analysis/fla/usub1/tp1-2.txt"\n\n');
end

reptextMANYFAST(fidtemp, fidfsf, keyname, value);
fclose(fidtemp);

%% part 2,3, glm

fidtemp1 = fopen([fsftempdir 'part2.txt'], 'r');
if fidtemp1 == -1
    error('Error in fla.m: unable to open %s\n',fidtemp1);
end

fidtemp2 = fopen([fsftempdir 'part3.txt'], 'r');
if fidtemp2 == -1
    error('Error in fla.m: unable to open %s\n',fidtemp2);
end

clear keyname; clear value;
for e1 = 1:nevs
    keyname{1}{1} = 'EVINDEX';   value{1}{1} = num2str(e1); %#ok<*AGROW>
    keyname{1}{2} = 'EVNAME';    value{1}{2} = evname{e1};
    keyname{1}{3} = 'EVFILE';    value{1}{3} = evfile{e1};
    keyname{1}{4} = 'TEMPDERIV'; value{1}{4} = tempderiv;
    
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
if fidtemp == -1
    error('Error in fla.m: unable to open %s\n',fidtemp);
end

reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);

%% part 5,6, real contrasts
fidtemp1 = fopen([fsftempdir 'part5.txt'], 'r');
if fidtemp1 == -1
    error('Error in fla.m: unable to open %s\n',fidtemp1);
end

fidtemp2 = fopen([fsftempdir 'part6.txt'], 'r');
if fidtemp2 == -1
    error('Error in fla.m: unable to open %s\n',fidtemp2);
end

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
    
    % # F-test 1 element 1
    % set fmri(ftest_real1.1) 0
    for f = 1:size(ftests,1)
        fprintf(fidfsf, ['# F-test ' num2str(f) ' element ' num2str(c) '\n']);
        fprintf(fidfsf, ['set fmri(ftest_real' num2str(f) '.' num2str(c) ') ' num2str(ftmat(f,c)) '\n\n']);
    end
end

fclose(fidtemp1);
fclose(fidtemp2);

%% part 7,8, original contrasts

fidtemp1 = fopen([fsftempdir 'part7.txt'], 'r');
if fidtemp1 == -1
    error('Error in fla.m: unable to open %s\n',fidtemp1);
end

fidtemp2 = fopen([fsftempdir 'part8.txt'], 'r');
if fidtemp2 == -1
    error('Error in fla.m: unable to open %s\n',fidtemp2);
end

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
if fidtemp == -1
    error('Error in fla.m: unable to open %s\n',fidtemp);
end

reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);

%% part10, f-test

fidtemp = fopen([fsftempdir 'part10.txt'], 'r');
if fidtemp == -1
    error('Error in fla.m: unable to open %s\n',fidtemp);
end

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
if fidtemp == -1
    error('Error in fla.m: unable to open %s\n',fidtemp);
end

reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);
fclose(fidfsf);

%% call feat!

if ~optInputs(varargin,'norun')
    
    if overwrite && strcmp(stage,'glm') && exist([featdir '.feat'], 'dir');
        fprintf('Deleting file: %s\n',[featdir '.feat']);
        unix(['rm -r ' featdir '.feat']);
    end
    
    keyboard;
    fprintf(['feat ' fsffile '\n']);
    unix_fsl(fsl_version, ['feat ' fsffile], varargin{:});
    
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




