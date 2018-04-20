function sla(exp,us,runtype,runnum,contr,model,varargin)

% function sla(s,runtype,runnum,contr,model,varargin)
% 
% main function for running second level analyses
% 
% EDIT: in future versions, should add runtype to gfeat and design

analysisdir = [params('rootdir') exp '/analysis/'];
flasubdir = [analysisdir 'fla/usub' num2str(us) '/'];
fsftempdir = [params('rootdir') 'fmri-analysis/code/fsf-templates/sla/'];

fwhm = read_smooth(exp, varargin{:});
fsfsubdir = [analysisdir 'fsf/sla/usub' num2str(us) '/'];
if ~exist(fsfsubdir, 'dir'); 
    mkdir(fsfsubdir); 
end

slasubdir = [analysisdir 'sla/usub' num2str(us) '/'];
if ~exist(slasubdir, 'dir'); 
    mkdir(slasubdir); 
end

gfeat = [slasubdir runtype '_' contr '_r' sprintf('%d',runnum) '_' model '_' num2str(fwhm*100,'%.0f') 'mm'];
if length(runnum) > 100
    gfeat = strrep(gfeat, ['_r' sprintf('%d',runnum)], ['_r' DataHash(runnum)]);
end
fsl_version = read_fsl_version(exp);

if exist([gfeat '.gfeat/cope1.feat/stats/zstat1.nii.gz'],'file') && ~optInputs(varargin, 'overwrite');
    return;
end

%% parameters

voxelthresh = '3.09';
clusterthresh = '0.05';

copefiles = cell(1,length(runnum));
for r = 1:length(runnum)
    featdir = [flasubdir runtype '_r' num2str(runnum(r)) '_' model  '_' num2str(fwhm*100,'%.0f') 'mm.feat/'];
    fid = fopen([featdir 'contrastnames.txt'],'r');
    tmp = textscan(fid,'%s\n'); fclose(fid); contrastnames = tmp{1};
    ind = strmatch(contr,contrastnames,'exact');
    if isempty(ind); error('Error: no matching contrast'); end
    if length(ind) > 1
        ind = ind(1);
        warning('Multiple contrast matches');
    end
    copefiles{r} = [featdir 'stats/cope' num2str(ind) '.nii.gz'];
end

fsffile = [fsfsubdir 'design_' runtype '_' contr '_r' sprintf('%d',runnum) '_' model  '_' num2str(fwhm*100,'%.0f') 'mm.fsf'];
if length(runnum) > 100
    fsffile = strrep(fsffile, ['_r' sprintf('%d',runnum)], ['_r' DataHash(runnum)]);
end
fidfsf = fopen(fsffile,'w');

%%

templatefile = [fsftempdir 'part1.txt'];
fidtemp = fopen(templatefile,'r');

clear keyname; clear value;
keyname{1}  = 'WATCHFEAT';        value{1}  = '0';
keyname{2}  = 'OUTDIRECTORY';     value{2}  = gfeat;
keyname{3}  = 'NCOPES';           value{3}  = num2str(length(copefiles));
keyname{4}  = 'CLUSTERTHRESH';    value{4}  = clusterthresh;
keyname{5}  = 'VOXELTHRESH';      value{5}  = voxelthresh;
if optInputs(varargin, 'subject_specific_standard')
    keyname{6}  = 'STANDARDBRAIN';    value{6} = [analysisdir 'preprocess/usub' num2str(us) '/struct_r1/brain'];
elseif optInputs(varargin, 'functional_reference_volume')
    allruns = read_runs(exp, us, runtype, varargin{:});
    keyname{6}  = 'STANDARDBRAIN';    value{6} = [analysisdir 'preprocess/usub' num2str(us) '/' runtype '_r' num2str(allruns(1)) '/example_func.nii.gz'];
else
    keyname{6}  = 'STANDARDBRAIN';    value{6} = '/software/fsl/fsl-4.1.3/data/standard/MNI152_T1_2mm_brain';
end

reptextMANYFAST(fidtemp, fidfsf, keyname, value);
fclose(fidtemp);

%%

templatefile = [fsftempdir 'part2.txt'];
fidtemp = fopen(templatefile,'r');

clear keyname; clear value;
for r = 1:length(runnum)
    keyname{1} = 'COPEINDEX';  value{1} = num2str(r); %#ok<*AGROW>
    keyname{2} = 'COPEFILE';   value{2} = copefiles{r};
    reptextMANYFAST(fidtemp, fidfsf, keyname, value);
    fseek(fidtemp,0,-1);
end

fclose(fidtemp);

%%

templatefile = [fsftempdir 'part3.txt'];
fidtemp = fopen(templatefile,'r');
reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);

%%

templatefile = [fsftempdir 'part4.txt'];
fidtemp = fopen(templatefile,'r');

clear keyname; clear value;
for r = 1:length(runnum)
    keyname{1} = 'COPEINDEX';  value{1} = num2str(r);
    reptextMANYFAST(fidtemp, fidfsf, keyname, value);
    fseek(fidtemp,0,-1);
end

fclose(fidtemp);
 
%%

templatefile = [fsftempdir 'part5.txt'];
fidtemp = fopen(templatefile,'r');

clear keyname; clear value;
for r = 1:length(runnum)
    keyname{1} = 'COPEINDEX';  value{1} = num2str(r);
    reptextMANYFAST(fidtemp, fidfsf, keyname, value);
    fseek(fidtemp,0,-1);
end

fclose(fidtemp);

%%

templatefile = [fsftempdir 'part6.txt'];
fidtemp = fopen(templatefile,'r');
reptextMANYFAST(fidtemp, fidfsf, {}, {});
fclose(fidtemp);
fclose(fidfsf);

%%
if exist([gfeat '.gfeat'], 'dir'); 
    fprintf('Deleting file: %s\n',[gfeat '.gfeat']);
    unix(['rm -r ' gfeat '.gfeat']); 
end

unix_fsl(fsl_version,['feat ' fsffile], varargin{:});

%%
filt_func_file = [gfeat '.gfeat/cope1.feat/filtered_func_data.nii.gz'];
if exist(filt_func_file,'file')
    delete(filt_func_file);
end

var_filt_func_file = [gfeat '.gfeat/cope1.feat/var_filtered_func_data.nii.gz'];
if exist(var_filt_func_file,'file')
    delete(var_filt_func_file);
end

res4d_file = [gfeat '.gfeat/cope1.feat/stats/res4d.nii.gz'];
if exist(res4d_file,'file')
    delete(res4d_file);
end

