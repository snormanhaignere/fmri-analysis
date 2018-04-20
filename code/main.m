function main(s)

%% to add new experiment
% read_runtypes
% read_conditions
% read_scanparams
% read_plotorder
% read_timings
% preprocAll
% flaAll
% read_fsl_version
% read_freesurfer_version
% read_pipeline_version
% read_contrasts
% read_smooth
% evgroup2evname
% check read_runs
% run convert2usubsAll for old studies
% read_reference_frame
% read_bet_fa
% recon directory

%%
% preprocAll
% read_fsl_version
% read_freesurfer_version
% flaAll

%% updated
% convertraw
% bet
% reorientBXH
% read_runtypes
% freesurfer
% stfAll
% read_models
% stf
% read_conditions
% read_scanparams

%% 

% to add subject
% fla/sub -> fla/usub
% sla/sub -> sla/usub
% runorders
% read runtypes, contrasts, conditions

%% New pipeline

preprocAll(exp,s);
flaAll(exp,s,'block','glm');
bbregisterAll(exp,s,'block');
flaAll(exp,s,'block','contrasts');
slaAll(s,'block');

%%

convertraw(exp,s);
reorientBXH(exp,s);
betAll(exp,s);
stfAll(exp,s);
flaAll(exp,s,'block','glm');
bbregisterAll(exp,s,'block');
flaAll(exp,s,'block','contrasts');
slaAll(s,'block');

%% subject 1

% files to update when starting new experiment
% bet: new fa values
% read_runtypes
% read_conditions
% read_models
% read_scanparams
% read_timings: format string used when reading in files
% stf: for how to handle different model types
% params: choose smoothing kernel
% read_plotorder
% read_plotformat

% read_contrasts
% fla: number of volumes, model
% copy fsf template files
% evgroup2evname: group to ev mapping

% preprocessing
convertraw(s);
reorientBXH(s);
bet(s);

% freesurfer reconstruction
freesurfer(s,'nohup')

% timing files
stfAll(s,'main');

% first level regression
flaAll(s,'block','glm');
% flaAll(s,'main','event','glm');

% Freesurfer needs to complete before next stage

% create and replace with bbregister transforms
bbregisterAll(s,'block');
% bbregisterAll(s,'event');

% first-level contrasts
flaAll(s,'block','contrasts');
% flaAll(s,'event','contrasts');

% second-level contrasts
slaAll(s,'block');
% slaAll(s,'event');

% highres versions of contrasts
flirt_sla2highresAll(s,'block');
flirt_sla2highresAll(s,'event');

% register anatomicals
format_atlasAll; % only needs to be run once at beginning of project
register_atlasAll(s);
register_anat_freesurferAll(s);
combine_atlasAll(s);

% up to hear in update

% roi analyses
roi_plotcondsAll([3 4 6 7 8 10 11],'main','anat_contrast','destrieux_pp-ics','meanperiodicity','event','weightorder',1);
roi_itemAll([3 4 6 7 8 10 11],'main','anat_contrast','destrieux_pp-ics','meanperiodicity','event','weightorder',1);
roianalyses; % all roi analyses with anatomical constraints

% auxilary stuff

% cluster the contrasts
cluster_highresAll(s,'block');
cluster_highresAll(s,'event');

% brains with all stims or conditions combined
allstimsbr(s);
allcondsbr(s,'zstat');
allcondsbr(s,'pe');

% clustering analysis
cluster_func(s,'pe');

% hand-pick rois, *MANUAL*
load_contrasts(s,'orchestra_vs_nonwords','block'); % sample
pick_clusters(s,'orchestra_vs_nonwords','block','atg'); % sample

% run roi analysis for hand-picked clusters
roi_categoryAll([3 4 6 7 8 10],'cluster','r-sts','meanperiodicity','event','weightorder',1);
roi_itemAll([3 4 6 7 8 10],'cluster','r-sts','meanperiodicity','event','weightorder',1);

tla(s,'orchestra_vs_nonwords','block');
tlaAll([3 4 6 7 8 10]);

register_greymatter_highresAll(s);
register_greymatter_mni;
tsnr_runsAll(3:9);
tsnr_submeanAll(3:9);

% bbregister_anat2featAll(s,'block');
% bbregister_anat2featAll(s,'event');
% 
% bbregister_flirtfilesAll(s,'block');
% bbregister_flirtfilesAll(s,'event');

% regswapAll(s,'block','bbreg');
% regswapAll(s,'event','bbreg');
% 
% featregapplyAll(s,'block');
% featregapplyAll(s,'event');

%%

svgdir = '~/pitch_dp/figures/';
roinames = {'mylabel_rh-te11','mylabel_lh-te11','mylabel_rh-te10','mylabel_lh-te10','mylabel_rh-te12','mylabel_lh-te12','mylabel_rh-pp','mylabel_lh-pp','mylabel_rh-pt','mylabel_lh-pt'};

for i = 1:length(roinames)
    roi_pscAll(roinames{i},'thr',0.5,'','',[],'sigav_plateau',5,'withsub_sem','simplot');
    ylim([0 1.5]);
    plot2svg([svgdir roinames{i} '.svg'],gcf,'svg');
end



%%

svgdir = '~/pitch_dp/figures/';
subjnum = [15, 5, 11, 8, 10,  7, 2, 6];

for i = 1:length(subjnum)
    roi_pscAll(subjnum(i),'mylabel_pc_fwhm8_p2','thr',0.5,'harm_vs_noise','nvox',[1 20],'sigav_plateau',5,'simplot');
    plot2svg([svgdir 's' num2str(subjnum(i)) '_pc_fwhm_p2.svg'],gcf,'svg');
end

%% Correlation between tonotopy and pitch selectivity

subjnum = setdiff(1:15,[4 12 14]);
nsub = length(subjnum);

contrasts = {'harm-lowfreq_vs_noise-lowfreq','harm-highfreq_vs_noise-highfreq','harm-unres_vs_noise-highfreq',...
    'harm-lowfreq-lowmask_vs_noise-lowfreq-lowmask','harm-highfreq-lowmask_vs_noise-highfreq-lowmask','harm-unres-lowmask_vs_noise-highfreq-lowmask'};
ncontrasts = length(contrasts);

roi = 'mylabel_stp';

rval_all = nan(nsub,ncontrasts);
zval_all = nan(nsub,ncontrasts);
for i = 1:ncontrasts
    [rval zval] = roi_corr_contrastAll(subjnum,'mylabel_stp','thr',0.5,{'highfreq_vs_lowfreq',contrasts{i}},'noplot','contype',{'neg','pos'});
    rval_all(:,i) = rval(:,2,1);
    zval_all(:,i) = zval(:,2,1);
end

rval_mean = mean(rval_all);
rval_sem = std(rval_all)/sqrt(nsub-1);

mybar(rval_mean,contrasts,contrasts,'errorbar',rval_sem,'ylim',[0 1]);
title('Tonotopy vs. Pitch Correlations in the Superior Temporal Plane');

%% Confusion matrix between resolved and unresolved harmonics

subjnum = setdiff(1:15,[4 12 14]);
nsub = length(subjnum);

contrasts = {'harm-lowfreq_vs_noise-lowfreq','harm-highfreq_vs_noise-highfreq','harm-unres_vs_noise-highfreq',...
    'harm-lowfreq-lowmask_vs_noise-lowfreq-lowmask','harm-highfreq-lowmask_vs_noise-highfreq-lowmask','harm-unres-lowmask_vs_noise-highfreq-lowmask'};

[rval zval] = roi_corr_contrastAll(subjnum,'mylabel_stp','thr',0.5,contrasts,'noplot');

res_to_res = mean([rval(:,1,1) rval(:,2,2)],2);
res_to_unres = mean([rval(:,1,3) rval(:,2,3)],2);
unres_to_unres = rval(:,3,3);
res_to_res_mask = mean([rval(:,4,4) rval(:,5,5)],2);
res_to_unres_mask = mean([rval(:,4,6) rval(:,5,6)],2);
unres_to_unres_mask = rval(:,6,6);

res_unres_all = [res_to_res res_to_unres unres_to_unres res_to_res_mask res_to_unres_mask unres_to_unres_mask];
res_unres_mean = mean(res_unres_all);
res_unres_sem = std(res_unres_all)/sqrt(nsub-1);

barnames = {'Res to Res','Unres to Res','Unres to Unres','Res to Res (Mask)','Unres to Res (Mask)', 'Unres to Unres (Mask)'};
mybar(res_unres_mean,barnames,barnames,'errorbar',res_unres_sem,'ylim',[0 1]);
title('Pitch Correlations in the Superior Temporal Plane');

%% 

subjnum = setdiff(1:15,[4 12 14]);
nsub = length(subjnum);

conds = read_conditions(1,'main');
nconds = length(conds);

nbins = 10;
psc = nan(nsub,nconds,nbins);

for i = 1:nbins;
    [~,~,psc(:,:,i)] = roi_pscAll(subjnum,'mylabel_stp','thr',0.5,'highfreq_vs_lowfreq','perc',[0 0.1],'sigav_plateau',5,'withsub_sem','negcontrast');
end

mybar(res_unres_mean,barnames,barnames,'errorbar',res_unres_sem,'ylim',[0 1]);


%%

nbins = 10;
selectivity_all = nan(4,nbins);
lstderr = nan(4,nbins);
ustderr = nan(4,nbins);
lconfint = nan(4,nbins);
uconfint = nan(4,nbins);
barnames = cell(1,nbins);
for i = 1:nbins
    barnames{i} = ['Bin ' num2str(i)];
    [selectivity_all(:,i) lstderr(:,i) ustderr(:,i) lconfint(:,i) uconfint(:,i)] = roi_pitchselectivity_bootstrap({'mylabel_stp'},'all','pos','nvox',[1 20*nbins],'loc2','highfreq_vs_lowfreq','pos','nvox',[1 20] + 20*(i-1),'noplot');
end


for i = 1:4
mybar(selectivity_all(i,:),barnames,barnames,'errorbar',lstderr(i,:),ustderr(i,:));
ylim([0 0.5]);
end

%%
% %% Confusion matrix between resolved and unresolved harmonics
% 
% subjnum = setdiff(1:15,[4 12 14]);
% nsub = length(subjnum);
% 
% contrasts = {'harm_vs_noise','harm-unres_vs_noise-highfreq',...
%     'harm-lowmask_vs_noise-lowmask','harm-unres-lowmask_vs_noise-highfreq-lowmask'};
% 
% [rval zval] = roi_corr_contrastAll(subjnum,'mylabel_stp','thr',0.5,contrasts,'noplot');
% 
% res_to_res = rval(:,1,1);
% res_to_unres = rval(:,1,2);
% unres_to_unres = rval(:,2,2);
% res_to_res_mask = rval(:,3,3);
% res_to_unres_mask = rval(:,3,4);
% unres_to_unres_mask = rval(:,4,4);
% 
% res_unres_all = [res_to_res res_to_unres unres_to_unres res_to_res_mask res_to_unres_mask unres_to_unres_mask];
% 
% barnames = {'Res to Res','Unres to Res' 'Unres to Unres','Res to Res (Mask)','Unres to Res (Mask)' 'Unres to Unres (Mask)'};
% mybar(mean(res_unres_all),barnames,barnames,'errorbar',stderr_withsub(res_unres_all));


