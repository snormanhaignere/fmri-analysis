function combine_runs(exp, us, runtype, rungroups, varargin)

fsl_version = read_fsl_version(exp,varargin{:});
for i = 1:length(rungroups)
    
    preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_combined_r' num2str(i) '/'];
    if ~exist(preprocdir,'dir');
        mkdir(preprocdir)
    end
    
    % combine functional data
    outputfile = [preprocdir 'motcorr.nii.gz'];
    if ~exist(outputfile, 'file') || optInputs(varargin, 'overwrite')
        cat_data = [];
        for j = 1:length(rungroups{i})
            
            fprintf('Run %d\n',rungroups{i}(j)); drawnow;

            mcflirt(exp,us,runtype,rungroups{i}(j),varargin{:}); % motion correct
            if j == 1;
                % copy over example_func file from first runtype
                if ~exist([preprocdir 'example_func.nii.gz'],'file') || optInputs(varargin, 'overwrite')
                    copyfile([params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(j)) '/example_func.nii.gz'], [preprocdir 'example_func.nii.gz'],'f');
                end
                
                % compute voxel means and store data matrix
                br = readmr([params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(j)) '/motcorr.nii.gz'],'NOPROGRESSBAR');
                dims = size(br.data);
                voxel_means = nan([dims(1:3), length(rungroups{i})]);
                voxel_means(:,:,:,j) = mean(br.data,4);
                cat_data = br.data - repmat(voxel_means(:,:,:,j), [1 1 1 dims(4)]);
                
            else
                % register motion-corrected volume to first runtype
                exfunc_to_move = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(j)) '/example_func.nii.gz'];
                exfunc_target = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(1)) '/example_func.nii.gz'];
                regmat = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(j)) '/r' num2str(rungroups{i}(j)) '_to_r' num2str(rungroups{i}(1)) '.mat'];
                
                motcorr_to_move = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(j)) '/motcorr.nii.gz'];
                motcorr_target = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(1)) '/motcorr.nii.gz'];
                motcorr_final = [params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(rungroups{i}(j)) '/motcorr_reg_to_r' num2str(rungroups{i}(1)) '.nii.gz'];
                
                if ~exist(motcorr_final, 'file') || optInputs(varargin, 'overwrite')
                    unix_fsl(fsl_version, ['flirt -in ' exfunc_to_move ' -ref ' exfunc_target ' -omat ' regmat]);
                    unix_fsl(fsl_version, ['flirt -interp trilinear -in ' motcorr_to_move ' -ref ' motcorr_target ' -applyxfm -init ' regmat ' -out ' motcorr_final]);
                end
                
                % compute voxel means and store data matrix
                br = readmr(motcorr_final,'NOPROGRESSBAR');
                dims = size(br.data);
                voxel_means(:,:,:,j) = mean(br.data,4);
                cat_data = cat(4, cat_data, br.data - repmat(voxel_means(:,:,:,j), [1 1 1 dims(4)]));
            end
            
        end
        
        % store data with voxel mean added back in
        br.data = cat_data + repmat(mean(voxel_means,4), [1 1 1 size(cat_data,4)]);
        br.info.dimensions(4).size = size(cat_data,4);
        bxhwrite(br,outputfile);
    end
    
    % combine timing files
    datadir = [params('rootdir') exp '/data/experiment/'];
    new_datafile = [datadir 'usub' num2str(us) '/seq/' runtype '_combined_r' num2str(i) '_block.par'];
    %     if ~exist(new_datafile,'file') || optInputs(varargin,'overwrite')
    fid = fopen(new_datafile, 'w');
    onset_index = 0;
    for j = 1:length(rungroups{i})
        b = read_timings(exp,us,runtype,rungroups{i}(j),varargin{:});
        for k = 1:length(b.onsets)
            if length(b.conds{k}) > 48
                error('Condition string too long');
            end
            fprintf(fid, '%8.2f%5d%8.2f%5d%50s\n', b.onsets(k) + onset_index, b.condition_indices(k), b.durs(k), 1, b.conds{k});
        end
        [~, ~, TR, ~, ~, ~, ~, nTR, ~] = read_scanparams(exp,us,runtype,'run',rungroups{i}(j),varargin(:));
        onset_index = onset_index + nTR*TR;
    end
    fclose(fid);
    
    %     end
    
    %     clear br;
    %     cat_data = [];
    %     for j = 1:length(old_runtypes)
    %         br = readmr([params('rootdir') exp '/analysis/preprocess/usub' num2str(us) '/' old_runtypes{j} '_r' num2str(runs(i)) '/smooth' params('smooth') 'mm_intnorm_hpfilt' num2str(params('hpcutoff')) '.nii.gz'],'NOPROGRESSBAR');
    %         cat_data = cat(4, cat_data, br.data);
    %     end
    %
    %     br_new = br;
    %     [blockdur, nulldur, TR, TA, stimdur, stim2scan, win, nTR, disdaqs] = read_scanparams(exp,us,new_runtype,varargin(:));
    %     br_new.info.dimensions(4).size = nTR;
    %     br_new.data = cat_data;
    
end