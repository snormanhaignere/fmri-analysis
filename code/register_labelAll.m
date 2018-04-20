function register_labelAll(exp,usubs,varargin)

% select labels
labels = mydir('~/freesurfer/fsaverage/label/');
labels = strrep(labels,'rh.','');
labels = strrep(labels,'lh.','');
inds = false(size(labels));
for i = 1:length(labels);
    if findstr(labels{i},'.label')
        inds(i) = true;
    end
end
labels = labels(inds);
labels = strrep(labels,'.label','');
labels = unique(labels);

% specify labels
if optInputs(varargin, 'labels')
    labels = varargin{optInputs(varargin, 'labels')+1};
end
fsl_version = read_fsl_version(exp,varargin{:});
analysisdir = [params('rootdir') exp '/analysis/'];
model = 'block';

for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp,usubs(i),'func');
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(labels)
        
        fprintf('%s, highres: usub %d, %s\n',exp,usubs(i),labels{j}); drawnow;
        
        if optInputs(varargin, 'rh')
            register_label2highres(exp,usubs(i),labels{j},'rh',varargin{:});
        elseif optInputs(varargin, 'lh')
            register_label2highres(exp,usubs(i),labels{j},'lh',varargin{:});
        else
            
            register_label2highres(exp,usubs(i),labels{j},'rh',varargin{:});
            register_label2highres(exp,usubs(i),labels{j},'lh',varargin{:});
            
            
            rightlabel = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_rh-' labels{j}  '_highres_highres2standard_freesurfer.nii.gz'];
            leftlabel = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_lh-' labels{j} '_highres_highres2standard_freesurfer.nii.gz'];
            bilabel = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_' labels{j} '_highres_highres2standard_freesurfer.nii.gz'];
            oldbilabel = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_' labels{j} '.nii.gz'];
            
            if ~exist(bilabel,'file') && exist(oldbilabel,'file')
                copyfile(oldbilabel,bilabel);
            end
            
            if ~exist(bilabel,'file') || optInputs(varargin, 'overwrite');
                unix_fsl(fsl_version, ['fslmaths ' rightlabel ' -add ' leftlabel ' ' bilabel]);
            end
            
            rightlabel_subsamp = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_rh-' labels{j}  '_highres_highres2standard_freesurfer_2mm.nii.gz'];
            leftlabel_subsamp = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_lh-' labels{j} '_highres_highres2standard_freesurfer_2mm.nii.gz'];
            bilabel_subsamp = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_' labels{j} '_highres_highres2standard_freesurfer_2mm.nii.gz'];
            
            if ~exist(rightlabel_subsamp,'file') || optInputs(varargin, 'overwrite')
                unix_fsl(fsl_version, ['fslmaths ' rightlabel ' -subsamp2 ' rightlabel_subsamp]);
            end
            
            if ~exist(leftlabel_subsamp,'file') || optInputs(varargin, 'overwrite')
                unix_fsl(fsl_version, ['fslmaths ' leftlabel ' -subsamp2 ' leftlabel_subsamp]);
            end
            
            if ~exist(bilabel_subsamp,'file') || optInputs(varargin, 'overwrite')
                unix_fsl(fsl_version, ['fslmaths ' bilabel ' -subsamp2 ' bilabel_subsamp]);
            end
            %
            %             rightlabel_subsamp = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_rh-' labels{j}  '_highres_highres2standard_freesurfer_2mm.nii.gz'];
            %             leftlabel_subsamp = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_lh-' labels{j} '_highres_highres2standard_freesurfer_2mm.nii.gz'];
            %             bilabel_subsamp = [analysisdir 'sla/usub' num2str(usubs(i)) '/masks/mylabel_' labels{j} '_highres_highres2standard_freesurfer_2mm.nii.gz'];
            %
            
        end
        
        
        for m = 1:length(runtypes)
            allruns = read_runs(exp,usubs(i),runtypes{m});
            for k = 1:length(allruns)
                fprintf('func: usub %d, %s, %s, run %d\n', usubs(i),runtypes{m},labels{j},allruns(k)); drawnow;
                if optInputs(varargin, 'rh')
                    register_label2func(exp,['mylabel_rh-' labels{j}],usubs(i),runtypes{m},allruns(k),model,'overwrite',varargin{:});
                elseif optInputs(varargin, 'lh')
                    register_label2func(exp,['mylabel_lh-' labels{j}],usubs(i),runtypes{m},allruns(k),model,'overwrite',varargin{:});
                else
                    register_label2func(exp,['mylabel_rh-' labels{j}],usubs(i),runtypes{m},allruns(k),model,'overwrite',varargin{:});
                    register_label2func(exp,['mylabel_lh-' labels{j}],usubs(i),runtypes{m},allruns(k),model,'overwrite',varargin{:});
                    register_label2func(exp,['mylabel_' labels{j}],usubs(i),runtypes{m},allruns(k),model,'overwrite',varargin{:});
                end
            end
        end
    end
end