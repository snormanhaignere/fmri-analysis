function reorientBXH(exp,usubs,varargin)

% function reorientBXH(usubs,varargin)
%
% reorients data converted by convertraw.m
% using the bxhreorient command written by syam gadde

datadir = [params('rootdir') exp '/data/'];
fsl_version = read_fsl_version(exp,varargin{:});

for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp,usubs(i),'raw',varargin{:});
    if optInputs(varargin, 'runtypes');
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(runtypes)
        
        runnum = read_runs(exp,usubs(i),runtypes{j},varargin{:});
        if optInputs(varargin,'runnum')
            runnum = varargin{optInputs(varargin,'runnum')+1};
        end
        
        for k = 1:length(runnum)
            
            niidir = [datadir 'brain/nifti/usub' num2str(usubs(i)) '/'];
            niifile = [niidir runtypes{j} '_r' num2str(runnum(k))];
            lasfile = [niifile '_LAS'];
            
            if ~exist([lasfile '.nii.gz'],'file') || optInputs(varargin, 'overwrite');
                
                
                if exist([lasfile '.nii.gz'],'file'); fprintf('Deleting file: %s\n', [lasfile '.nii.gz']); delete([lasfile '.nii.gz']); end
                if exist([lasfile '.bxh'],'file'); fprintf('Deleting file: %s\n', [lasfile '.bxh']); delete([lasfile '.bxh']); end
                
                %                 unix('rm temp.*');
                id = [exp, num2str(usubs(i)), runtypes{j}, num2str(runnum(k))];
                fprintf(['fslwrapbxh ' niidir '\n']);
                unix(['fslwrapbxh ' niidir]);
                fprintf(['bxhreorient --orientation=LAS ' niifile '.bxh ' niidir 'temp' id '.bxh \n']);
                unix(['bxhreorient --orientation=LAS ' niifile '.bxh ' niidir 'temp' id '.bxh']);
                fprintf(['bxh2analyze -s --niigz ' niidir 'temp' id '.bxh ' lasfile '\n']);
                unix(['bxh2analyze -s --niigz ' niidir 'temp' id '.bxh ' lasfile]);
                fprintf(['rm ' niidir 'temp' id '.* \n']);
                unix(['rm ' niidir 'temp'  id '.*']);
            end
            
            % preprocessing directory
            preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(usubs(i)) '/' runtypes{j} '_r' num2str(runnum(k)) '/'];
            if ~exist(preprocdir,'dir');
                mkdir(preprocdir);
            end
            
            % input file
            preprocessfile = [preprocdir 'raw.nii.gz'];
            if ~exist(preprocessfile, 'file') || optInputs(varargin, 'overwrite')
                copyfile([lasfile '.nii.gz'], preprocessfile, 'f');
            end
            
        end
    end
end

