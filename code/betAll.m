function betAll(exp,usubs,varargin)

% function bet(usubs,varargin)
%
% bet skull-strips the reoriented data using
% fsl's skull-stripping utility
% fa values must be specified for each new subject

datadir = [params('rootdir') exp '/data/'];

if optInputs(varargin, 'usubs');
    usubs = varargin{optInputs(varargin,'usubs')+1};
end

for i = 1:length(usubs)
    
    % fa values
    switch exp
        otherwise
            fafunc = 0.08;
            fastruct = 0.15;
    end
    
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
            lasfile = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS'];
            
            if strcmp(runtypes{j},'struct');
                favalue = fastruct;
            else
                favalue = fafunc;
            end
            
            if ~exist([lasfile '_brain.nii.gz'],'file') || optInputs(varargin, 'overwrite');
                
                fprintf(['bet ' lasfile ' ' lasfile '_brain' ' -R  -f ' num2str(favalue) ' -g 0 -n -m' '\n']);
                unix(['bet ' lasfile ' ' lasfile '_brain' ' -R  -f ' num2str(favalue) ' -g 0 -n -m']);
                
                fprintf(['fslmaths ' lasfile ' -mas ' lasfile '_brain_mask' ' ' lasfile '_brain' '\n']);
                unix(['fslmaths ' lasfile ' -mas ' lasfile '_brain_mask' ' ' lasfile '_brain']);
                
            end
        end
    end
end

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