function copyflirt(exp,us,runtype,r,model,varargin)

% function copyflirt(s,runtype,r,model)
% 
% saves a backup copy of the original flirt registration files
% e.g. ensures that bbregister doesn't overwrite these

analysisdir = [params('rootdir') exp '/analysis/'];

fwhm = read_smooth(exp, varargin{:});
featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(fwhm*100,'%.0f') 'mm.feat/'];

% files to swap based on the new example_func2highres registration
swapfiles = {'example_func2highres1.png','example_func2highres2.png','example_func2highres.mat','example_func2highres.nii.gz','example_func2highres.png',...
    'example_func2standard1.png','example_func2standard2.png','example_func2standard.mat','example_func2standard.nii.gz','example_func2standard.png',...
    'highres2example_func.mat','standard2example_func.mat'};

flirtdir = [featdir 'reg/flirt/'];

if ~exist(flirtdir,'dir');
    
    mkdir(flirtdir);
    
    for i = 1:length(swapfiles);
        
        fprintf('Moving old %s',swapfiles{i});
        
        src = [featdir 'reg/' swapfiles{i}];
        dest = [flirtdir swapfiles{i}];
        unix(['cp ' src ' ' dest]);
        
    end
    
end