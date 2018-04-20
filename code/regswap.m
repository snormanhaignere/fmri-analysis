function regswap(exp,us,runtype,r,model,regtype,varargin)

analysisdir = [params('rootdir') exp '/analysis/'];
fwhm = read_smooth(exp,varargin{:});
featdir = [analysisdir 'fla/usub' num2str(us) '/' runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/'];

fsl_version = read_fsl_version(exp);

switch fsl_version
    case '5.0'
        % files to swap based on the new example_func2highres registration
        swapfiles = {'example_func2highres.mat','example_func2highres.nii.gz','example_func2highres.png',...
            'example_func2standard1.png','example_func2standard.mat','example_func2standard.nii.gz','example_func2standard.png',...
            'highres2example_func.mat','standard2example_func.mat'};
    otherwise
        % files to swap based on the new example_func2highres registration
        swapfiles = {'example_func2highres1.png','example_func2highres2.png','example_func2highres.mat','example_func2highres.nii.gz','example_func2highres.png',...
            'example_func2standard1.png','example_func2standard2.png','example_func2standard.mat','example_func2standard.nii.gz','example_func2standard.png',...
            'highres2example_func.mat','standard2example_func.mat'};
end

%% Moving in new bbregs

copyflirt(exp,us,runtype,r,model); % always make sure initial flirts have been copied

flirtdir = [featdir 'reg/flirt/'];
regdir = [featdir 'reg/' regtype '/']; % might be the same as flirtdir if regtype is flirt

for i = 1:length(swapfiles)
    
    %     if ~exist([flirtdir swapfiles{i}],'file');
    %         error(['Error: could not find file ' flirtdir swapfiles{i}]);
    %     end;
    
    fprintf('Swapping in new %s',swapfiles{i});
    
    src = [regdir swapfiles{i}];
    dest = [featdir 'reg/' swapfiles{i}];
    unix(['cp -f ' src ' ' dest]);

end


%% scraps

% 
% 
% if strcmp(regtype,'flirt') || ~exist(flirtdir,'dir');
%     
%     mkdir(flirtdir);
%     for i = 1:length(swapfiles);
%         
%         fprintf('Moving old %s',swapfiles{i});
%         
%         src = [featdir 'reg/' swapfiles{i}];
%         dest = [flirtdir swapfiles{i}];
%         unix(['mv ' src ' ' dest]);
%         
%     end
% end