function sla_surfAll(exp,usubs,fwhm,varargin)

% function sla_surfAll(usubs,model,varargin)
%
% wrapper function for sla_surf.m

model = 'block';
if ~optInputs(varargin, 'subject_specific_standard')
    varargin = [varargin, {'fsaverage'}];
end

correction = 'cache';
if optInputs(varargin,'correction')
    correction = varargin{optInputs(varargin,'correction')+1};
end

con_sign = 'pos';
if optInputs(varargin,'con_sign')
    con_sign = varargin{optInputs(varargin,'con_sign')+1};
end

voxthresh = 3;
if optInputs(varargin,'voxthresh')
    voxthresh = varargin{optInputs(varargin,'voxthresh')+1};
end

nsim = 10000;
if optInputs(varargin,'nsim')
    nsim = varargin{optInputs(varargin,'nsim')+1};
end

hemis = {'rh','lh'};
if optInputs(varargin, 'rh');
    hemis = {'rh'};
end
if optInputs(varargin, 'lh');
    hemis = {'lh'};
end

for i = 1:length(usubs);
    
    runtypes = read_runtypes(exp,usubs(i),'func');
    if optInputs(varargin,'runtypes');
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for m = 1:length(runtypes)
        
        contrasts = read_contrasts(exp,usubs(i),runtypes{m},model);
        if optInputs(varargin,'contrasts');
            contrasts = varargin{optInputs(varargin,'contrasts')+1};
        end
        
        for j = 1:length(contrasts)
            
            allruns = read_runs(exp,usubs(i),runtypes{m},'contrast',contrasts{j},model);
            runset = {allruns};
            
            % loop through sets of runs
            for k = 1:length(runset)
                
                if isempty(runset{k});
                    warning(['Skipping at least 1 run set of ' contrasts{j}]); %#ok<WNTAG>
                    continue;
                end
                
                for q = 1:length(hemis)
                    fprintf('sub %d, %s, runs %s, %s, %s\n',usubs(i),runtypes{m},sprintf('%d',runset{k}),contrasts{j},hemis{q});
                    if optInputs(varargin, 'monkey')
                        sla_surf_monkey(exp,usubs(i),contrasts{j},fwhm,correction,con_sign,voxthresh,nsim,hemis{q},runtypes{m},model,varargin{:});
                    else
                        sla_surf(exp,usubs(i),contrasts{j},fwhm,correction,con_sign,voxthresh,nsim,hemis{q},runtypes{m},model,varargin{:});
                    end
                end
            end
        end
    end
end
