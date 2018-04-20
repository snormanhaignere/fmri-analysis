function register_handpick_roiAll(rois,exp,usubs,varargin)

model = 'block';

for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp,usubs(i),'func');
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(rois)
        for m = 1:length(runtypes)
            allruns = read_runs(exp,usubs(i),runtypes{m});
            for k = 1:length(allruns)
                fprintf('func: usub %d, %s, %s, run %d\n', usubs(i),runtypes{m},rois{j},allruns(k)); drawnow;
                register_handpick_roi(exp,rois{j},usubs(i),runtypes{m},allruns(k),model,varargin{:});
            end
        end
    end
end