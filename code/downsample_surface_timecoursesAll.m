function downsample_surface_timecoursesAll(exp, usubs, input_fname, grid_spacing_mm, plot_figures, varargin)
% function downsample_surface_timecoursesAll(exp, usubs, input_fname, grid_spacing_mm, plot_figures, varargin)
% 
% Wrapper function for downsample_surface_timecourses.m.
% 
% Example:
% exp = 'pitch_localizer_monkey';
% usubs = 158;
% input_fname = 'motcorr_smooth286mm';
% grid_spacing_mm = 2.86*0.5;
% plot_figures = 1;
% downsample_surface_timecoursesAll(exp, usubs, input_fname, grid_spacing_mm, plot_figures, 'runtypes', {'pitchloc2_combined_split'}, 'runnum', 1:2, 'monkey')

for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp,usubs(i),'func');
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(runtypes)
        
        runnum = read_runs(exp,usubs(i),runtypes{j},varargin{:});
        if optInputs(varargin, 'runnum');
            runnum = varargin{optInputs(varargin,'runnum')+1};
        end
        
        for k = 1:length(runnum)
            fprintf('%s, subject %d, %s, run %d\n',exp,usubs(i),runtypes{j},runnum(k));
            downsample_surface_timecourses(exp, usubs(i), runtypes{j}, runnum(k), input_fname, grid_spacing_mm, plot_figures, varargin{:})
        end
    end
end