function y = params(type,varargin)

% y = params(type,varargin)
% stores paramaters that are used throughout the analysis

switch type
    case 'smooth'
        error('No longer specifying smoothing kernel here.');
        
    case 'zthresh'
        y = 3.09;
                
    case 'tempderiv'
        if optInputs(varargin,'block')
            y = '0';
        elseif optInputs(varargin,'event')
            y = '0';
        end
        
        if optInputs(varargin,'mspec')
            y = '1';
        end
        
    case 'cluster_percmask'
        y = 0.25;

    case 'cluster_maxdiff'
        y = 0.75;
        
    case 'rootdir'
        x = strrep(which('params.m'),'scripts','');
        if strfind(x, '/mindhive/nklab/u/svnh/');
            y = '/mindhive/nklab/u/svnh/';
        elseif strfind(x, '/mindhive/mcdermott/dropbox/SYNCtoBOOTHS')
            y = '/mindhive/nklab/u/svnh/';
        elseif strfind(x,'/net/munin.mit.edu/vol/nklab/u/svnh/')
            y = '/mindhive/nklab/u/svnh/';
        elseif strfind(x, '/Users/svnh2/Desktop/');
            y = '/Users/svnh2/Desktop/projects/';
        elseif strfind(x, '/Users/svnh2/Dropbox (MIT)/mcdexp-svnh/');
            y = '/Users/svnh2/Desktop/projects/';
        elseif strfind(x, '/Volumes/nklab/u/svnh/');
            y = '/Volumes/nklab/u/svnh/';
        else
            fprintf('No valid root directory.\n'); drawnow;
            keyboard;
        end
    case 'hpcutoff'
        y = 250;
end

if optInputs(varargin,'str')
    y = strrep(num2str(y),'.','');
end

if optInputs(varargin,'num')
    y = str2num(y);
end