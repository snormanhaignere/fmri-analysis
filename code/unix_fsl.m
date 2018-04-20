function [y, full_command] = unix_fsl(version, command, varargin)

if optInputs(varargin, 'nohup')
    command = ['nohup ' command ' &'];
end

switch version
    case '4.1'
        full_command = ['/usr/bin/fsl4.1-' command];
    case '5.0'
        full_command = ['/usr/bin/fsl5.0-' command];
    otherwise
        full_command = command;
end

fprintf([full_command '\n']);
[~,y] = unix(full_command);