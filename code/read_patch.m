function patch = read_patch(pname,varargin)

if ~strcmp(pname(end-3:end), '.asc')
    pname_asc = [pname '.asc'];
    if ~exist(pname_asc,'file') || optInputs(varargin,'overwrite');
        unix_freesurfer(['mris_convert -p ' pname ' ' pname_asc]);
    end
else
    pname_asc = pname;
end

fid = fopen(pname_asc, 'r');

% Dump first line
patch.h = fgets(fid);

% Nr of vertices and faces
S = zeros(1, 3);
S(1) = 1;
[S(2:3)] = fscanf(fid, '%d', 2);

% Read vertices and its indices
x = textscan(fid, '%d%s\n%f%f%f\n', S(2));
patch.vnums = sign(x{1}).*(abs(x{1})-1);
patch.vras = [x{3} x{4} x{5}];

% Read faces and its indices
x = textscan(fid, '%d\n%d%d%d\n', S(3));
patch.fnums = x{1};
patch.fverts = [x{2} x{3} x{4}];

fclose(fid);