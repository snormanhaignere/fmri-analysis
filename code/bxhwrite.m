function bxhwrite(mrstruct,writefile)

% bxhwrite(mrstruct,writefile)
% 
% writes mrstructs to file using writemr and bxh2analyze

warning('OFF');
writefile = strrep(writefile, '.nii.gz','');

% remove all similar files
if ~isempty(dir([writefile '.*']))
    unix(['rm ' writefile '.*']);
end

% writemr
bxhfile = [writefile '_img.bxh'];
writemr(mrstruct,bxhfile,'BXH','OVERWRITE');

% bxh2analyze
unix(['bxh2analyze --niigz ' bxhfile ' ' writefile]);
unix(['rm ' writefile '_img*']);

warning('ON');