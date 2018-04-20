function outputfiles = mydir(inputdir)

%%
filenames = dir(inputdir);

filecell = cell(length(filenames),1);
for fileindex = 1:length(filenames)
    filecell{fileindex} = filenames(fileindex).name;
end

% remove hidden files
outputfiles = filecell(setxor(strmatch('.',filecell), (1:length(filecell))));
