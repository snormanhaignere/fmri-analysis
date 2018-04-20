function reptextMANYFAST(fidr,fidw,list1,list2,varargin)

while( true )
    readline = fgets(fidr);
    if readline == -1
        break;
    end
    writeline = readline;
    for listindex = 1:length(list1)
        writeline = strrep(writeline, list1{listindex}, list2{listindex}); % replace
    end
    fprintf(fidw,'%s',writeline);
end