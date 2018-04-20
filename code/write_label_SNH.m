function write_label_SNH(lb,lname,varargin)

if ~exist(lname,'file') || optInputs(varargin,'overwrite');
    nverts = length(lb.vnums);
    fid = fopen(lname,'w');
    fprintf(fid,'%s\n',strtrim(lb.h));
    fprintf(fid,'%d\n',nverts);
    x = [double(lb.vnums) lb.vras double(zeros(size(lb.vnums)))];
    fprintf(fid,'%d  %.3f  %.3f  %.3f  %f\n',x');
    fclose(fid);
end