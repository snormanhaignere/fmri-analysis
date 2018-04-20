function inds = strmatchMANY(str2find,strpool)

inds = [];
for strindex = 1:length(str2find)
    inds = [inds;strmatch(str2find{strindex},strpool,'exact')];
end