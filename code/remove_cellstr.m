function c_pruned = remove_cellstr(c, cellstr)

% remove_string_from_cell(c, str)
% 
% Removes set of strings from an arbitrary cell array
% 
% -- Example --
% 
% remove_cellstr({2,'leave','remove 1', 'remove 2'}, {'remove 1', 'remove 2'})

remove_element = false(size(c));
for i = 1:numel(c)
   if ischar(c{i}) && any(ismember(cellstr, c{i}))
       remove_element(i) = true;
   end
end
c_pruned = c(~remove_element);