function D = ifelse(A,B,C)

% D = ifelse(A,B,C)
% 
% if A is true, returns B, if not true, returns C
% 
% compresses if/else statement into a single function
% particularly useful for generating file strings with many paramaters
% much less retarded function name than old bool2ans...

if A
    D = B;
else
    D = C;
end
