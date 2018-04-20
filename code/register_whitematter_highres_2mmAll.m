function register_whitematter_highres_2mmAll(exp,usubs,varargin)

for i = 1:length(usubs)
    fprintf('sub %d\n',usubs(i));
    register_whitematter_highres_2mm(exp,usubs(i),varargin{:});
end