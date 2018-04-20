function obr = brainops(ibr,op,opval,varargin)

switch op
    
    case 'thr'
        
        ibr(ibr==0) = NaN;
        obr = zeros(size(ibr));
        obr(ibr>opval) = 1;
        fprintf('Selecting %d voxels\n',sum(obr(:))); drawnow;
        
    case 'negthr'
        
        ibr(ibr==0) = NaN;
        obr = zeros(size(ibr));
        obr(ibr<opval) = 1;
        fprintf('Selecting %d voxels\n',sum(obr(:))); drawnow;
        
    case 'power'
        
        obr = ibr.^opval;
        
    case 'nvox'
        
        % ensure masked values have lowest values
        ibr(ibr==0) = -inf;
                
        % sort input values
        [~, sortinds] = sort(ibr(:),'descend');
        
        % select desired voxels
        inds = sortinds(opval(1):opval(2));

        % check that masked values were never chosen
        if any(ibr(inds)==-inf);
            fprintf('Error in brainops: Selected vertices should not be zero.');
            keyboard;
        end
        
        % create output brain with chosen voxels set to 1
        obr = zeros(size(ibr));
        obr(inds) = 1;
        
    case 'perc'
        
        % ensure masked values have lowest values
        ibr(ibr==0) = -inf;
        
        % sort input values
        [~, sortinds] = sort(ibr(:),'descend');

        % select desired voxels
        masksize = varargin{1};
        nvox = [max(round(masksize*opval(1)),1), round(masksize*opval(2))];
        inds = sortinds(nvox(1):nvox(2));
        fprintf('Selecting voxels %d through %d\n',nvox(1),nvox(2)); drawnow;

        % check that masked values were never chosen
        if any(ibr(inds)==-inf);
            fprintf('Error in brainops: Selected vertices should not be zero.');
            keyboard;
        end

        % create output brain with chosen voxels set to 1
        obr = zeros(size(ibr));
        obr(inds) = 1;
        
end