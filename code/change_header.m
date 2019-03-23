function change_header(nii_oldheader, nii_newheader, varargin)

% This script makes it possible to change the resolution of a file. There are
% two changes. First the size of the voxels is changed (dx, dy, dz) and second
% the transformation matrices which take you from voxel to spatial coordinates
% and back are scaled up or down. This is done twice to modify both q-form and
% s-form matrices (see:
% http://gru.stanford.edu/doku.php/mrtools/coordinatetransforms).

% 2018-04-22: Last updated, fixed a bug which caused dx, dy, dz to only be
% updated when the original resolution was 1 mm.

old_header_file = strrep(nii_oldheader, '.nii.gz', '.xml');
new_header_file = strrep(nii_newheader, '.nii.gz', '.xml');
unix_fsl('5.0', ['fslhd -x ' nii_oldheader ' > ' old_header_file]);

fidr = fopen(old_header_file, 'r');
fidw = fopen(new_header_file, 'w');

while( true )
    readline = fgets(fidr);
    if readline == -1
        break;
    end
    writeline = readline;
    writeline = strrep(writeline, 'sto_ijk matrix = ','sto_ijk_matrix = ');
    if optInputs(varargin, 'falseheader')
        
        resolution_factor = varargin{optInputs(varargin, 'falseheader')+1};
        
        dims = {'dx', 'dy', 'dz'};
        for i = 1:length(dims)
            
            % find string of a format like dx = 'NUMBER';
            match = regexp(writeline, [dims{i} ' = \''[\d\.]*\'''], 'match'); 
            if ~isempty(match)
                assert(length(match)==1);
                
                % figure out the original resolution
                native_res = str2double(regexp(match{1}, '[\d\.]*', 'match'));
                
                % multiply by factor
                faked_res = native_res*resolution_factor;
                
                % replace with new number
                writeline = strrep(writeline, match{1}, regexprep(match{1}, '[\d\.]*', num2str(faked_res)));
            end
        end
        
        % q-form transformation matrix, which takes you from voxel coordinates
        % (i.e. indexes of the matrix) to coordinates of the coordinates in the
        % magnet, http://gru.stanford.edu/doku.php/mrtools/coordinatetransforms
        if ~isempty(strfind(writeline, 'qto_xyz_matrix = '));
            xi = strfind(writeline, '''');
            xi = xi(1)+1:xi(2)-1;
            vec = str2num(writeline(xi));
            mat = reshape(vec,4,4)';
            mat(1:3,1:3) = resolution_factor*mat(1:3,1:3);
            mat = mat';
            vec = mat(:)';
            writeline = [writeline(1:xi(1)-1), num2str(vec,'%.6f ') '''' writeline(end)];         
        end
        
        % s-form transformation matrix, which takes you from voxel coordinates
        % to coordinates of a standard space (e.g. MNI)
        if ~isempty(strfind(writeline,'sto_xyz_matrix = '));
            xi = strfind(writeline, '''');
            xi = xi(1)+1:xi(2)-1;
            vec = str2num(writeline(xi));
            mat = reshape(vec,4,4)';
            mat(1:3,1:3) = resolution_factor*mat(1:3,1:3);
            mat = mat';
            vec = mat(:)';
            writeline = [writeline(1:xi(1)-1), num2str(vec,'%.6f ') '''' writeline(end)];         
        end
        
        % inverse of qto_xyz_matrix
        if ~isempty(strfind(writeline, 'qto_ijk_matrix = '));
            xi = strfind(writeline, '''');
            xi = xi(1)+1:xi(2)-1;
            vec = str2num(writeline(xi));
            mat = reshape(vec,4,4)';
            mat(1:3,1:3) = (1/resolution_factor)*mat(1:3,1:3);
            mat = mat';
            vec = mat(:)';
            writeline = [writeline(1:xi(1)-1), num2str(vec,'%.6f ') '''' writeline(end)];         
        end
        
        % inverse of sto_xyz_matrix
        if ~isempty(strfind(writeline, 'sto_ijk_matrix = '));
            xi = strfind(writeline, '''');
            xi = xi(1)+1:xi(2)-1;
            vec = str2num(writeline(xi));
            mat = reshape(vec,4,4)';
            mat(1:3,1:3) = (1/resolution_factor)*mat(1:3,1:3);
            mat = mat';
            vec = mat(:)';
            writeline = [writeline(1:xi(1)-1), num2str(vec,'%.6f ') '''' writeline(end)];         
        end
                
    end
        
    
    %     for listindex = 1:length(list1)
    %         writeline = strrep(writeline, list1{listindex}, list2{listindex}); % replace
    %     end
    fprintf(fidw,'%s',writeline);
end
fclose(fidr); fclose(fidw);
copyfile(nii_oldheader, nii_newheader, 'f');
unix_fsl('5.0', ['fslcreatehd ' new_header_file ' ' nii_newheader]);

