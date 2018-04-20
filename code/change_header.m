function change_header(nii_oldheader, nii_newheader, varargin)

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
        
        writeline = strrep(writeline, 'dx = ''1''', ['dx = ''' num2str(resolution_factor) '''']);
        writeline = strrep(writeline, 'dy = ''1''', ['dy = ''' num2str(resolution_factor) '''']);
        writeline = strrep(writeline, 'dz = ''1''', ['dz = ''' num2str(resolution_factor) '''']);
        
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

