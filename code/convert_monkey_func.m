function convert_monkey_func(exp, usubs, varargin)

% function convert_monkey_func(exp, usubs, varargin)
% changes the header file information to make the voxels 2.8571 bigger than
% those actually measured, voxel size change facilitates analysis with
% packages designed for humans

% directory setup
datadir = [params('rootdir') exp '/data/'];
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

for i = 1:length(usubs)
    
    runtypes = read_runtypes(exp, usubs(i), 'raw_func', varargin{:});
    if optInputs(varargin,'runtypes')
        runtypes = varargin{optInputs(varargin,'runtypes')+1};
    end
    
    for j = 1:length(runtypes)
        
        try
            [runnum, dicomid, scanid] = read_runs(exp,usubs(i),runtypes{j},varargin{:});
        catch me
            print_error_message(me);
            keyboard;
        end
        if optInputs(varargin,'runnum')
            x = varargin{optInputs(varargin,'runnum')+1};
            [~,xi] = intersect(runnum,x);
            runnum = runnum(xi);
            dicomid = dicomid(xi);
            scanid = scanid(xi);
        end
        
        for k = 1:length(runnum)
            
            niidir = [datadir 'brain/nifti/usub' num2str(usubs(i)) '/'];
            if ~exist(niidir,'dir');
                mkdir(niidir);
            end
            
            fprintf('%s, sub %d, %s, run %d\n',exp,usubs(i),runtypes{j},runnum(k));
            
            sphinx_flag = '';
            if optInputs(varargin,'sphinx_correct')
                sphinx_flag = ' --sphinx ';
            end
            
            % copy file to nifti directory and gzip
            raw_directory = [datadir 'brain/raw' ...
                '/usub' num2str(usubs(i)) '_scan' num2str(scanid(k))];
            
            switch usubs(i)
                case {157, 158, 170} % Harvard format
                    nii_scanner_oldloc =  [raw_directory ...
                        '/' strrep(sprintf('%3d', dicomid(k)),' ','0') '/f.nii'];
                case {372, 373} % NIH format
                    single_run_directory = mydir(...
                        raw_directory, ['_E' num2str(dicomid(k)) '_']);
                    assert(length(single_run_directory) == 1);
                    single_run_directory = single_run_directory{1};
                    nii_scanner_oldloc = ...
                        [raw_directory '/' single_run_directory '/MRIm0001'];
                otherwise
                    error(['Need to specify whether us %d' ...
                        'is Harvard or NIH format'], usubs(i));
            end
            
            nii_scanner_newloc = [niidir runtypes{j} '_r' num2str(runnum(k)) '.nii.gz'];
            if ~exist(nii_scanner_newloc,'file') || optInputs(varargin, 'overwrite');
                unix_freesurfer_version(freesurfer_version, ['mri_convert ' sphinx_flag ' ' nii_scanner_oldloc ' ' nii_scanner_newloc]);
            end
            
            % create the false header with voxels multiplied by 2.8571
            nii_falseheader = [niidir runtypes{j} '_r' num2str(runnum(k)) '_falseheader.nii.gz'];
            if ~exist(nii_falseheader,'file') || optInputs(varargin, 'overwrite');
                switch usubs(i)
                    case {170, 372, 373}
                        change_header(nii_scanner_newloc, nii_falseheader, 'falseheader',2);
                    case {157, 158}
                        change_header(nii_scanner_newloc, nii_falseheader, 'falseheader',2.8571);
                    otherwise
                        error('Need to specify the voxel size of the false header');
                end
            end
            
            % file with orientation information more easily read by FSL
            % nii_RAS = [niidir runtypes{j} '_r' num2str(runnum(k)) '_falseheader_reorient.nii.gz'];
            % if ~exist(nii_RAS,'file') || optInputs(varargin, 'overwrite');
            %   unix_fsl(fsl_version, ['fslreorient2std '  nii_falseheader '  ' nii_RAS]);
            % end
            
            % copy file to preprocessing directory
            preprocdir = [params('rootdir') exp '/analysis/preprocess/usub' num2str(usubs(i)) '/' runtypes{j} '_r' num2str(runnum(k)) '/'];
            if ~exist(preprocdir,'dir');
                mkdir(preprocdir);
            end
            nii_raw = [preprocdir 'raw.nii.gz'];
            if ~exist(nii_raw,'file') || optInputs(varargin, 'overwrite');
                copyfile(nii_falseheader, nii_raw, 'f');
            end
            
            %             x = load_nii(nii_RAS);
            %
            %             nii_LAS_falseheader = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS_falseheader.nii.gz'];
            %             if ~exist(nii_LAS_falseheader,'file') || optInputs(varargin, 'overwrite');
            %                 correctheader = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS_correctheader.xml'];
            %                 falseheader1 = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS_falseheader1.xml'];
            %                 falseheader2 = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS_falseheader2.xml'];
            %                 unix_fsl(fsl_version, ['fslhd -x ' nii_LAS ' > ' correctheader]);
            %                 fid1 = fopen(correctheader, 'r');
            %                 fid2 = fopen(falseheader1,'w');
            %                 keyname{1} = 'dx = ''1'''; value{1} = 'dx = ''2.8571''';
            %                 keyname{2} = 'dy = ''1'''; value{2} = 'dy = ''2.8571''';
            %                 keyname{3} = 'dz = ''1'''; value{3} = 'dz = ''2.8571''';
            %                 reptextMANYFAST(fid1, fid2, keyname, value);
            %                 fclose(fid1); fclose(fid2);
            %                 unix(['cat ' falseheader1 ' | grep -v ''^[ 	]*#'' | grep -v ''^[ 	]*$'' > ' falseheader2])
            %                 %
            %                 if exist(nii_LAS_falseheader, 'file');
            %                     delete(nii_LAS_falseheader);
            %                 end
            %                 copyfile(nii_LAS, nii_LAS_falseheader, 'f');
            %                 %                 unix_fsl(fsl_version, ['fslcreatehd 96 96 33 130 2.8571 2.8571 2.8571 3.4 0 0 0 4 ' nii_LAS_falseheader]);
            %             end
        end
    end
end