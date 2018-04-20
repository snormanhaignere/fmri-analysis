function convert_monkey_struct(exp, usubs, varargin)

% directory setup
datadir = [params('rootdir') exp '/data/'];
fsl_version = read_fsl_version(exp,varargin{:});
freesurfer_version = read_freesurfer_version(exp,varargin{:});

for i = 1:length(usubs)
    
    
    %     monkey_name = read_monkey_name(usubs(i));
    %     struct_RAS = [params('rootdir') 'monkey_anatomicals/' monkey_name '/' monkey_name '.nii'];
    %     struct_LAS = [preproc_directory 'brain.nii.gz'];
    %     struct_RAS
    %     if ~exist(struct_LAS, 'file') || optInputs(varargin, 'overwrite')
    %         unix_freesurfer_version(freesurfer_version, ['mri_convert --out_orientation LAS --out_type nii '  struct_RAS '  ' struct_LAS]);
    %     end
    
    niidir = [datadir 'brain/nifti/usub' num2str(usubs(i)) '/'];
    if ~exist(niidir,'dir');
        mkdir(niidir);
    end
        
    % convert to niigz
    mgz_RAS = [params('rootdir') 'freesurfer/' exp '_us' num2str(usubs) '/mri/brain.mgz'];
    nii_RAS = [niidir 'struct_r1_falseheader.nii.gz'];
    if ~exist(nii_RAS,'file') || optInputs(varargin, 'overwrite');
        unix_freesurfer_version(freesurfer_version, ['mri_convert ' mgz_RAS ' ' nii_RAS]);
    end
        
    %     nii_LAS = [niidir 'struct_r1_LAS_falseheader.nii.gz'];
    %     if ~exist(nii_LAS,'file') || optInputs(varargin, 'overwrite');
    %         unix_fsl(fsl_version, ['fslreorient2std '  nii_RAS '  ' nii_LAS]);
    %         falseheader = [niidir 'struct_r1_LAS_falseheader.xml'];
    %         unix_fsl(fsl_version, ['fslhd -x ' nii_LAS ' > ' falseheader]);
    %         unix_fsl(fsl_version, ['fslcreatehd ' falseheader ' ' nii_LAS]);
    %     end
    
    preproc_directory = [params('rootdir') '/' exp '/analysis/preprocess/usub' num2str(usubs(i)) '/struct_r1/'];
    if ~exist(preproc_directory, 'dir')
        mkdir(preproc_directory);
    end
    
    % copy to preprocessing directory
    preproc_brain = [preproc_directory 'brain.nii.gz'];
    if ~exist(preproc_brain,'file') || optInputs(varargin, 'overwrite');
        copyfile(nii_RAS, preproc_brain, 'f');
    end
    
    % copy over brain mask
    mask_mgz_RAS = [params('rootdir') 'freesurfer/' exp '_us' num2str(usubs) '/mri/brainmask.mgz'];
    mask_nii_RAS = [preproc_directory 'brain_mask.nii.gz'];
    if ~exist(mask_nii_RAS,'file') || optInputs(varargin, 'overwrite');
        unix_freesurfer_version(freesurfer_version, ['mri_convert ' mask_mgz_RAS ' ' mask_nii_RAS]);
    end
    
    
    
    %     nii_LAS_falseheader = [niidir 'struct_r1_LAS_falseheader.nii.gz'];
    %     if ~exist(nii_LAS_falseheader,'file') || optInputs(varargin, 'overwrite');
    %         %                 correctheader = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS_correctheader.xml'];
    %         %                 falseheader1 = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS_falseheader1.xml'];
    %         %                 falseheader2 = [niidir runtypes{j} '_r' num2str(runnum(k)) '_LAS_falseheader2.xml'];
    %         %                 unix_fsl(fsl_version, ['fslhd -x ' nii_LAS ' > ' correctheader]);
    %         %                 fid1 = fopen(correctheader, 'r');
    %         %                 fid2 = fopen(falseheader1,'w');
    %         %                 keyname{1} = 'dx = ''1'''; value{1} = 'dx = ''2.8571''';
    %         %                 keyname{2} = 'dy = ''1'''; value{2} = 'dy = ''2.8571''';
    %         %                 keyname{3} = 'dz = ''1'''; value{3} = 'dz = ''2.8571''';
    %         %                 reptextMANYFAST(fid1, fid2, keyname, value);
    %         %                 fclose(fid1); fclose(fid2);
    %         %                 unix(['cat ' falseheader1 ' | grep -v ''^[ 	]*#'' | grep -v ''^[ 	]*$'' > ' falseheader2])
    %         %
    %         if exist(nii_LAS_falseheader, 'file');
    %             delete(nii_LAS_falseheader);
    %         end
    %         copyfile(nii_LAS, nii_LAS_falseheader, 'f');
    %         unix_fsl(fsl_version, ['fslcreatehd 96 96 33 130 2.8571 2.8571 2.8571 3.4 0 0 0 4 ' nii_LAS_falseheader]);
    %     end
end