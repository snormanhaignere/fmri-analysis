function unix_freesurfer_version(freesurfer_version, command)

% 2017-11-13: Updated to use openmind, Sam NH

if exist('/cm/shared/openmind/freesurfer', 'dir')
    setupstr = ['export FREESURFER_HOME=/cm/shared/openmind/freesurfer/' freesurfer_version  '; source ${FREESURFER_HOME}/SetUpFreeSurfer.sh; export SUBJECTS_DIR=~/freesurfer; '];
else
    switch computer
        case 'GLNXA64'
            setupstr = ['export FREESURFER_HOME=/software/Freesurfer/' freesurfer_version  '; source ${FREESURFER_HOME}/SetUpFreeSurfer.sh; export SUBJECTS_DIR=~/freesurfer; '];
        case 'MACI64'
            setupstr = 'export FREESURFER_HOME=/Applications/freesurfer; source ${FREESURFER_HOME}/SetUpFreeSurfer.sh; export SUBJECTS_DIR=/Users/svnh2/Desktop/projects/freesurfer; ';
    end
end
unix([setupstr command]);