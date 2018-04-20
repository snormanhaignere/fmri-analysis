function [runnum, dicomid, scanid] = read_runs(exp, us, runtype, varargin)

% function [runnum, dicomid, scanid] = read_runs(exp, us, runtype, varargin)
% returns all of the run numbers, dicom numbers, and scan numbers for a given
% runtype, subject and experiment
% data is read from file svnh/EXP/data/brain/runorders/RUNTYPE_usUS.txt

if strcmp(exp, 'amusia')
    runnum = 1;
    dicomid = NaN;
    scanid = 1;
    return;
end

runordersdir = [params('rootdir') exp '/data/brain/runorders/'];
analysisdir = [params('rootdir') exp '/analysis/'];
fid = fopen([runordersdir runtype '_us' num2str(us) '.txt'],'r');

if fid == -1;
    runnum = [];
    dicomid = [];
    return;
end


switch exp
    
    case {'pitch_natural', 'pitch_speech'}
        
        tmp1 = fscanf(fid,'%d%d');
        tmp2 = reshape(tmp1,2,length(tmp1)/2)';
        runnum = tmp2(:,1)';
        dicomid = tmp2(:,2)';
        scanid = [];
        fclose(fid);
        
    case {'mspec'}
        
        dicomid = fscanf(fid,'%d');
        fclose(fid);
        runnum = 1:length(dicomid);
        
        % subjects 3, 5, and 6 had bad first runs
        if any(us == [35 37 12]) && strcmp(runtype,'main');
            runnum = runnum(2:end);
            dicomid = dicomid(2:end);
        end
        scanid = [];
                
  otherwise
        
        x = fscanf(fid,'%d%d%d');
        y = reshape(x,3,length(x)/3)';
        
        if optInputs(varargin,'scans');
            scans = varargin{ optInputs(varargin,'scans') + 1 };
            inds = any(repmat(y(:,3), 1, length(scans)) == repmat(scans, size(y,1),1),2);
            y = y(inds,:);
        end
        
        runnum = y(:,1)';
        dicomid = y(:,2)';
        scanid = y(:,3)';
        
        
        fclose(fid);
        
end

if optInputs(varargin, 'oddruns')
    xi = mod(1:length(runnum),2) == 1;
    runnum = runnum(xi);
    dicomid = dicomid(xi);
    scanid = scanid(xi);
elseif optInputs(varargin, 'evenruns')
    xi = mod(1:length(runnum),2) == 0;
    runnum = runnum(xi);
    dicomid = dicomid(xi);
    scanid = scanid(xi);
elseif optInputs(varargin, 'runnum')
    xi = ismember(runnum,varargin{optInputs(varargin, 'runnum')+1});
    runnum = runnum(xi);
    dicomid = dicomid(xi);
    scanid = scanid(xi);
elseif optInputs(varargin, 'first_Nruns')
    N = length(runnum);
    firstN_runs = varargin{optInputs(varargin, 'first_Nruns')+1};
    runnum = runnum(1:min(N,firstN_runs));
    dicomid = dicomid(1:min(N,firstN_runs));
    scanid = scanid(1:min(N,firstN_runs));
end


% for roi analysis ignore runs with not enough TRs
if strcmp(exp, 'pitch_dp_v2') && optInputs(varargin, 'roi_test_runs') && us == 3
    bad_runs = ismember(runnum, [1 2]);
    runnum(bad_runs) = [];
    dicomid(bad_runs) = [];
    scanid(bad_runs) = [];
end

if optInputs(varargin, 'no_contrast_matching_for_run_selection')
    return;
end

fwhm = read_smooth(exp, us, runtype, varargin{:});

% only return runs that have a particular contrast
% only relevant/possible after running first level analysis
if optInputs(varargin,'contrast')
    
    contr = varargin{optInputs(varargin,'contrast')+1};
    model = varargin{optInputs(varargin,'contrast')+2};
    
    flasubdir = [analysisdir 'fla/usub' num2str(us) '/'];
    
    runnum_contrast = [];
    dicomid_contrast = [];
    scanid_contrast = [];
    for r = runnum
        featdir = [flasubdir runtype '_r' num2str(r) '_' model '_' num2str(100*fwhm, '%.0f') 'mm.feat/'];
        fid = fopen([featdir 'contrastnames.txt'],'r');
        try
            tmp = textscan(fid,'%s\n'); fclose(fid); cname = tmp{1};
        catch
            keyboard
        end
        ind = strmatch(contr,cname,'exact');
        if ~isempty(ind);
            runnum_contrast = [runnum_contrast r];  %#ok<AGROW>
            dicomid_contrast = [dicomid_contrast dicomid(r==runnum)];  %#ok<AGROW>
            if ~isempty(scanid)
                scanid_contrast = [scanid_contrast scanid(r==runnum)]; %#ok<AGROW>
            end
        end
    end
    
    runnum = runnum_contrast;
    dicomid = dicomid_contrast;
    scanid = scanid_contrast;
    
end

% 
% dicomid = fscanf(fid,'%d');
% fclose(fid);
% runnum = 1:length(dicomid);
% 
% % subjects 3, 5, and 6 had bad first runs
% if any(s == [3 5 6]) && strcmp(runtype,'main');
%     runnum = runnum(2:end);
%     dicomid = dicomid(2:end);
% end