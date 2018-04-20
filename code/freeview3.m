function freeview3(subjid, hemi, varargin)

% annotfile = ['/mindhive/nklab/u/svnh/freesurfer/fsaverage/label/' hemi '.stp-pm2al-5lab-recolor.annot'];
inflated_file = ['~/freesurfer/' subjid '/surf/' hemi '.inflated_dist2'];
if ~exist(inflated_file,'file')
    inflated_file = ['~/freesurfer/' subjid '/surf/' hemi '.inflated'];
end
if strcmp(hemi,'rh')
    cam = [200 20 20];
elseif strcmp(hemi,'lh')
    cam = [-20 -20 15];
end

% subdirectory with useful functions
if isempty(strfind(path, 'export_fig'))
    addpath(genpath([pwd '/export_fig']));
end


% keyboard;

% [vertices, label, ct] = read_annotation(annot_file);
% % label(label == 0) = NaN;
% % ct.numEntries = ct.numEntries-1;
% % ct.struct_names = ct.struct_names(2:end);
% % ct.table = ct.table(2:end,:);
% ct.table(1,1:3) = 1*ones(1,3);
% % label(label==0) = 2960768;
% % label = round(rand(size(label))*1000);
% annot_file_new = strrep(annot_file, '.annot','_formatted.annot');
% delete(annot_file_new);
% write_annotation(annot_file_new, vertices, label, ct);
% [~,annot_file_new_abs] = unix(['ls -d ' annot_file_new]);

% freeview_call = ['freeview --surface ' inflated_file ':annot=' strtrim(annot_file_new_abs) ':edgethickness=0' ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' --viewsize 800 600 --zoom 1.2 &'];
% unix_freesurfer_version('5.3.0', freeview_call);


freeview_call = ['freeview --surface ' inflated_file];

if optInputs(varargin, 'overlay');
    overlay_file = varargin{optInputs(varargin, 'overlay')+1};
    freeview_call = [freeview_call ':overlay=' overlay_file];
end

if optInputs(varargin, 'overlay_threshold')
    overlay_threshold = varargin{optInputs(varargin, 'overlay_threshold')+1};
    freeview_call = [freeview_call ':overlay_threshold=' num2str(overlay_threshold(1)) ',' num2str(overlay_threshold(2)) ',' num2str(overlay_threshold(3))];
end

if optInputs(varargin, 'piecewise')
    freeview_call = [freeview_call ':overlay_method=piecewise'];
end

if optInputs(varargin, 'label');
    label_file = varargin{optInputs(varargin, 'label')+1};
    [~,label_file_abs] = unix(['ls -d ' label_file]);
    freeview_call = [freeview_call ':label=' strtrim(label_file_abs) ':label_outline=true'];
end

if optInputs(varargin, 'annot');
    annot_file = varargin{optInputs(varargin, 'annot')+1};
    %     [~,annot_file_abs] = unix(['ls -d ' annot_file]);
    tmpfile = [pwd '/tmp_' hemi  '.annot'];
    copyfile(annot_file, tmpfile, 'f'); 
    freeview_call = [freeview_call ':annot=' tmpfile];
    [~,label,ct] = read_annotation(annot_file);
    figure;
    set(gcf, 'Position', [0 0 200 600]);
    labels_used = unique(label(label~=0));
    ypositions = nan(1,length(labels_used));
    for i = 1:length(labels_used)
        xi = find(ct.table(:,5) == labels_used(i));
        try
        barh(xi-1, 1, 'FaceColor', ct.table(xi,1:3)/255)
        catch
            keyboard
        end
        if i == 1;
            hold on;
        end
        ypositions(i) = xi-1;
    end
    ylim([0, max(ypositions)+1]);
    set(gca, 'YTick', 1:max(ypositions));
    ylabel('Cluster Number');
    title('colorbar');
    if optInputs(varargin, 'screenshot')
        export_fig(strrep(varargin{optInputs(varargin, 'screenshot')+1},'.png','_ct.eps'),'-eps','-nocrop');
    end
end

x = [];
% if optInputs(varargin, 'piecewise')
%     x = [x, ' --overlay_method=piecewise '];
% end

if optInputs(varargin, 'screenshot')
    x = [x, ' --screenshot ' varargin{optInputs(varargin, 'screenshot')+1} '_' hemi '.png'];
else
    x = [x, ' &'];
end
freeview_call = [freeview_call ':edgethickness=0' ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' --viewport 3d --viewsize 800 600 --zoom 1.2 ' x];
unix_freesurfer_version('5.3.0', freeview_call);


% if optInputs(varargin, 'label')
%     label_file = varargin{optInputs(varargin, 'label')+1};
%     [~,label_file_abs] = unix(['ls -d ' label_file]);
%     unix_freesurfer5p3(['freeview --surface ' inflated_file ':overlay=' overlay_file ':overlay_threshold=' num2str(fminmidmax(1)) ',' num2str(fminmidmax(2)) ',' num2str(fminmidmax(3)) ':label=' strtrim(label_file_abs) ':label_outline=true' ':edgethickness=0' ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' --viewsize 800 600 --zoom 1.2 &']);
% elseif optInputs(varargin, 'annot')
%     unix_freesurfer5p3(['freeview --surface ' inflated_file ':overlay=' overlay_file ':overlay_threshold=' num2str(fminmidmax(1)) ',' num2str(fminmidmax(2)) ',' num2str(fminmidmax(3)) ':annot=' annotfile ':label_outline=true' ':edgethickness=0' ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' --zoom 1.2 &']);
% else
%     unix_freesurfer5p3(['freeview --surface ' inflated_file ':overlay=' overlay_file ':overlay_threshold=' num2str(fminmidmax(1)) ',' num2str(fminmidmax(2)) ',' num2str(fminmidmax(3)) ':edgethickness=0' ' -cam Azimuth ' num2str(cam(1)) ' Roll ' num2str(cam(2)) ' Elevation ' num2str(cam(3)) ' --viewsize 800 600 --zoom 1.2 &']); %1400 700
% end