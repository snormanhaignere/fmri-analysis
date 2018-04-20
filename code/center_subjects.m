function [v_centered, grandmean, offset_matrix] = center_subjects(v,si)

usubs = unique(si);
grandmean = nanmean(v,2);
offset_matrix = nan(size(v));
for i = 1:length(usubs)
    xi = si == usubs(i);
    offset_matrix(:,xi) = mean(v(:,xi),2)*ones(1,sum(xi)) - grandmean*ones(1,sum(xi));
end
v_centered = v - offset_matrix;