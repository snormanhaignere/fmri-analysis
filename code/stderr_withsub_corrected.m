function [sem, data_zeromean_rows] = stderr_withsub_corrected(data)
% [sem, data_subnorm] = stderr_withsub_corrected(data)
% 
% Computes standard errors for columns of the matrix data, removing mean differences
% between rows. This is often useful when there are consistent relative differences 
% between conditions across subjects, but large overall differenes in the mean across conditions.
% Because removing differences in the row means will always remove variance, even if there 
% is no "true" difference in the row means from which the data are sampled, a correction term
% is needed to make the standard errors unbiased. 
% 
% Data is a [M x N x ..] matrix. 
% sem is a [1 x N x ...] matrix with the standard error for each column of the input matrix. 
% data_zeromean_rows is the data matrix with demeaned columns, from which the standard errors are computed
% 
% Example
% submean = 1000*rand(10000,1);
% condmean = 100*rand(1,2);
% data = submean * ones(1,2) + ones(10000,1) * condmean + 100*randn(10000,2);
% analytic_stderr = 100/sqrt(10000)
% across_subject_stderr = std(data)/sqrt(10000-1)
% within_subjects_stderr_corrected = stderr_withsub_corrected(data)
% [~,data_zeromean_rows] = stderr_withsub_corrected(data);
% within_subjects_stderr_uncorrected = std(data_zeromean_rows)/sqrt(10000-1)
% 
% Last Modified by Sam Norman-Haignere on 1/12/2015

% matrix dimensions
dims = size(data);
ndims = length(dims);

% number of independent measurements for each column
nsub = dims(1);

% number of columns
nconds = prod(dims(2:end));

% correction factor to make estimates unbiased
correction_factor = sqrt(1/(nconds*(nconds-1)/(nconds^2)));

% row means
row_means = sumdims(data, 2:ndims)/nconds;

% data with zero-mean rows
data_zeromean_rows = data - repmat(row_means, [1 dims(2:end)]);

% corrected standard error measure
sem = correction_factor * std(data_zeromean_rows,[],1)/sqrt(nsub-1);

function outMat = sumdims(inMat, dim)
% sums matrix inMat across dimensions in vector dim

for xx = 1:length(dim)
    %     curDim = dim(xx)
    %     outMat = mean(inMat, curDim)
    %     dim = setdiff(dim, curDim);
    %     dim(dim>curDim) = dim(dim>curDim)-1;
    inMat = nansum(inMat, dim(xx));
end
outMat = inMat;