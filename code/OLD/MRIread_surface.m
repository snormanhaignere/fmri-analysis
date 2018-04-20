function S = MRIread_surface(surface_rh, surface_lh)
% function S = MRIread_surface(surface_rh, surface_lh)
% 
% Wrapper function for the matlab script "MRIread" used to read Freesurfer image files.
% Designed to be used with surface data, and returns a [M x N] matrix, with one column per surface vertex.
% Each vertex can have multiple numbers representing different features of that point.
% Assumes there is data from both the left and right hemisphere.
% 
% Example
% surface_rh = '/mindhive/nklab/u/svnh/fmri-analysis/test_data/rh.naturalsound_example_subject.mgz';
% surface_lh = '/mindhive/nklab/u/svnh/fmri-analysis/test_data/lh.naturalsound_example_subject.mgz';
% S = MRIread_surface(surface_rh,surface_lh);
% 
% Last modified by Sam Norman-Haignere on 1/5/2014

% directory this file is contained in
source_directory = strrep(which('MRIread_surface.m'),'MRIread_surface.m',''); 

% directory of matlab scripts for manipulating freesurfer files
addpath([source_directory '/fs']);

% read in surface data
rh = MRIread(surface_rh);
lh = MRIread(surface_lh);
S = [squeeze(rh.vol)', squeeze(lh.vol)'];