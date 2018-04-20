function preproc(exp, us, runtype, r, varargin)

if optInputs(varargin, 'monkey')
    steps = {'mcflirt','slicetimer','bet_func','ipfsl','detrend'};
elseif optInputs(varargin, 'stop_after_bet')
    steps = {'mcflirt','slicetimer','bet_func'};
else
    steps = {'mcflirt','slicetimer','bet_func','susan','detrend','hpfilt'};
end

if optInputs(varargin, 'steps')
    steps = varargin{optInputs(varargin, 'steps')+1};
end

for i = 1:length(steps)
    switch steps{i}
        case 'mcflirt'
            mcflirt(exp,us,runtype,r,varargin{:});
        case 'slicetimer'
            slicetimer(exp,us,runtype,r,varargin{:});
        case 'bet_func'
            bet_func(exp,us,runtype,r,varargin{:});
        case 'ipfsl'
            ipfsl(exp,us,runtype,r,varargin{:});
        case 'susan'
            susan(exp,us,runtype,r,varargin{:});
        case 'hpfilt'
            hpfilt(exp,us,runtype,r,varargin{:});
        case 'detrend'
            detrend(exp,us,runtype,r,varargin{:});
        case 'detrend_v2' % detrends the slicetime corrected data rather than the smoothed data
            detrend_v2(exp,us,runtype,r,varargin{:});
        case 'resample_slicetimecorr2surface'
            resample_slicetimecorr2surface(exp,us,runtype,r,'rh',varargin{:});
            resample_slicetimecorr2surface(exp,us,runtype,r,'lh',varargin{:});
        case 'hpfilt_surf'
            hpfilt_surf(exp,us,runtype,r,'rh',varargin{:});
            hpfilt_surf(exp,us,runtype,r,'lh',varargin{:});
        case 'resample_detrend2surface'
            resample_detrend2surface(exp,us,runtype,r,'rh',varargin{:});
            resample_detrend2surface(exp,us,runtype,r,'lh',varargin{:});
        case 'smooth_surf'
            smooth_surf(exp,us,runtype,r,'rh',varargin{:});
            smooth_surf(exp,us,runtype,r,'lh',varargin{:});
        case 'combine_surf_runs_demean'
            rungroups = varargin{optInputs(varargin,'rungroups')+1};
            combine_surf_runs_demean(exp,us,runtype,rungroups,'rh',varargin{:});
            combine_surf_runs_demean(exp,us,runtype,rungroups,'lh',varargin{:});
        case 'combine_surf_runs_nulldemean'
            rungroups = varargin{optInputs(varargin,'rungroups')+1};
            combine_surf_runs_nulldemean(exp,us,runtype,rungroups,'rh',varargin{:});
            combine_surf_runs_nulldemean(exp,us,runtype,rungroups,'lh',varargin{:});
    end
end

% end
% if optInputs(varargin, 'mcflirt-slicetimecorr')
%     mcflirt(exp,us,runtype,r,varargin{:});
%     slicetimer(exp,us,runtype,r,varargin{:});
% elseif
%
% else
%     mcflirt(exp,us,runtype,r,varargin{:});
%     slicetimer(exp,us,runtype,r,varargin{:});
%     if optInputs(varargin, 'overwrite_bet')
%         bet_func(exp,us,runtype,r,varargin{:},'overwrite');
%     else
%     end
%
%     if optInputs(varargin, 'stop_after_bet')
%         return;
%     end
%
%     if optInputs(varargin, 'monkey')
%     else
%     end
%
%     if optInputs(varargin, 'monkey')
%     else
%     end
%
% end