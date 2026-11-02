function [mask, events] = ldes_outages(n, dt_h, freq_per_yr, dur_h, mode, startHours)
% LDES_OUTAGES  Grid-outage mask for n time steps of dt_h hours (one year).
%   freq_per_yr  expected outages per year
%   dur_h        duration of each outage (h)
%   mode  'random'  number of outages ~ Poisson(freq_per_yr), start times uniform over
%                   the year (call rng(seed) before this for reproducibility)
%         'even'    round(freq_per_yr) outages spread evenly through the year
%         'given'   start times from startHours (h from Jan 1 00:00)
% events: [startStep, endStep] per outage (an outage running past Dec 31 is cut off).
if nargin < 5 || isempty(mode), mode = 'random'; end
nDur = max(1, round(dur_h/dt_h));
switch lower(mode)
    case 'random'
        k  = poissonDraw(freq_per_yr);
        s0 = sort(floor(rand(k,1)*n)) + 1;
    case 'even'
        k  = round(freq_per_yr);
        if k == 0 && freq_per_yr > 0
            warning('freq %.2f/yr rounds to 0 outages in "even" mode; use "random" with many realizations.', freq_per_yr);
        end
        s0 = round(((1:k)' - 0.5) * n/max(k,1)) + 1;
    case 'given'
        s0 = round(startHours(:)/dt_h) + 1;
    otherwise
        error('Unknown outage mode "%s"', mode);
end
mask = false(n,1);
events = zeros(numel(s0), 2);
for j = 1:numel(s0)
    e = min(s0(j) + nDur - 1, n);
    mask(s0(j):e) = true;
    events(j,:) = [s0(j) e];
end
end

function k = poissonDraw(lambda)
    % Knuth's method (no Statistics Toolbox needed); fine for lambda < ~50
    L = exp(-lambda);  k = 0;  pr = rand;
    while pr > L
        k = k + 1;  pr = pr*rand;
    end
end
