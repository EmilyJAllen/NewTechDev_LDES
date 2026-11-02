function [load_kW, info] = ldes_load_profile(cfg)
% LDES_LOAD_PROFILE  Annual 15-minute data-center load (kW) from a CSV.
%
% cfg.file        CSV path (required).
% cfg.column      load column: name, index, or [] (last numeric column)
% cfg.units       'MW' (default) or 'kW'
% cfg.scaleTo_MW  [] = use as-is; a number rescales the profile so its
%                 cfg.scaleBasis ('mean' default, or 'peak') equals this data-center size
% cfg.dt_h        time step (0.25 h)
%
% Accepted lengths: 35040/35136 (a year of 15-min steps), 8760/8784 (hourly, repeated
% to 15-min), 96 (one average day, tiled to 365 days), 672 (one week, tiled).

if ~isfield(cfg,'file'),        cfg.file = '';        end
if ~isfield(cfg,'column'),      cfg.column = [];      end
if ~isfield(cfg,'units'),       cfg.units = 'MW';     end
if ~isfield(cfg,'scaleTo_MW'),  cfg.scaleTo_MW = [];  end
if ~isfield(cfg,'scaleBasis'),  cfg.scaleBasis = 'mean'; end
if ~isfield(cfg,'dt_h'),        cfg.dt_h = 0.25;      end

nYear = round(8760/cfg.dt_h);
if isempty(cfg.file) || exist(cfg.file, 'file') ~= 2
    error('Load CSV "%s" not found. Put it next to ldes_main.m or give the full path.', cfg.file);
else
    [x, col] = ldes_read_csv_column(cfg.file, cfg.column);
    units = cfg.units;  info.source = cfg.file;  info.column = col;
end

n = numel(x);  x = x(:);
stepsPerDay = round(24/cfg.dt_h);
if n == nYear || n == nYear + stepsPerDay
    % full year (leap year kept as is)
elseif n == 8760 || n == 8784
    x = repelem_(x, round(1/cfg.dt_h));
    info.note = 'hourly data repeated to 15-min';
elseif n == stepsPerDay
    x = repmat(x, 365, 1);  info.note = 'one day tiled to a year';
elseif n == 7*stepsPerDay
    x = repmat(x, 53, 1);  x = x(1:nYear);  info.note = 'one week tiled to a year';
else
    error('Load profile has %d rows; expected %d (year of 15-min), 8760 (hourly), %d (day) or %d (week).', ...
          n, nYear, stepsPerDay, 7*stepsPerDay);
end

switch upper(units)
    case 'MW', load_kW = 1000*x;
    case 'KW', load_kW = x;
    otherwise, error('Unknown units "%s"', units);
end

if ~isempty(cfg.scaleTo_MW)
    switch lower(cfg.scaleBasis)
        case 'peak', ref = max(load_kW);
        otherwise,   ref = mean(load_kW);
    end
    load_kW = load_kW * (1000*cfg.scaleTo_MW/ref);
end
info.n = numel(load_kW);  info.mean_MW = mean(load_kW)/1e3;  info.peak_MW = max(load_kW)/1e3;
info.annual_MWh = sum(load_kW)*cfg.dt_h/1e3;
end

function y = repelem_(x, k)
    y = reshape(repmat(x(:)', k, 1), [], 1);
end
