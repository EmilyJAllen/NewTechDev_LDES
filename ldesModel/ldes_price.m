function [price, info] = ldes_price(cfg, n, dt_h, weights)
% LDES_PRICE  Electricity price ($/kWh) for each of n steps of dt_h hours.
%
% cfg.mode
%   'pjm2025'  Synthetic 15-min real-time price series CALIBRATED to the PJM 2025
%              State of the Market statistics (Monitoring Analytics, 2025 SOM, Vol. 1):
%                load-weighted average real-time LMP  $50.73/MWh
%                median                               $38.11/MWh
%                standard deviation                   $58.74/MWh
%              The series reproduces those three numbers exactly (load-weighted with the
%              data-center load). Its TIMING (daily/seasonal/weekend pattern, multi-day
%              weather swings) is stylized, not historical: the SOM gives the distribution,
%              not the hourly series. For real timing, export PJM real-time hourly LMPs
%              (Data Miner 2 feed rt_hrl_lmps) and use 'csv'.
%   'csv'      column cfg.column of cfg.file; hourly or 15-min; cfg.units '$/MWh' or '$/kWh'
%   'flat'     cfg.flat $/kWh every step
%   'vector'   cfg.vector (length n), $/kWh
%
% weights: load (any units) for the load-weighted mean/std in 'pjm2025' (default flat).
if nargin < 4 || isempty(weights), weights = ones(n,1); end
info = struct('mode', cfg.mode);
switch lower(cfg.mode)
    case 'pjm2025'
        price = pjmSynth(cfg, n, dt_h, weights(:));
        w = weights(:)/sum(weights);
        info.mean_MWh   = 1000*sum(w.*price);
        info.median_MWh = 1000*median(price);
        info.std_MWh    = 1000*sqrt(sum(w.*(price - sum(w.*price)).^2));
        info.max_MWh    = 1000*max(price);
    case 'flat'
        price = cfg.flat*ones(n,1);
    case 'csv'
        x = ldes_read_csv_column(cfg.file, cfg.column);
        k = round(1/dt_h);
        if numel(x)*k == n || numel(x)*k == n + 24*k     % hourly prices
            x = reshape(repmat(x(:)', k, 1), [], 1);
        end
        x = x(1:min(end,n));
        if numel(x) ~= n, error('Price column has %d rows, need %d.', numel(x), n); end
        price = x(:);
        if ~isfield(cfg,'units') || strcmpi(cfg.units, '$/MWh'), price = price/1000; end
    case 'vector'
        price = cfg.vector(:);
        if numel(price) ~= n, error('Price vector has %d entries, need %d.', numel(price), n); end
    otherwise
        error('Unknown price mode "%s"', cfg.mode);
end
end

function price = pjmSynth(cfg, n, dt_h, wload)
    % Stylized shape -> standardized score z (median 0) -> monotone transform fitted so the
    % load-weighted mean and std hit the targets; the median is exact by construction.
    rs = rng;  rng(cfg.seed);                          % fixed seed: same series every run
    t   = (0:n-1)'*dt_h;
    hod = mod(t, 24);   doy = t/24;   day = floor(doy);
    dow = mod(day + cfg.startDow, 7);                  % 0 = Monday
    bump = @(h, c, s) exp(-0.5*(min(abs(h - c), 24 - abs(h - c))/s).^2);
    win  = max(cos(2*pi*(doy - 20)/365), 0);           % winter weight (peaks ~Jan 20)
    sum_ = max(-cos(2*pi*(doy - 20)/365), 0);          % summer weight (peaks ~Jul 21)
    daily = 0.8*bump(hod, 18.5, 2.5) + 0.6*win.*bump(hod, 8, 1.5) ...
          + 0.8*sum_.*bump(hod, 16.5, 3) - 0.7*bump(hod, 3.5, 3);
    season  = 0.5*cos(4*pi*(doy - 20)/365);            % high Jan & Jul, low Apr & Oct
    weekend = -0.5*(dow >= 5);
    nd = day(end) + 1;                                 % multi-day weather swings, AR(1) by day
    dz = zeros(nd,1);  e = randn(nd,1);
    for d = 2:nd, dz(d) = 0.75*dz(d-1) + sqrt(1 - 0.75^2)*e(d); end
    phi = exp(-dt_h/4);  fz = zeros(n,1);  e = randn(n,1);   % intra-day noise, ~4 h memory
    for k = 2:n, fz(k) = phi*fz(k-1) + sqrt(1 - phi^2)*e(k); end
    zr = daily + season + weekend + 0.9*dz(day + 1) + 0.6*fz;
    rng(rs);
    z  = (zr - median(zr))/std(zr);

    w   = wload/sum(wload);  med = cfg.median_MWh;
    P   = @(a) med*exp(abs(a(1))*z + abs(a(2))*max(z,0).^2);
    mom = @(x) [sum(w.*x), sqrt(sum(w.*(x - sum(w.*x)).^2))];
    err = @(a) sum(((mom(P(a)) - [cfg.mean_MWh cfg.std_MWh])./[cfg.mean_MWh cfg.std_MWh]).^2);
    a   = fminsearch(err, [0.6 0.05], optimset('TolX',1e-10,'TolFun',1e-14,'MaxFunEvals',4000,'MaxIter',4000));
    price = P(a)/1000;                                  % $/kWh
end
