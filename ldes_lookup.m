function [Pel, Pcool, PlossOpen, PlossClosed] = ldes_lookup(sur, E, scale)
% LDES_LOOKUP  Surrogate power flows (kW) at battery energy E (kWh above Tmin).
% scale = number of block-equivalents (battery capacity / one-block capacity).
% Linear interpolation on the uniform grid in sur; clamped at the grid ends.
if nargin < 3, scale = 1; end
Eb = E/scale;
n  = numel(sur.E_kWh);
s  = (Eb - sur.E_kWh(1))/sur.dE;
i  = min(max(floor(s), 0), n-2) + 1;
w  = min(max(s - (i-1), 0), 1);
Pel         = scale*(sur.Pel_kW(i)         + w.*(sur.Pel_kW(i+1)         - sur.Pel_kW(i)));
Pcool       = scale*(sur.Pcool_kW(i)       + w.*(sur.Pcool_kW(i+1)       - sur.Pcool_kW(i)));
PlossOpen   = scale*(sur.PlossOpen_kW(i)   + w.*(sur.PlossOpen_kW(i+1)   - sur.PlossOpen_kW(i)));
PlossClosed = scale*(sur.PlossClosed_kW(i) + w.*(sur.PlossClosed_kW(i+1) - sur.PlossClosed_kW(i)));
end
