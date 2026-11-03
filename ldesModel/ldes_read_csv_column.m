function [x, colName] = ldes_read_csv_column(file, column)
% LDES_READ_CSV_COLUMN  Read one numeric column from a CSV (MATLAB and Octave).
%   column: header name (case-insensitive), column index, or [] = the LAST column
%           that is mostly numeric (timestamps and text columns are skipped).
% A header row is detected automatically.
if nargin < 2, column = []; end
txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');
lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));
% Quoted fields may contain the delimiter (e.g. "12:00 AM Jan 01, 2026"): blank those
% delimiters out so the columns don't shift.
for it = 1:5
    new = regexprep(lines, '"([^"]*?)[,;]([^"]*)"', '"$1 $2"');
    if isequal(new, lines), break; end
    lines = new;
end
delim = ',';
if ~contains_(lines{1}, ',') && contains_(lines{1}, ';'), delim = ';'; end
if ~contains_(lines{1}, delim) && contains_(lines{1}, sprintf('\t')), delim = sprintf('\t'); end

first = strsplit(lines{1}, delim);
hasHeader = any(isnan(str2double(first)) & ~cellfun(@isempty, strtrim(first)));
if hasHeader
    hdr = strtrim(strrep(first, '"', ''));  body = lines(2:end);
else
    hdr = arrayfun(@(k) sprintf('col%d', k), 1:numel(first), 'UniformOutput', false);  body = lines;
end
nc = numel(hdr);
cells = cellfun(@(s) strsplit(s, delim), body, 'UniformOutput', false);
cells = cellfun(@(c) [c, repmat({''}, 1, max(nc - numel(c), 0))], cells, 'UniformOutput', false);
M = cellfun(@(c) str2double(strrep(c(1:nc), '"', '')), cells, 'UniformOutput', false);
M = vertcat(M{:});

if isempty(column)
    numFrac = mean(~isnan(M), 1);
    k = find(numFrac > 0.95, 1, 'last');
    if isempty(k), error('No numeric column found in %s', file); end
elseif ischar(column)
    k = find(strcmpi(hdr, column), 1);
    if isempty(k), error('Column "%s" not found in %s. Columns: %s', column, file, strjoin(hdr, ', ')); end
else
    k = column;
end
x = M(:, k);  colName = hdr{k};
if any(isnan(x))
    warning('%d non-numeric/blank rows in column "%s" filled by interpolation', nnz(isnan(x)), colName);
    good = ~isnan(x);  idx = (1:numel(x))';
    x(~good) = interp1(idx(good), x(good), idx(~good), 'nearest', 'extrap');
end
end

function tf = contains_(s, pat)
    tf = ~isempty(strfind(s, pat));
end
