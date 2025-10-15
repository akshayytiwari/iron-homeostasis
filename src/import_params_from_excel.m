function params = import_params_from_excel(xlsxfile, sheet)
% IMPORT_PARAMS_FROM_EXCEL Read parameter spreadsheet exported from your Excel.
%   params = import_params_from_excel('parameters_2mo_09Nov2024_supplement_foods_tuned.xlsx')
% Expects a sheet with two columns: name | value (or parameter table similar to your CSV).
% If your spreadsheet layout differs, adapt the reading or export to CSV with (name,value).

if nargin < 2, sheet = 1; end

% try to read as table (supports both xlsx and csv saved as .csv)
try
    T = readtable(xlsxfile,'Sheet',sheet,'ReadVariableNames',false);
catch
    error('Cannot read file %s. Make sure file exists in current folder or provide a CSV.', xlsxfile);
end

% Expect first column = param name, second = numeric value
names = T{:,1};
vals  = T{:,2};

params = struct();
for i=1:length(names)
    % make name safe for struct field (replace spaces, parentheses etc.)
    fname = matlab.lang.makeValidName(string(names{i}));
    params.(fname) = double(vals(i));
end

% Quick mapping for some convenience fields used in the ODE function.
% (If your excel already has those exact names, mapping is not necessary.)
% Example: if excel stores 'k_L' or 'k L', ensure consistent field name
% You can add more mapping here if necessary.
end
