function ds = readdata(file)
    warning('Off', 'stats:dataset:subsasgn:DefaultValuesAdded');
    fprintf('Temporarily turning off warnings for default value addition to dataset\n');
    fid = fopen(file, 'r');             % Open text file
    temp = textscan(fid, '%s', 'Delimiter', '\n');   % Scan till first return
    rows = size(temp{1}, 1) - 1;  % Count number of data rows (substract header row)
    fprintf('Number of data points: %i\n', rows);
    frewind(fid);
    flat = textscan(fid, '%s'); % Read all fields
    cols = size(flat{1}, 1) / (rows+1); % Count number of observables
    fprintf('Number of input/observables: %i\n', cols);
    ds = dataset(); % generate empty dataset
    names = regexprep(flat{1}(1:cols), '\.','_'); %generate input/observable names, replace . with _

    ds.(names{1}) = flat{1}{cols+1};
    for i = 2 : cols
        switch flat{1}{cols+i}
            case 'TRUE'
                ds.(names{i}) = true;
            case 'FALSE'
                ds.(names{i}) = false;
            case 'NA'
                ds.(names{i}) = NaN;
            otherwise
                ds.(names{i}) = str2double(flat{1}{cols+i});
        end
    end
    for j = 3 : rows
        ds.(names{1})(j-1,:) = flat{1}{(j-1)*cols+1};
        for i = 2 : cols
            switch flat{1}{(j-1)*cols+i}
                case 'TRUE'
                    ds.(names{i})(j-1,:) = true;
                case 'FALSE'
                    ds.(names{i})(j-1,:) = false;
                case 'NA'
                    ds.(names{i})(j-1,:) = NaN;
                otherwise
                    ds.(names{i})(j-1,:) = str2double(flat{1}{(j-1)*cols+i});
            end
        end
    end
    warning('On', 'stats:dataset:subsasgn:DefaultValuesAdded');    
end