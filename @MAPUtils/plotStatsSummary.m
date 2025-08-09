function plotStatsSummary(fig, data, varargin)

% support multiple datasets
if iscell(data)
    M = numel(data);
else
    M = 1;
    data = {data};
end

armsName_ = cell(1, M);
for m = 1:M
    armsName_{m} = data{m}(1).name;
end

% Get options
for nVar = 1:2:length(varargin)
    switch lower(varargin{nVar})
        case 'armsname'
            armsName_ = varargin{nVar+1};
    end
end

yMax = 0.87;
yMin = min(0.17, max(0.04, (4-M)*(0.17 - 0.03)/(4-2) + 0.03));
lenY = 0.38;

for m = 1:M
    name = armsName_{m};
    if contains(name, '#')
        name = name(1:strfind(name, '#')-1);
    end
    
    bgStats = uibuttongroup(fig);
    bgStats.Position = [0.84, yMin + (M - m) * min(lenY, (yMax - yMin)/M) + max(0, (yMax - yMin - M * lenY)/2), 0.15, min(lenY, (yMax - yMin)/M)]; % x, y, lenX, lenY
    bgStats.BackgroundColor = [0.9, 0.9, 0.9];
    bgStats.Tag = ['Outcomes: ', name];
    
    if length(data{m}) == 1
        fields = cell(length(data{m})+1, 2);
    else
        fields = cell(length(data{m})+3, 2);
    end
    
    table = nan(length(data{m}), 4);
            
    rownames = cell(1, numel(data{m}));
    for k = 1:numel(rownames)
        if contains(data{m}(k).name, '#')
            rownames{k} = ['D', data{m}(k).name(strfind(data{m}(k).name, '#'):end)];
        elseif ~isempty(data{m}(k).patientName)
            rownames{k} = data{m}(k).patientName;
        else
            rownames{k} = num2str(k);
        end
    end
    
    [~, idxInTab] = sort(rownames);
    
    for k = 1:numel(rownames)
        rowname = rownames{k};
        
        table(idxInTab(k), :) = [ ...
            MAPUtils.sensorTime(data{m}(k).glucoseInterp, data{m}(k).stepTime), ...
            MAPUtils.timeIn(data{m}(k).glucoseInterp, 3.9, 10+1e-5, data{m}(k).stepTime), ...
            MAPUtils.timeIn(data{m}(k).glucoseInterp, 0, 3.9, data{m}(k).stepTime), ...
            nansum(data{m}(k).bolusInsulin) + nansum(data{m}(k).basalInsulin) * data{m}(k).stepTime / 60, ...
            ];
        
        rowdata = { ...
            sprintf('%4.1f%%', table(idxInTab(k), 1)); ...
            sprintf('%4.1f%%', table(idxInTab(k), 2)); ...
            sprintf('%4.1f%%', table(idxInTab(k), 3))};
        if table(idxInTab(k), 4) > 0
            rowdata(end+1) = {sprintf('%4.1fU', table(idxInTab(k), 4))};
        end
        
        fields(idxInTab(k)+1, :) = {rowname, rowdata};
    end
    
    if length(data{m}) > 1
        idx2remove = find(table(:, 1) < 10 | isnan(table(:, 1)));
        table(idx2remove, :) = [];
        fields(idx2remove+1, :) = [];
        
        rowMean = { ...
            sprintf('%4.1f%%', nanmean(table(:, 1))); ...
            sprintf('%4.1f%%', nanmean(table(:, 2))); ...
            sprintf('%4.1f%%', nanmean(table(:, 3)))};
        if nanmean(table(:, 4)) > 0
            rowMean(end+1) = {sprintf('%4.1fU', nanmean(table(:, 4)))};
        end
        fields(end-1, :) = {'Mean', rowMean};
        
        rowMean = { ...
            sprintf('%4.1f%%', nanstd(table(:, 1))); ...
            sprintf('%4.1f%%', nanstd(table(:, 2))); ...
            sprintf('%4.1f%%', nanstd(table(:, 3)))};
        if nanstd(table(:, 4)) > 0
            rowMean(end+1) = {sprintf('%4.1fU', nanstd(table(:, 4)))};
        end
        fields(end, :) = {'SD', rowMean};
    end
    
    if nanmean(table(:, 4)) > 0
        fields(1, :) = {'Title', {'Sensor'; 'Target'; '< 3.9 (mmol/L)'; 'Insulin'}};
    else
        fields(1, :) = {'Title', {'Sensor'; 'Target'; '< 3.9 (mmol/L)'}};
    end
    fields = flipud(fields);
    
    MAPUtils.plotTable( ...
        bgStats, ...
        fields, ...
        'TitleFont', 10, ...
        'YPad', 0.005, ...
        'YOffset', 0.01, ...
        'XOffset', 0.01, ...
        'XPad', 0.005, ...
        'XTitleLength', 0.2, ...
        'XTitleFont', 10, ...
        'TableFont', 10);
end
