function sheet = toSheet(obj, varargin)
sheet = cell(0);

row = 1;
sheet{row, 1} = 'Intervention';
sheet{row, 2} = obj.name;

detailsTab = obj.details();
for out = 1:size(detailsTab, 1)
    row = row + 1;
    sheet{row, 1} = detailsTab.Row{out};
    if iscell(detailsTab{out, 1})
        if isnumeric(detailsTab{out, 1}{1})
            sheet{row, 2} = sprintf('%6.3f', detailsTab{out, 1}{1});
        else
            sheet{row, 2} = detailsTab{out, 1}{1};
        end
    else
        sheet{row, 2} = sprintf('%6.3f', detailsTab{out, 1});
    end
end

outcomeTab = obj.outcomes(varargin{:});
for out = 1:size(outcomeTab, 1)
    row = row + 1;
    sheet{row, 1} = outcomeTab.Row{out};
    sheet{row, 2} = sprintf('%6.3f', outcomeTab{out, 1});
end

len = length(obj.timeStamp);
row = row + 1;
col = 0;

if nansum(obj.timeStamp) > 0
    col = col + 1;
    sheet{row, col} = 'Time Stamp (dd-mmm-yyyy HH:MM)';
    sheet{row+1, col} = 'timeStamp';
    sheet(row+2:row+len+1, col) = cellstr(datestr(obj.timeStamp, 'dd-mmm-yyyy HH:MM'));
    
    col = col + 1;
    sheet{row, col} = 'Time (min)';
    sheet{row+1, col} = 'time';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5d', c), num2cell(obj.time), 'UniformOutput', false);
end

if nansum(obj.sensorGlucose) > 0
    col = col + 1;
    sheet{row, col} = 'Sensor Glucose (mmol/L)';
    sheet{row+1, col} = 'sensorGlucose';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.3f', c), num2cell(obj.sensorGlucose), 'UniformOutput', false);
end

if nansum(obj.glucoseInterp) > 0
    col = col + 1;
    sheet{row, col} = 'Interpolated Glucose (mmol/L)';
    sheet{row+1, col} = 'glucoseInterp';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.3f', c), num2cell(obj.glucoseInterp), 'UniformOutput', false);
end

if nansum(obj.basalInsulin) > 0
    col = col + 1;
    sheet{row, col} = 'Insulin Basal (U/h)';
    sheet{row+1, col} = 'basalInsulin';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.3f', c), num2cell(obj.basalInsulin), 'UniformOutput', false);
end

if nansum(obj.basalInjection) > 0
    col = col + 1;
    sheet{row, col} = 'Insulin Basal (U)';
    sheet{row+1, col} = 'basalInjection';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.3f', c), num2cell(obj.basalInjection), 'UniformOutput', false);
end

if nansum(obj.bolusInsulin) > 0
    col = col + 1;
    sheet{row, col} = 'Insulin Bolus (U)';
    sheet{row+1, col} = 'bolusInsulin';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.3f', c), num2cell(full(obj.bolusInsulin)), 'UniformOutput', false);
end

if nansum(obj.carbs) > 0
    col = col + 1;
    sheet{row, col} = 'Meals Carbohydrates (g)';
    sheet{row+1, col} = 'carbs';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.1f', c), num2cell(obj.carbs), 'UniformOutput', false);
end

if nansum(obj.bolusGlucagon) > 0
    col = col + 1;
    sheet{row, col} = 'Glucagon Bolus (U)';
    sheet{row+1, col} = 'bolusGlucagon';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.1f', c), num2cell(obj.bolusGlucagon), 'UniformOutput', false);
end

if nansum(obj.pumpBasals) > 0
    col = col + 1;
    sheet{row, col} = 'Pump Basal Insulin (U/h)';
    sheet{row+1, col} = 'pumpBasals';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.3f', c), num2cell(obj.pumpBasals), 'UniformOutput', false);
end

if nansum(obj.carbFactors) > 0
    col = col + 1;
    sheet{row, col} = 'Carb Factor (g/U)';
    sheet{row+1, col} = 'carbFactors';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.1f', c), num2cell(obj.carbFactors), 'UniformOutput', false);
end

if nansum(obj.steps) > 0
    col = col + 1;
    sheet{row, col} = 'Steps (#)';
    sheet{row+1, col} = 'steps';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.1f', c), num2cell(obj.steps), 'UniformOutput', false);
end

if nansum(obj.heartrate) > 0
    col = col + 1;
    sheet{row, col} = 'Heart Rate (#/min)';
    sheet{row+1, col} = 'heartrate';
    sheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%5.1f', c), num2cell(obj.heartrate), 'UniformOutput', false);
end
end