function outcomeTab = glucoseOutcomes(obj, interval_)
if nargin < 2 || isempty(interval_)
    indices_ = obj.intervalToIndices([0, 0]);
else
    indices_ = obj.intervalToIndices(interval_);
end

row = 0;

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time < 2.8 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time < 50 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(-inf, 2.8, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time < 3.0 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time < 54 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(-inf, 3.0, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time < 3.3 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time < 60 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(-inf, 3.3, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time < 3.9 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time < 70 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(-inf, 3.9, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time 3.9-10.0 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time 70-180 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(3.9, 10.0+1e-5, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time 3.9-7.8 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time 70-140 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(3.9, 7.8+1e-5, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time > 7.8 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time > 140 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(7.8, inf, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time > 10.0 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time > 180 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(10.0+1e-5, inf, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time > 13.9 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time > 250 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(13.9, inf, indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Time > 16.7 mmol/L (%)';
else
    outcomeName{row, 1} = 'Time > 300 mg/dL (%)';
end
outcomeValue(row, 1) = obj.getTimeIn(16.7, inf, indices_);

row = row + 1;
outcomeName{row, 1} = 'Low Glucose Index';
outcomeValue(row, 1) = obj.getLBGI(indices_);

row = row + 1;
outcomeName{row, 1} = 'High Glucose Index';
outcomeValue(row, 1) = obj.getHBGI(indices_);

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Mean sensor glucose (mmol/L)';
    outcomeValue(row, 1) = nanmean(obj.sensorGlucose(indices_));
else
    outcomeName{row, 1} = 'Mean sensor glucose (mg/dL)';
    outcomeValue(row, 1) = nanmean(obj.sensorGlucose(indices_)) * 18.018;
end

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'SD of sensor glucose (mmol/L)';
    outcomeValue(row, 1) = nanstd(obj.sensorGlucose(indices_));
else
    outcomeName{row, 1} = 'SD of sensor glucose (mg/dL)';
    outcomeValue(row, 1) = nanstd(obj.sensorGlucose(indices_)) * 18.018;
end

row = row + 1;
outcomeName{row, 1} = 'CV of sensor glucose (%)';
outcomeValue(row, 1) = 100 * nanstd(obj.sensorGlucose(indices_)) / nanmean(obj.sensorGlucose(indices_));

row = row + 1;
if strcmp(obj.units, 'uk')
    outcomeName{row, 1} = 'Starting glucose (mmol/L)';
else
    outcomeName{row, 1} = 'Starting glucose (mg/dL)';
end
if isempty(find(~isnan(obj.sensorGlucose(indices_)), 1))
    outcomeValue(row, 1) = NaN;
else
    if strcmp(obj.units, 'uk')
        outcomeValue(row, 1) = obj.sensorGlucose(indices_(find(~isnan(obj.sensorGlucose(indices_)), 1)));
    else
        outcomeValue(row, 1) = obj.sensorGlucose(indices_(find(~isnan(obj.sensorGlucose(indices_)), 1))) * 18.018;
    end
end

outcomeTab = table(outcomeValue, 'RowNames', outcomeName);
outcomeTab.Properties.VariableNames = {obj.name};
outcomeTab.Properties.VariableDescriptions = {'Data Name'};
end
