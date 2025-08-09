function outcomeTab = insulinOutcomes(obj, interval_)
if nargin < 2 || isempty(interval_)
    indices_ = obj.intervalToIndices([0, 0]);
else
    indices_ = obj.intervalToIndices(interval_);
end

if nargin < 2 || isempty(interval_)
    unitOfDay = 24 * 60;
else
    unitOfDay = mod(diff(interval_*60), 24*60);
    if unitOfDay == 0
        unitOfDay = 24 * 60;
    end
end

row = 0;
row = row + 1;
outcomeName{row, 1} = 'Basal Insulin (U)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
elseif days(obj.durationOutcome) < 2
    outcomeValue(row, 1) = obj.getTDBasal(indices_);
else
    outcomeName{row, 1} = [outcomeName{row, 1}(1:end - 1), '/day)'];
    outcomeValue(row, 1) = obj.getTDBasal(indices_) / (sum(~isnan(obj.basalInsulin)) * obj.stepTime / unitOfDay);
end

row = row + 1;
outcomeName{row, 1} = 'Basal Insulin Percentage (%)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
else
    outcomeValue(row, 1) = 100 * obj.getTDBasal(indices_) / obj.getTDD(indices_);
end

row = row + 1;
outcomeName{row, 1} = 'Bolus Insulin (U)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
elseif days(obj.durationOutcome) < 2
    outcomeValue(row, 1) = obj.getTDBolus(indices_);
else
    outcomeName{row, 1} = [outcomeName{row, 1}(1:end - 1), '/day)'];
    outcomeValue(row, 1) = obj.getTDBolus(indices_) / (sum(~isnan(obj.basalInsulin)) * obj.stepTime / unitOfDay);
end

row = row + 1;
outcomeName{row, 1} = 'Bolus Insulin Percentage (%)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
else
    outcomeValue(row, 1) = 100 * obj.getTDBolus(indices_) / obj.getTDD(indices_);
end

row = row + 1;
outcomeName{row, 1} = 'Total Insulin (U)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
elseif days(obj.durationOutcome) < 2
    outcomeValue(row, 1) = obj.getTDD(indices_);
else
    outcomeName{row, 1} = [outcomeName{row, 1}(1:end - 1), '/day)'];
    outcomeValue(row, 1) = obj.getTDD(indices_) / (sum(~isnan(obj.basalInsulin)) * obj.stepTime / unitOfDay);
end

outcomeTab = table(outcomeValue, 'RowNames', outcomeName);
outcomeTab.Properties.VariableNames = {obj.name};
outcomeTab.Properties.VariableDescriptions = {'Data Name'};
end
