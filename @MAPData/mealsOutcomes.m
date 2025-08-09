function outcomeTab = mealsOutcomes(obj, interval_)
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
outcomeName{row, 1} = 'Total Carbs (g)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
elseif days(obj.durationOutcome) < 2
    outcomeValue(row, 1) = obj.getCarbs(indices_);
else
    outcomeName{row, 1} = [outcomeName{row, 1}(1:end-1) '/day)'];
    outcomeValue(row, 1) = obj.getCarbs(indices_) / (sum(~isnan(obj.basalInsulin)) * obj.stepTime / unitOfDay);
end

[amountTreats, nbrTreats] = obj.getTreats(indices_);
row = row + 1;
outcomeName{row, 1} = 'Total Treats (g)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
elseif days(obj.durationOutcome) < 2
    outcomeValue(row, 1) = amountTreats;
else
    outcomeName{row, 1} = [outcomeName{row, 1}(1:end-1) '/day)'];
    outcomeValue(row, 1) = amountTreats / (sum(~isnan(obj.basalInsulin)) * obj.stepTime / unitOfDay);
end

row = row + 1;
outcomeName{row, 1} = 'Number of Treats (#)';
if isempty(obj.durationOutcome)
    outcomeValue(row, 1) = NaN;
elseif days(obj.durationOutcome) < 2
    outcomeValue(row, 1) = nbrTreats;
else
    outcomeName{row, 1} = [outcomeName{row, 1}(1:end-1) '/day)'];
    outcomeValue(row, 1) = nbrTreats / (sum(~isnan(obj.basalInsulin)) * obj.stepTime / unitOfDay);
end

if nargin < 2 || isempty(interval_)
    meals = {'breakfast', 'lunch', 'dinner'};
    for m = 1:length(meals)
        info = obj.getMealInfo(meals{m});
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' (#)'];
        if isempty(obj.durationOutcome)
            outcomeValue(row, 1) = NaN;
        else
            if days(obj.durationOutcome) < 2
                outcomeValue(row, 1) = nansum(info.idx > 0);
            else
                outcomeName{row, 1} = [outcomeName{row, 1}(1:end-1) '/day)'];
                outcomeValue(row, 1) = nansum(info.idx > 0) / (sum(~isnan(obj.basalInsulin)) * obj.stepTime / unitOfDay);;
            end
        end
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' Meal (g)'];
        outcomeValue(row, 1) = nanmean(info.amount);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' carb Factor (g/U)'];
        outcomeValue(row, 1) = nanmean(info.carbFactor);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 1h iAUC (h.mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.AUC1h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 2h iAUC (h.mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.AUC2h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 3h iAUC (h.mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.AUC3h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 4h iAUC (h.mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.AUC4h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 1h inc increase (mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.incGlucose1h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 2h inc increase (mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.incGlucose2h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 3h inc increase (mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.incGlucose3h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' 4h inc increase (mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.incGlucose4h);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' max inc increase (1h->4h) (mmol/L)'];
        outcomeValue(row, 1) = nanmean(info.incGlucoseMax);
        
        row = row + 1;
        outcomeName{row, 1} = [meals{m}, ' hypos (1h->4h) (#)'];
        outcomeValue(row, 1) = nanmean(info.hypos1h_4h);
    end
end

outcomeTab = table(outcomeValue, 'RowNames', outcomeName);
outcomeTab.Properties.VariableNames = {obj.name};
outcomeTab.Properties.VariableDescriptions = {'Data Name'};
end
