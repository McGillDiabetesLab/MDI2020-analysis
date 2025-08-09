function out = getMealInfo(obj, meal, patientPerPatient)
if numel(obj) > 1
    if nargin < 3
        patientPerPatient = false;
    end
    if nargin < 2
        out_ = arrayfun(@(c)(c.getMealInfo()), obj);
    else
        out_ = arrayfun(@(c)(c.getMealInfo(meal)), obj);
    end
    if patientPerPatient
        out = out_;
    else
        out = struct();
        for p = 1:numel(out_)
            if ~isfield(out, 'ID')
                out.ID = repmat(string(obj(p).patientName), size(string(out_(p).name)));
            else
                out.ID = [out.ID; repmat(string(obj(p).patientName), size(string(out_(p).name)))];
            end
            fn = fieldnames(out_(p));
            for i = 1:length(fn)
                if ~isfield(out, fn{i})
                    out.(fn{i}) = out_(p).(fn{i});
                else
                    out.(fn{i}) = [out.(fn{i}); out_(p).(fn{i})];
                end
            end
        end
    end
    return;
end

if nargin < 2
    meal = 'all';
end

meals = {'breakfast', 'lunch', 'dinner'};
if isnumeric(meal)
    if meal < 4
        meal = meals{meal};
    else
        meal = 'all';
    end
end

meal = lower(meal);

if strcmpi(meal, 'all')
    out = struct();
    for m = 1:numel(meals)
        out_ = obj.getMealInfo(meals{m});
        fn = fieldnames(out_);
        for i = 1:length(fn)
            if ~isfield(out, fn{i})
                out.(fn{i}) = out_.(fn{i});
            else
                out.(fn{i}) = [out.(fn{i}); out_.(fn{i})];
            end
        end
    end
    return;
end

time_ = obj.time + obj.startTime;
interval = find(mod(time_, 1440) >= obj.([meal, 'Interval'])(1) & mod(time_, 1440) <= obj.([meal, 'Interval'])(2));
if ~isnan(obj.outcomeStartTime)
    interval(obj.time(interval) < obj.outcomeStartTime-obj.startTime) = [];
end
if ~isnan(obj.outcomeEndTime)
    interval(obj.time(interval) >= obj.outcomeEndTime-obj.startTime) = [];
end
if ~isempty(obj.carbs)
    idxMealAll = interval(obj.carbs(interval) > obj.snacksThreshold);
    if obj.onlyCountBolusedMeals
        bolusThresh = [-1; 0; 1; 2];
        for idx = idxMealAll(:)'
            idxBolusThresh = idx + bolusThresh;
            idxBolusThresh(idxBolusThresh <= 0) = [];
            idxBolusThresh(idxBolusThresh > length(obj.bolusInsulin)) = [];
            if nansum(obj.bolusInsulin(idxBolusThresh)) == 0
                idxMealAll(idxMealAll == idx) = NaN;
            end
        end
        idxMealAll(isnan(idxMealAll)) = [];
    end
else
    idxMealAll= [];
end

if isempty(idxMealAll)
    out = struct( ...
        'name', string(meal), ...
        'idx', NaN, ...
        'timeStamp', NaT, ...
        'amount', NaN, ...
        'glucose', NaN, ...
        'carbFactor', NaN, ...
        'AUC1h', NaN, ...
        'AUC2h', NaN, ...
        'AUC3h', NaN, ...
        'AUC4h', NaN, ...
        'incGlucose1h', NaN, ...
        'incGlucose2h', NaN, ...
        'incGlucose3h', NaN, ...
        'incGlucose4h', NaN, ...
        'incGlucoseMax', NaN, ...
        'hypos1h_4h', NaN, ...
        'data', ...
        struct( ...
        'name', obj.name, ...
        'time', [], ...
        'startTime', [], ...
        'duration', [], ...
        'timeStamp', [], ...
        'sensorGlucose', [], ...
        'glucoseInterp', [], ...
        'pumpBasals', [], ...
        'basalInsulin', [], ...
        'bolusInsulin', [], ...
        'carbFactors', [], ...
        'carbs', [], ...
        'carbsActual', [], ...
        'treats', [], ...
        'tdd', [], ...
        'iobBolus', []));
    return;
end

% regroup meals within 1 hour
veryCloseMealsIdx = [(flipud(diff(flipud(idxMealAll))) < -1 * 60 / obj.stepTime); true];
idxMeal = idxMealAll(veryCloseMealsIdx);
carbs = obj.carbs;
if any(~veryCloseMealsIdx)
    for n = find(~veryCloseMealsIdx)'
        carbs(idxMealAll(n)) = 0;
        nextMealIdx = n - 1 + find(veryCloseMealsIdx(n:end) == 1, 1);
        carbs(idxMealAll(nextMealIdx)) = carbs(idxMealAll(nextMealIdx)) + obj.carbs(idxMealAll(n));
    end
end

% only analyze bigger meal
if obj.onlyCountBiggerMeal
    idxMealsToRemove = [];
    
    for n = idxMeal(:)'
        idxOfMealsInSameInterval = find(idxMeal - n < diff(obj.([meal, 'Interval']))/obj.stepTime & idxMeal - n > 0) ;
        
        if isempty(idxOfMealsInSameInterval)
            continue;
        end
        
        idxOfMealsInSameInterval = [n; idxMeal(idxOfMealsInSameInterval)];
        
        [~, idxBiggerMeal] = max(carbs(idxOfMealsInSameInterval));
        idxMealsToRemove = [idxMealsToRemove; setdiff(idxOfMealsInSameInterval, idxOfMealsInSameInterval(idxBiggerMeal))];
    end
    
    if ~isempty(idxMealsToRemove)
        idxMeal = setdiff(idxMeal, unique(idxMealsToRemove));
    end
end

out.name = string(repmat(meal, size(idxMeal)));
out.idx = idxMeal;
out.timeStamp = datetime(obj.timeStamp(idxMeal), 'ConvertFrom', 'datenum');
out.amount = carbs(idxMeal);
out.glucose = obj.glucoseInterp(idxMeal);
out.carbFactor = obj.carbFactors(idxMeal);

out.AUC1h = MAPUtils.AUC(obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(0:1 * 60 / obj.stepTime), length(time_))));
out.AUC2h = MAPUtils.AUC(obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(0:2 * 60 / obj.stepTime), length(time_))));
out.AUC3h = MAPUtils.AUC(obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(0:3 * 60 / obj.stepTime), length(time_))));
out.AUC4h = MAPUtils.AUC(obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(0:4 * 60 / obj.stepTime), length(time_))));

out.incGlucose1h = obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(1 * 60 / obj.stepTime), length(time_))) - obj.glucoseInterp(idxMeal);
out.incGlucose2h = obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(2 * 60 / obj.stepTime), length(time_))) - obj.glucoseInterp(idxMeal);
out.incGlucose3h = obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(3 * 60 / obj.stepTime), length(time_))) - obj.glucoseInterp(idxMeal);
out.incGlucose4h = obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(4 * 60 / obj.stepTime), length(time_))) - obj.glucoseInterp(idxMeal);

out.incGlucoseMax = obj.glucoseInterp(min(idxMeal+ones(size(idxMeal))*(0:4 * 60 / obj.stepTime), length(time_))) - obj.glucoseInterp(idxMeal);
if size(out.incGlucoseMax, 2) == 1
    out.incGlucoseMax = max(out.incGlucoseMax);
else
    out.incGlucoseMax = max(out.incGlucoseMax, [], 2);
end

if nansum(obj.treats) > 0
    if size(idxMeal, 1) == 1
        out.hypos1h_4h = nansum(obj.treats(min(idxMeal+(1 * 60 / obj.stepTime:4 * 60 / obj.stepTime), length(obj.treats))) > 0);
    else
        out.hypos1h_4h = nansum(obj.treats(min(idxMeal+ones(size(idxMeal))*(1 * 60 / obj.stepTime:4 * 60 / obj.stepTime), length(obj.treats))) > 0, 2);
    end
else
    out.hypos1h_4h = nan(numel(out.name), 1);
end

offset = obj.mealInfoSnapshot / obj.stepTime;
for k = length(idxMeal):-1:1
    out.data(k, :).name = obj.name;
    out.data(k, :).time = (0:1:(offset(2) - offset(1)))'*obj.stepTime;
    out.data(k, :).startTime = mod(obj.time(out.idx) + obj.startTime, 1440) + obj.mealInfoSnapshot(1);
    out.data(k, :).duration = minutes(diff(obj.mealInfoSnapshot));
    out.data(k, :).timeStamp = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.timeStamp(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    if nansum(obj.sensorGlucose) > 0
        out.data(k, :).sensorGlucose = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.sensorGlucose(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).sensorGlucose = [];
    end
    if nansum(obj.glucoseInterp) > 0
        out.data(k, :).glucoseInterp = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.glucoseInterp(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).glucoseInterp = [];
    end
    if nansum(obj.pumpBasals) > 0
        out.data(k, :).pumpBasals = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.pumpBasals(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).pumpBasals = [];
    end
    if nansum(obj.basalInsulin) > 0
        out.data(k, :).basalInsulin = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.basalInsulin(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).basalInsulin = [];
    end
    if nansum(obj.bolusInsulin) > 0
        out.data(k, :).bolusInsulin = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.bolusInsulin(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).bolusInsulin = [];
    end
    if nansum(obj.carbFactors) > 0
        out.data(k, :).carbFactors = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.carbFactors(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).carbFactors = [];
    end
    if nansum(obj.carbs) > 0
        out.data(k, :).carbs = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.carbs(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).carbs = [];
    end
    if nansum(obj.carbsActual) > 0
        out.data(k, :).carbsActual = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.carbsActual(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).carbsActual = [];
    end
    if nansum(obj.treats) > 0
        out.data(k, :).treats = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.treats(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).treats = [];
    end
    if nansum(obj.iobBolus) > 0
        out.data(k, :).iobBolus = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.iobBolus(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).iobBolus = [];
    end
    if nansum(obj.tdd) > 0
        out.data(k, :).tdd = [nan(1-min(idxMeal(k)+offset(1), 1), 1); obj.tdd(max(idxMeal(k)+offset(1), 1):min(idxMeal(k)+offset(2), length(obj.time))); nan(max(idxMeal(k)+offset(2), length(obj.time))-length(obj.time), 1)];
    else
        out.data(k, :).tdd = [];
    end
end
end
