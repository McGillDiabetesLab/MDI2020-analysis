function outcomeTab = compareMean(obj, varargin)
arms_ = obj.arms;
interval_ = [];
set_ = 'outcomes';
baseline_ = '';
var_ = true;
normality_ = true;

% Get options
for nVar = 1:2:length(varargin)
    switch lower(varargin{nVar})
        case 'interval'
            interval_ = varargin{nVar+1};
        case 'arms'
            arms_ = varargin{nVar+1};
        case 'set'
            set_ = varargin{nVar+1};
        case 'baseline'
            baseline_ = varargin{nVar+1};
        case 'variability'
            var_ = varargin{nVar+1};
        case 'normality'
            normality_ = varargin{nVar+1};
    end
end

if ~isempty(baseline_) && numel(baseline_) ~= numel(arms_)
    error('[compareMean] Baseline should have the same size as arms')
end

if length(arms_) > 2
    error('[compareMean] Only support up to two arms.')
end

outcome = struct([]);
if strcmpi(obj.type, 'crossover')
    periods = {'P1_', 'P2_'};
    
    for p = 1:numel(obj.patients)
        if numel(intersect({obj.patients(p).data.name}, arms_)) ~= 2
            continue;
        end
        
        out_ = obj.patients(p).(set_)(varargin{:}, 'interval', interval_, 'arms', arms_);
        
        startDate_ = [obj.patients(p).getData(arms_{1}).startDate, obj.patients(p).getData(arms_{2}).startDate];
        
        [~, idxPeriod] = sort(startDate_);
        for aIdx = 1:length(arms_)
            outcome(p).([periods{idxPeriod(aIdx)}, arms_{aIdx}]) = out_.(arms_{aIdx});
        end
    end
    
    outcomeTab = table;
    if length(fieldnames(outcome)) ~= 4
        outcomeTab.NOT_ENOUGH_DATA = nan(size(out_.Properties.RowNames));
        outcomeTab.Properties.RowNames = out_.Properties.RowNames;
        return;
    end
    
    for aIdx = 1:length(arms_)
        data_ = [[outcome.(['P1_', arms_{aIdx}])], [outcome.(['P2_', arms_{aIdx}])]];
        outcomeTab.([arms_{aIdx}, '_Mean_n_', sprintf('%02d', size(data_, 2))]) = nanmean(data_, 2);
        if var_
            outcomeTab.([arms_{aIdx}, '_SD']) = nanstd(data_, 0, 2);
        end
    end
    
    for fn = fieldnames(outcome)'
        outcomeTab.([fn{1}, '_Mean_n_', sprintf('%02d', size([outcome.(fn{1})], 2))]) = nanmean([outcome.(fn{1})], 2);
        if var_
            outcomeTab.([fn{1}, '_SD']) = nanstd([outcome.(fn{1})], 0, 2);
        end
    end
    
    outcomeTab.(['Seq1_', arms_{1}, '_', arms_{2}]) = nanmean([outcome.(['P1_', arms_{1}])]-[outcome.(['P2_', arms_{2}])], 2);
    outcomeTab.(['Seq2_', arms_{2}, '_', arms_{1}]) = nanmean([outcome.(['P1_', arms_{2}])]-[outcome.(['P2_', arms_{1}])], 2);
    
    [~, p] = ttest2( ...
        ([outcome.(['P1_', arms_{1}])] + [outcome.(['P2_', arms_{2}])])', ...
        ([outcome.(['P1_', arms_{2}])] + [outcome.(['P2_', arms_{1}])])');
    
    outcomeTab.Carryover_P_Value = p(:);
    
    [~, p, ci] = ttest2( ...
        0.5*([outcome.(['P1_', arms_{1}])] - [outcome.(['P2_', arms_{2}])])', ...
        0.5*([outcome.(['P1_', arms_{2}])] - [outcome.(['P2_', arms_{1}])])');
    
    outcomeTab.([arms_{1}, '_', arms_{2}]) = 0.5 * (outcomeTab.(['Seq1_', arms_{1}, '_', arms_{2}]) - outcomeTab.(['Seq2_', arms_{2}, '_', arms_{1}]));
    outcomeTab.CI_005 = ci(1, :)';
    outcomeTab.CI_095 = ci(2, :)';
    outcomeTab.P_Value = p(:);
elseif strcmpi(obj.type, 'matched') %matched is crossover with ttest2
    k = 0;
    for p = 1:numel(obj.patients)
        if numel(intersect({obj.patients(p).data.name}, arms_)) ~= 2
            continue;
        end
        
        k = k + 1;
        out_ = obj.patients(p).(set_)(varargin{:}, 'interval', interval_, 'arms', arms_);
        
        for aIdx = 1:length(arms_)
            outcome(k).(arms_{aIdx}) = out_.(arms_{aIdx});
        end
    end
    
    outcomeTab = table;
    if length(fieldnames(outcome)) ~= 2 || length(outcome) == 1
        outcomeTab.NOT_ENOUGH_DATA = nan(size(out_.Properties.RowNames));
        outcomeTab.Properties.RowNames = out_.Properties.RowNames;
        return;
    end
    
    for aIdx = 1:length(arms_)
        data_ = [outcome.(arms_{aIdx})];
        outcomeTab.([arms_{aIdx}, '_Mean']) = nanmean(data_, 2);
        if var_
            outcomeTab.([arms_{aIdx}, '_SD']) = nanstd(data_, 0, 2);
        end
    end
    
    [~, p, ci] = ttest( ...
        [outcome.(arms_{1})]', ...
        [outcome.(arms_{2})]');
    
    if normality_
        outcomeTab.([arms_{1}, '_', arms_{2}, '_Normality']) = MAPUtils.swtest(([outcome.(arms_{1})] - [outcome.(arms_{2})])')';
    end
    
    outcomeTab.([arms_{1}, '_', arms_{2}, '_Diff']) = mean([outcome.(arms_{1})]-[outcome.(arms_{2})], 2);
    outcomeTab.([arms_{1}, '_', arms_{2}, '_CI_2_5']) = ci(1, :)';
    outcomeTab.([arms_{1}, '_', arms_{2}, '_CI_97_5']) = ci(2, :)';
    outcomeTab.P_Value = p(:);
else % assume that each patient has one outcome arm
    for p = 1:numel(obj.patients)
        aIdx = find(contains(arms_, obj.patients(p).arms));
        if isempty(baseline_)
            out_ = obj.patients(p).(set_)(varargin{:}, 'interval', interval_, 'arms', arms_{aIdx});
        else
            out_ = obj.patients(p).(set_)(varargin{:}, 'interval', interval_, 'arms', {baseline_{aIdx}, arms_{aIdx}});
            outcome(p).(baseline_{aIdx}) = out_.(baseline_{aIdx});
        end
        outcome(p).(arms_{aIdx}) = out_.(arms_{aIdx});
    end
    
    outcomeTab = table;
    if (~isempty(baseline_) && length(fieldnames(outcome)) ~= 4) || (isempty(baseline_) && length(fieldnames(outcome)) ~= 2)
        outcomeTab.NOT_ENOUGH_DATA = nan(size(out_.Properties.RowNames));
        outcomeTab.Properties.RowNames = out_.Properties.RowNames;
        return;
    end
    
    for aIdx = 1:length(arms_)
        if ~isempty(baseline_)
            data_ = [outcome.(baseline_{aIdx})];
            outcomeTab.([baseline_{aIdx}, '_Mean_n_', sprintf('%02d', size(data_, 2))]) = nanmean(data_, 2);
            if var_
                outcomeTab.([baseline_{aIdx}, '_SD']) = std(data_, 0, 2);
            end
        end
        data_ = [outcome.(arms_{aIdx})];
        outcomeTab.([arms_{aIdx}, '_Mean_n_', sprintf('%02d', size(data_, 2))]) = nanmean(data_, 2);
        if var_
            outcomeTab.([arms_{aIdx}, '_SD']) = nanstd(data_, 0, 2);
        end
        if isempty(baseline_)
            if normality_
                outcomeTab.([arms_{aIdx}, '_Normal']) = MAPUtils.swtest([outcome.(arms_{aIdx})]')';
            end
        else
            [~, p] = ttest( ...
                [outcome.(arms_{aIdx})]'-[outcome.(baseline_{aIdx})]');
            
            outcomeTab.([arms_{aIdx}, '_', baseline_{aIdx}]) = nanmean([outcome.(arms_{aIdx})]-[outcome.(baseline_{aIdx})], 2);
            if normality_
                outcomeTab.([arms_{aIdx}, '_', baseline_{aIdx}, '_Normal']) = MAPUtils.swtest([outcome.(arms_{aIdx})]'-[outcome.(baseline_{aIdx})]')';
            end
            outcomeTab.([arms_{aIdx}, '_', baseline_{aIdx}, '_P_Value']) = p(:);
        end
    end
    if ~isempty(baseline_)
        if normality_
            outcomeTab.([arms_{1}, '_', baseline_{1}, 'vs', arms_{2}, '_', baseline_{2}, '_Normality']) = (outcomeTab.([arms_{1}, '_', baseline_{1}, '_Normal']) > 0.05) & (outcomeTab.([arms_{2}, '_', baseline_{2}, '_Normal']) > 0.05);
        end
        
        [~, p, ci] = ttest2( ...
            [outcome.(arms_{1})]'-[outcome.(baseline_{1})]', ...
            [outcome.(arms_{2})]'-[outcome.(baseline_{2})]');
        
        outcomeTab.([arms_{1}, '_', baseline_{1}, 'vs', arms_{2}, '_', baseline_{2}]) = nanmean([outcome.(arms_{1})]-[outcome.(baseline_{1})], 2) ...
            -nanmean([outcome.(arms_{2})]-[outcome.(baseline_{2})], 2);
    else
        if normality_
            outcomeTab.([arms_{1}, 'vs', arms_{2}, '_Normality']) = (outcomeTab.([arms_{1}, '_Normal']) > 0.05) & (outcomeTab.([arms_{2}, '_Normal']) > 0.05);
        end
        
        [~, p, ci] = ttest2( ...
            [outcome.(arms_{1})]', ...
            [outcome.(arms_{2})]');
        
        outcomeTab.([arms_{1}, 'vs', arms_{2}]) = nanmean([outcome.(arms_{1})], 2) - nanmean([outcome.(arms_{2})], 2);
    end
    outcomeTab.CI_005 = ci(1, :)';
    outcomeTab.CI_095 = ci(2, :)';
    outcomeTab.P_Value = p(:);
end
outcomeTab.Properties.RowNames = out_.Properties.RowNames;
end
