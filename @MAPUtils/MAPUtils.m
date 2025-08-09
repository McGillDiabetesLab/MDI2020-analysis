classdef MAPUtils < matlab.mixin.Copyable
    %MAPUTILS Utils functions for MAP Data
    
    properties (GetAccess = public, SetAccess = public)
    end
    
    methods (Static)
        [ax, lgd] = plotData(ax, data, varargin)
        [ax, lgd] = plotSummary(ax, data, varargin)
        plotStats(fig, data)
        plotStatsSummary(fig, data, varargin)
        plotTable(ui, cellTable, varargin)
        
        %% merge
        function dMAP = merge(dMAP1, dMAP2)
            if nargin < 2
                if ~isempty(dMAP1)
                    dMAP = MAPData('struct', MAPUtils.structMerge(dMAP1.toStruct('format', 'array'), dMAP1(1).stepTime));
                else
                    dMAP = copy(dMAP1);
                end
            else
                if isempty(dMAP1)
                    dMAP = copy(dMAP2);
                elseif isempty(dMAP2)
                    dMAP = copy(dMAP1);
                else
                    dMAP = MAPData('struct', MAPUtils.structMerge([dMAP1.toStruct('format', 'array'), dMAP2.toStruct('format', 'array')], dMAP1.stepTime));
                    dMAP.name = sprintf('%s_%s', dMAP1.name, dMAP2.name);
                end
           end
        end
        
        %% structMerge
        function s = structMerge(ss, stepTime)
            if nargin < 2
                stepTime = 10;
            end
            
            if numel(ss) == 1
                s = ss(1);
                return;
            end
            
            if ~isfield(ss(1), 'timeStamp')
                if ~isfield(ss(1), 'time')
                    error('[MAPData][structMerge] Required timestamp/time field.');
                else
                    minTime = inf;
                    maxTime = -1;
                    for n = 1:numel(ss)
                        if isempty(ss(n).time)
                            continue;
                        end
                        minTime = min(minTime, min(ss(n).time));
                        maxTime = max(maxTime, max(ss(n).time));
                    end
                    
                    s.time = (minTime:stepTime:maxTime)';
                    
                    for fn = fieldnames(ss(n))'
                        if isstruct(s.(fn{1}))
                            continue;
                        end
                        if strcmp(fn{1}, 'time')
                            continue;
                        end
                        if any(arrayfun(@(a)(length(a.(fn{1})) ~= length(a.time)), ss))
                            s.(fn{1}) = ss(1).(fn{1});
                            continue;
                        end
                        s.(fn{1}) = nan(size(s.time));
                        for n = 1:numel(ss)
                            time = ss(n).time(~isnan(ss(n).(fn{1})));
                            if isempty(time)
                                continue;
                            end
                            val = ss(n).(fn{1})(~isnan(ss(n).(fn{1})));
                            idxR = any(abs(time-s.time') < stepTime/(24 * 60)/2, 2);
                            idxL = any(abs(s.time-time') < stepTime/(24 * 60)/2, 2);
                            s.(fn{1})(idxL) = val(idxR);
                        end
                    end
                end
            else
                startDateNum = inf;
                endDateNum = -1;
                for n = 1:numel(ss)
                    if isempty(ss(n).timeStamp)
                        continue;
                    end
                    startDateNum = min(startDateNum, min(ss(n).timeStamp));
                    endDateNum = max(endDateNum, max(ss(n).timeStamp));
                end
                s.timeStamp = (((startDateNum * 24 * 60):stepTime:(endDateNum * 24 * 60))') / (24 * 60);
                
                for fn = fieldnames(ss(1))'
                    if strcmp(fn{1}, 'timeStamp')
                        continue;
                    end
                    if any(arrayfun(@(a)isstruct(a.(fn{1})), ss))
                        continue;
                    end
                    if any(arrayfun(@(a)(length(a.(fn{1})) ~= length(a.timeStamp)), ss))
                        s.(fn{1}) = ss(1).(fn{1});
                        continue;
                    end
                    s.(fn{1}) = nan(size(s.timeStamp));
                    for n = 1:numel(ss)
                        if strcmp(fn{1}, 'timestampString')
                            continue
                        end
                        
                        timeStamp = ss(n).timeStamp(~isnan(ss(n).(fn{1})));
                        if isempty(timeStamp)
                            continue;
                        end
                        val = ss(n).(fn{1})(~isnan(ss(n).(fn{1})));
                        idxR = any(abs(timeStamp-s.timeStamp') < stepTime/(24 * 60)/2, 2);
                        idxL = any(abs(s.timeStamp-timeStamp') < stepTime/(24 * 60)/2, 2);
                        s.(fn{1})(idxL) = val(idxR);
                    end
                end
                
                if isfield(s, 'name')
                    s.name(strfind(s.name, '#'):end) = '';
                end
                
                s = MAPUtils.sortStruct(s);
            end
        end
        
        %% structResize
        function s = structResize(s, startDate, endDate, stepTime)
            if isempty(s) || isnat(startDate) || isnat(endDate)
                return;
            end
            
            if ~isfield(s, 'timeStamp') || isempty(s.timeStamp)
                return;
            end
            
            timeStamp = (((datenum(startDate) * 24 * 60):stepTime:(datenum(endDate) * 24 * 60))') / (24 * 60);
            
            for fn = fieldnames(s)'
                if strcmp(fn{1}, 'timeStamp')
                    continue;
                end
                
                if length(s.(fn{1})) ~= length(s.timeStamp)
                    continue;
                end
                
                temp = s.(fn{1});
                s.(fn{1}) = nan(size(timeStamp));
                idxR = any(abs(s.timeStamp-timeStamp') < stepTime/(24 * 60)/2, 2);
                idxL = any(abs(timeStamp-s.timeStamp') < stepTime/(24 * 60)/2, 2);
                s.(fn{1})(idxL) = temp(idxR);
            end
            
            s.timeStamp = timeStamp;
        end
        
        %% sortStruct
        function s = sortStruct(s)
            [sortedTimeStamp, iSort] = sort(s.timeStamp);
            [~, iUnique] = unique(round(sortedTimeStamp*60*24));
            n = length(s.timeStamp);
            
            for fn = fieldnames(s)'
                if length(s.(fn{1})) ~= length(s.timeStamp)
                    continue;
                end
                
                s.(fn{1}) = s.(fn{1})(iSort(any(iUnique == 1:n)));
            end
        end
        
        %% struct2array
        function v = struct2array(s, t, dt)
            if nargin < 3
                dt = 10;
            end
            time = s.time(:);
            value = s.value(:);
            
            if s.time(1) > 0
                time = [0; time];
                value = [value(end); value];
            end
            
            if s.time(end) < 1430
                time = [time; 1430];
                value = [value; value(end)];
            end
            
            v = interp1(time, value, (0:dt:(1440 - dt)), 'previous');
            v = v(floor(mod(t, 1440)/dt)+1);
            v = v(:);
        end
        
        %% concatSheets
        function shOut = concatSheets(sh1, sh2)
            [sA1, sA2] = size(sh1);
            [sB1, sB2] = size(sh2);
            shOut(1:sB1, sA2+1:sA2+sB2) = sh2;
            shOut(1:sA1, 1:sA2) = sh1;
        end
        
        %% roundToXmin
        function t = roundToXmin(t, x, rel)
            if isempty(t)
                return;
            end
            if ~isdatetime(t) && all(isnan(t))
                return;
            end
            if isdatetime(t) && all(isnat(t))
                return;
            end
            if nargin < 3
                rel = 1;
            end
            if rel
                if ~isdatetime(t)
                    tStart = t(find(~isnan(t), 1, 'first'));
                else
                    tStart = t(find(~isnat(t), 1, 'first'));
                end
                t = t - tStart;
            end
            if all(isdatetime(t)) && all(~isnat(t))
                t = datetime(datestr(x*round(datenum(t)*60*24/x)/60/24, 'dd-mmm-yyyy HH:MM:SS'));
            elseif isduration(t)
                t = minutes(x*round(minutes(t)*60*24/x)/60/24);
            elseif isnumeric(t)
                t = x * round(t*60*24/x) / 60 / 24;
            end
            if rel
                t = t + tStart;
            end
        end
        
        %% floorTo10min
        function t = floorToXmin(t, x, rel)
            if isempty(t)
                return;
            end
            if ~isdatetime(t) && all(isnan(t))
                return;
            end
            if isdatetime(t) && all(isnat(t))
                return;
            end
            if nargin < 3
                rel = 1;
            end
            if rel
                if ~isdatetime(t)
                    tStart = t(find(~isnan(t), 1, 'first'));
                else
                    tStart = t(find(~isnat(t), 1, 'first'));
                end
                t = t - tStart;
            end
            if isdatetime(t) && ~isnat(t)
                t = datetime(datestr(x*floor(datenum(t)*60*24/x)/60/24, 'dd-mmm-yyyy HH:MM:SS'));
            elseif isduration(t)
                t = minutes(x*floor(minutes(t)*60*24/x)/60/24);
            elseif isnumeric(t)
                t = x * floor(t*60*24/x) / 60 / 24;
            end
            if rel
                t = t + tStart;
            end
        end
        
        %% ceilTo10min
        function t = ceilToXmin(t, x, rel)
            if isempty(t)
                return;
            end
            if ~isdatetime(t) && all(isnan(t))
                return;
            end
            if isdatetime(t) && all(isnat(t))
                return;
            end
            if nargin < 3
                rel = 1;
            end
            if rel
                if ~isdatetime(t)
                    tStart = t(find(~isnan(t), 1, 'first'));
                else
                    tStart = t(find(~isnat(t), 1, 'first'));
                end
                t = t - tStart;
            end
            if isdatetime(t) && ~isnat(t)
                t = datetime(datestr(x*ceil(datenum(t)*60*24/x)/60/24, 'dd-mmm-yyyy HH:MM:SS'));
            elseif isduration(t)
                t = minutes(x*floor(minutes(t)*60*24/x)/60/24);
            elseif isnumeric(t)
                t = x * ceil(t*60*24/x) / 60 / 24;
            end
            if rel
                t = t + tStart;
            end
        end
        
        %% stof
        function f = toDouble(s)
            if isnumeric(s)
                f = s;
            elseif iscell(s)
                f = str2double(s);
            else
                f = NaN;
            end
        end
        
        %% smooth
        function gOut = smooth(gIn, action, glucoseInterpHolesSize, glucoseInterpHolesSizeAfterMealsOrBoluses, stepTime)
            gOut = fillmissing(gIn, 'pchip');
            
            missing = find(isnan(gIn));
            
            holes = double([1; (abs(diff(missing)) > 1)]);
            holeIdx = 0;
            for k = 1:length(holes)
                if holes(k) > 0
                    holeIdx = holeIdx + 1;
                end
                holes(k) = holeIdx;
            end
            if any(missing == 1)
                gOut(missing(holes == holes(missing == 1))) = NaN;
            end
            if any(missing == length(gIn))
                gOut(missing(holes == holes(missing == length(gIn)))) = NaN;
            end
            for idx = 1:max(holes)
                holeSize = sum(holes == idx) * stepTime;
                if holeSize < glucoseInterpHolesSize
                    continue;
                end
                if ~isempty(intersect(missing(holes == idx), find(action))) || ...
                        holeSize >= glucoseInterpHolesSizeAfterMealsOrBoluses
                    gOut(missing(holes == idx)) = NaN;
                end
            end
        end
        
        %% sensorTime
        function [perc, time] = sensorTime(glucose, stepTime)
            if nargin < 2
                stepTime = 10;
            end
            
            glucose = glucose(:);
            time = sum(~isnan(glucose)) * stepTime;
            perc = 100 * sum(~isnan(glucose)) / length(glucose);
        end
        
        %% timeIn
        function [perc, time, out] = timeIn(glucose, fromThresh, toThresh, stepTime)
            if nargin < 4
                stepTime = 10;
            end
            
            glucose = glucose(:);
            out = false(size(glucose));
            idxIn = find(glucose >= fromThresh & glucose < toThresh);
            out(idxIn) = true;
            time = length(idxIn) * stepTime;
            perc = 100 * length(idxIn) / sum(~isnan(glucose));
        end
        
        %% hypoCount
        function [count, out] = hypoCount(treats, glucose, thresh, duration, stepTime)
            if nargin < 3
                thresh = 3.9; %mmol/L
            end
            if nargin < 4
                duration = 20; %minutes
            end
            if nargin < 5
                stepTime = 10; %minutes
            end
            count = 0;
            out = false(size(glucose));
            
            if any(treats > 0)
                idxHypoZone = find(glucose < max(max(glucose(treats > 0)), thresh)+0.1);
            else
                idxHypoZone = find(glucose < max(5.0, thresh)+0.1);
            end
            
            if isempty(idxHypoZone)
                count = count + nansum(treats > 0);
                out = treats > 0;
            else
                hypoZones = double([1; (abs(diff(idxHypoZone)) > 1)]);
                hypoZoneIdx = 0;
                for k = 1:length(hypoZones)
                    if hypoZones(k) > 0
                        hypoZoneIdx = hypoZoneIdx + 1;
                    end
                    hypoZones(k) = hypoZoneIdx;
                end
                for idx = 1:max(hypoZones)
                    if nansum(treats(idxHypoZone(hypoZones == idx))) > 0
                        count = count + nansum(treats(idxHypoZone(hypoZones == idx)) > 0);
                        out(nonzeros(idxHypoZone(hypoZones == idx).*(treats(idxHypoZone(hypoZones == idx)) > 0))) = true;
                    else
                        hypoDuration = 0;
                        for k = idxHypoZone(hypoZones == idx)'
                            if glucose(k) < thresh
                                hypoDuration = hypoDuration + stepTime;
                            else
                                hypoDuration = 0;
                            end
                            if hypoDuration >= duration && ...
                                    hypoDuration < duration + stepTime% only count once
                                count = count + 1;
                                out(k) = true;
                            end
                        end
                    end
                end
            end
        end
        
        %% timeInCL
        function [percInCL, timeInCL] = timeInCL(closedLoopActive, stepTime)
            if nargin < 2
                stepTime = 10;
            end
            timeInCL = sum(closedLoopActive == 1) * stepTime;
            percInCL = 100 * sum(closedLoopActive == 1) / length(closedLoopActive);
        end
        
        %% AUC
        function auc = AUC(glucose, stepTime) % assumes a row array as input
            if nargin < 2
                stepTime = 10;
            end
            
            if length(glucose) == numel(glucose)
                glucose = glucose(:)';
            end
            
            glucose = glucose - glucose(:, 1);
            
            nanGlucose = isnan(glucose);
            row2nan = (sum(nanGlucose, 2) > 0.5 * size(glucose, 2));
            
            glucose = fillmissing(glucose', 'pchip')';
            glucose(glucose < 0) = 0;
            auc = nansum(glucose, 2) * stepTime / 60;
            
            auc(row2nan, :) = NaN;
        end
        
        %% Shapiro-Wilk and Shapiro-Francia normality tests.
        % Code adapted from https://www.mathworks.com/matlabcentral/fileexchange/13964-shapiro-wilk-and-shapiro-francia-normality-tests
        function p = swtest(X)
            if any(size(X) == 1) && ismatrix(X)
                X = X(:);
            end
            p = nan(1, size(X, 2));
            for k = 1:size(X, 2)
                x = X(~isnan(X(:, k)), k);
                if length(x) < 3
                    p(k) = NaN;
                    continue;
                end
                
                % First, calculate the a's for weights as a function of the m's
                % See Royston (1992, p. 117) and Royston (1993b, p. 38) for details
                % in the approximation.
                x = sort(x); % Sort the vector X in ascending order.
                n = length(x);
                mtilde = norminv(((1:n)' - 3 / 8)/(n + 1 / 4));
                weights = zeros(n, 1); % Preallocate the weights.
                if kurtosis(x) > 3
                    % The Shapiro-Francia test is better for leptokurtic samples.
                    weights = 1 / sqrt(mtilde'*mtilde) * mtilde;
                    %
                    % The Shapiro-Francia statistic W' is calculated to avoid excessive
                    % rounding errors for W' close to 1 (a potential problem in very
                    % large samples).
                    %
                    W = (weights' * x)^2 / ((x - mean(x))' * (x - mean(x)));
                    % Royston (1993a, p. 183):
                    nu = log(n);
                    u1 = log(nu) - nu;
                    u2 = log(nu) + 2 / nu;
                    mu = -1.2725 + (1.0521 * u1);
                    sigma = 1.0308 - (0.26758 * u2);
                    newSFstatistic = log(1-W);
                    %
                    % Compute the normalized Shapiro-Francia statistic and its p-value.
                    %
                    NormalSFstatistic = (newSFstatistic - mu) / sigma;
                    % Computes the p-value, Royston (1993a, p. 183).
                    p(k) = 1 - normcdf(NormalSFstatistic, 0, 1);
                else
                    % The Shapiro-Wilk test is better for platykurtic samples.
                    c = 1 / sqrt(mtilde'*mtilde) * mtilde;
                    u = 1 / sqrt(n);
                    % Royston (1992, p. 117) and Royston (1993b, p. 38):
                    PolyCoef_1 = [-2.706056, 4.434685, -2.071190, -0.147981, 0.221157, c(n)];
                    PolyCoef_2 = [-3.582633, 5.682633, -1.752461, -0.293762, 0.042981, c(n-1)];
                    % Royston (1992, p. 118) and Royston (1993b, p. 40, Table 1)
                    PolyCoef_3 = [-0.0006714, 0.0250540, -0.39978, 0.54400];
                    PolyCoef_4 = [-0.0020322, 0.0627670, -0.77857, 1.38220];
                    PolyCoef_5 = [0.00389150, -0.083751, -0.31082, -1.5861];
                    PolyCoef_6 = [0.00303020, -0.082676, -0.48030];
                    PolyCoef_7 = [0.459, -2.273];
                    weights(n) = polyval(PolyCoef_1, u);
                    weights(1) = -weights(n);
                    if n > 5
                        weights(n-1) = polyval(PolyCoef_2, u);
                        weights(2) = -weights(n-1);
                        
                        count = 3;
                        phi = (mtilde' * mtilde - 2 * mtilde(n)^2 - 2 * mtilde(n-1)^2) / ...
                            (1 - 2 * weights(n)^2 - 2 * weights(n-1)^2);
                    else
                        count = 2;
                        phi = (mtilde' * mtilde - 2 * mtilde(n)^2) / ...
                            (1 - 2 * weights(n)^2);
                    end
                    % Special attention when n = 3 (this is a special case).
                    if n == 3
                        % Royston (1992, p. 117)
                        weights(1) = 1 / sqrt(2);
                        weights(n) = -weights(1);
                        phi = 1;
                    end
                    %
                    % The vector 'WEIGHTS' obtained next corresponds to the same coefficients
                    % listed by Shapiro-Wilk in their original test for small samples.
                    %
                    weights(count : n-count+1) = mtilde(count : n-count+1) / sqrt(phi);
                    %
                    % The Shapiro-Wilk statistic W is calculated to avoid excessive rounding
                    % errors for W close to 1 (a potential problem in very large samples).
                    %
                    if x == mean(x)
                        p(k) = NaN;
                        continue;
                    end
                    
                    W = (weights' * x)^2 / ((x - mean(x))' * (x - mean(x)));
                    %
                    % Calculate the normalized W and its significance level (exact for
                    % n = 3). Royston (1992, p. 118) and Royston (1993b, p. 40, Table 1).
                    %
                    newn = log(n);
                    if (n >= 4) && (n <= 11)
                        
                        mu = polyval(PolyCoef_3, n);
                        sigma = exp(polyval(PolyCoef_4, n));
                        gam = polyval(PolyCoef_7, n);
                        
                        newSWstatistic = -log(gam-log(1-W));
                        
                    elseif n > 11
                        
                        mu = polyval(PolyCoef_5, newn);
                        sigma = exp(polyval(PolyCoef_6, newn));
                        
                        newSWstatistic = log(1-W);
                        
                    elseif n == 3
                        mu = 0;
                        sigma = 1;
                        newSWstatistic = 0;
                    end
                    %
                    % Compute the normalized Shapiro-Wilk statistic and its p-value.
                    %
                    NormalSWstatistic = (newSWstatistic - mu) / sigma;
                    
                    % NormalSWstatistic is referred to the upper tail of N(0,1),
                    % Royston (1992, p. 119).
                    p(k) = 1 - normcdf(NormalSWstatistic, 0, 1);
                    
                    % Special attention when n = 3 (this is a special case).
                    if n == 3
                        p(k) = 6 / pi * (asin(sqrt(W)) - asin(sqrt(3/4)));
                        % Royston (1982a, p. 121)
                    end
                end
            end
        end
        
        %% utest2 (just runs ranksum on matrix col vs col
        function [p, stat] = utest2(X1, X2, varargin)
            if any(size(X1) == 1) && ismatrix(X1)
                X1 = X1(:);
            end
            if any(size(X2) == 1) && ismatrix(X2)
                X2 = X2(:);
            end
            p = nan(1, size(X1, 2));
            for k = 1:size(X1, 2)
                x1 = X1(~isnan(X1(:, k)), k);
                x2 = X2(~isnan(X2(:, k)), k);
                if isempty(x1) || isempty(x2)
                    p(k) = NaN;
                    stat.median(1, k) = NaN;
                    stat.ci(:, k) = [NaN; NaN];
                    continue;
                end
                if min(length(x1), length(x2)) < 20 && length(x1) + length(x2) < 40
                    p(k) = ranksum(x1, x2, 'method', 'exact', varargin{:});
                else
                    p(k) = ranksum(x1, x2, 'method', 'approximate', varargin{:});
                end
                
                n = length(x1);
                m = length(x2);
                diff = x1' - x2;
                [diffSorted, idx] = sort(diff(:));
                stat.median(:, k) = median(diffSorted);
                R(idx) = 1:numel(diff);
                
                R = reshape(R, [n, m]);
                C = fix(m*n/2-norminv(1-0.05/2)*realsqrt(m*n*(m + n + 1)/12));
                if C > 0
                    stat.ci(:, k) = [diff(R == C); diff(R == n*m-C+1)];
                else
                    stat.ci(:, k) = [NaN; NaN];
                end
                
                clear R;
            end
        end
        
        %% utest
        function [p, stat] = utest(X1, X2)
            if nargin < 2
                X2 = zeros(size(X1));
            end
            if any(size(X1) == 1) && ismatrix(X1)
                X1 = X1(:);
            end
            if any(size(X2) == 1) && ismatrix(X2)
                X2 = X2(:);
            end
            p = nan(1, size(X1, 2));
            for k = 1:size(X1, 2)
                x1 = X1(~isnan(X1(:, k)) & ~isnan(X2(:, k)), k);
                x2 = X2(~isnan(X1(:, k)) & ~isnan(X2(:, k)), k);
                if isempty(x1) || isempty(x2)
                    p(k) = NaN;
                    stat.median(1, k) = NaN;
                    stat.ci(:, k) = [NaN; NaN];
                    continue;
                end
                if min(length(x1), length(x2)) < 30 && length(x1) + length(x2) < 60
                    p(k) = signrank(x1, x2, 'method', 'exact');
                else
                    p(k) = signrank(x1, x2, 'method', 'approximate');
                end
                
                dff = sort(x1(:)'-x2(:)'); %difference between x1 and x2
%                 dff(dff == 0) = []; %eliminate null variations
                n = length(dff); %number of ranks
                if n == 0 %if all variations are null variations exit function
                    stat.median(1, k) = NaN;
                    stat.ci(:, k) = [NaN; NaN];
                    continue;
                end
                
                %Ranks of absolute value of samples differences with sign
                [I, J] = ndgrid(dff, dff);
                d = triu(I+J) ./ 2; %Walsh averages triangular matrix
                d(tril(ones(size(d)), -1) == 1) = NaN;
                ld = sort(d(~isnan(d))); %linearization of Walsh averages matrix
%                 ld = sort(d(d ~= 0)); %linearization of Walsh averages matrix
                stat.median(1, k) = median(ld); %Hodges-Lehmann estimator
                if n > 15
                    A = n * (n + 1) / 4;
                    B = realsqrt(n*(n + 1)*(2 * n + 1)/24);
                    Za = -realsqrt(2) .* erfcinv(2.*(1 - 0.05 / 2));
                    T = fix(A-Za.*B);
                else
                    TC = [0, 0, 0, 0, 0, 0, 2, 3, 5, 8, 10, 13, 17, 21, 25];
                    T = TC(n);
                end
                stat.ci(:, k) = [ld(T+1); ld(end-T)];
                stat.res = ld;
                
                clear I J ld
            end
        end
        
        function bolus = totalDailyBolus(data, indices)
            if nargin < 2
                indices = 1:length(data.bolusInsulin);
            end
            if ~isempty(data.bolusInsulin) && ~all(isnan(data.bolusInsulin))
                bolus = nansum(data.bolusInsulin(indices));
            else
                bolus = NaN;
            end
        end
        
        function basal = totalDailyBasal(data, indices)
            if nargin < 2
                indices = 1:length(data.basalInsulin);
            end
            basal = 0;
            if ~isempty(data.basalInsulin)
                basal = basal + nansum(data.basalInsulin(indices)) * data.stepTime / 60;
            end
            if ~isempty(data.basalInjection)
                basal = basal + nansum(data.basalInjection(indices));
            end
        end
        
        function tdd = totalDailyDose(data, indices)
            if nargin < 2
                indices = 1:length(data.basalInsulin);
            end
            tdd = MAPUtils.totalDailyBasal(data, indices) + MAPUtils.totalDailyBolus(data, indices);
        end
        
        function b = xlscol(a)
            base = 26;
            if iscell(a)
                b = cellfun(@xlscol, a, 'UniformOutput', false); % handles mixed case too
            elseif ischar(a)
                if contains(a, ':') % i.e. if is a range
                    b = cellfun(@xlscol, regexp(a, ':', 'split'));
                else % if isempty(strfind(a, ':')) % i.e. if not a range
                    b = a(isletter(a)); % get rid of numbers and symbols
                    if isempty(b)
                        b = {[]};
                    else % if ~isempty(a);
                        b = double(upper(b)) - 64; % convert ASCII to number from 1 to 26
                        n = length(b); % number of characters
                        b = b * base.^((n - 1):-1:0)';
                    end % if isempty(a)
                end % if ~isempty(strfind(a, ':')) % i.e. if is a range
            elseif isnumeric(a) && numel(a) ~= 1
                b = arrayfun(@xlscol, a, 'UniformOutput', false);
            else % if isnumeric(a) && numel(a) == 1
                n = ceil(log(a)/log(base)); % estimate number of digits
                d = cumsum(base.^(0:n + 1)); % offset
                n = find(a >= d, 1, 'last'); % actual number of digits
                d = d(n:-1:1); % reverse and shorten
                r = mod(floor((a - d)./base.^(n - 1:-1:0)), base) + 1; % modulus
                b = char(r+64); % convert number to ASCII
            end % if iscell(a)
            % attempt to "flatten" cell, by removing nesting
            if iscell(b) && (iscell([b{:}]) || isnumeric([b{:}]))
                b = [b{:}];
            end % if iscell(b) && (iscell([b{:}]) || isnumeric([ba{:}]))
        end
        
        function customizeExcel(filename, sheet)
            % Open file
            objExcel = actxserver('Excel.Application');
            objExcel.DisplayAlerts = false;
            objExcel.Workbooks.Open(fullfile(pwd, filename)); % Full path is necessary!
            
            % Check if empty 'Sheet' exist
            [~, sheets] = xlsfinfo(filename);
            sheetExist = any(contains(sheets, 'Sheet'));
            if sheetExist
                % Delete sheets.
                try
                    % Throws an error if the sheets do not exist.
                    objExcel.ActiveWorkbook.Worksheets.Item(['Sheet', '1']).Delete;
                    objExcel.ActiveWorkbook.Worksheets.Item(['Sheet', '2']).Delete;
                    objExcel.ActiveWorkbook.Worksheets.Item(['Sheet', '3']).Delete;
                catch
                    % Do nothing.
                end
                % Save
                objExcel.ActiveWorkbook.Save;
            end
            
            % Post processing: Merge, Bold, Italic
            eSheets = objExcel.ActiveWorkbook.Worksheets;
            for sheetname = fieldnames(sheet)'
                eSheet = eSheets.get('Item', sheetname{1});
                eSheet.Activate;
                
                if contains(sheetname{1}, 'outcome', 'IgnoreCase', true)
                    % Freeze column B3
                    objExcel.ActiveWindow.SplitRow = 2;
                    objExcel.ActiveWindow.SplitColumn = 1;
                    objExcel.ActiveWindow.FreezePanes = 1;
                end
                
                if contains(sheetname{1}, 'summary', 'IgnoreCase', true)
                    % Freeze column B9
                    objExcel.ActiveWindow.SplitRow = 8;
                    objExcel.ActiveWindow.SplitColumn = 1;
                    objExcel.ActiveWindow.FreezePanes = 1;
                end
                
                if contains(sheetname{1}, 'compare', 'IgnoreCase', true)
                    % Freeze column B2
                    objExcel.ActiveWindow.SplitRow = 1;
                    objExcel.ActiveWindow.SplitColumn = 1;
                    objExcel.ActiveWindow.FreezePanes = 1;
                end
                
                if contains(sheetname{1}, 'demographics', 'IgnoreCase', true)
                    % Freeze column B4
                    objExcel.ActiveWindow.SplitRow = 3;
                    objExcel.ActiveWindow.SplitColumn = 1;
                    objExcel.ActiveWindow.FreezePanes = 1;
                end
                
                if contains(sheetname{1}, 'lme', 'IgnoreCase', true)
                    % Freeze column B4
                    objExcel.ActiveWindow.SplitRow = 1;
                    objExcel.ActiveWindow.SplitColumn = 4;
                    objExcel.ActiveWindow.FreezePanes = 1;
                end
                
                for n = 1:size(sheet.(sheetname{1}), 1)
                    for m = 1:size(sheet.(sheetname{1}), 2)
                        if isempty(sheet.(sheetname{1}){n, m})
                            continue;
                        end
                        
                        if ~ischar(sheet.(sheetname{1}){n, m})
                            continue;
                        end
                        
                        if ~contains(sheetname{1}, 'compare', 'IgnoreCase', true)
                            if ~isempty(strfind(sheet.(sheetname{1}){n, m}, '=') ) && (strfind(sheet.(sheetname{1}){n, m}, '=') == 1)
                                cells = eSheet.get('Range', [MAPUtils.xlscol(m), num2str(n), ':', MAPUtils.xlscol(m), num2str(n)]);
                                set(cells.Font, 'Bold', true);
                            end
                        end
                        
                        if strfind(sheet.(sheetname{1}){n, m}, '##merge') == 1
                            mergeIdx = textscan(sheet.(sheetname{1}){n, m}, '##merge{%d}{%d}', 1);
                            if ~isempty(mergeIdx) && numel(mergeIdx) == 2
                                colOffset = double(mergeIdx{1});
                                rowOffset = double(mergeIdx{2});
                                
                                cells = eSheet.get('Range', [MAPUtils.xlscol(m+colOffset), num2str(n+rowOffset), ':', MAPUtils.xlscol(m), num2str(n)]);
                                cells.MergeCells = true;
                                
                                % this is a hack, TODO, apply these with
                                % some kind of flag
                                cells.ColumnWidth = 10.0;
                                cells.HorizontalAlignment = -4108;
                                set(cells.Font, 'Bold', true);
                            end
                        end
                    end
                end
            end
            
            % Clean up
            objExcel.ActiveWorkbook.Save;
            objExcel.ActiveWorkbook.Close;
            objExcel.Quit;
            objExcel.delete;
        end
    end
    
    methods (Access = public)
        
        %% Empty constructor
        function obj = MAPUtils()
        end
        
        %% Data manipulation functions
        function [perc, time] = getSensorTime(obj, indices, useGlucoseInterp)
            if numel(obj) > 1
                if nargin < 2
                    perc = arrayfun(@(c)(c.getSensorTime()), obj);
                elseif nargin < 3
                    perc = arrayfun(@(c)(c.getSensorTime(indices)), obj);
                else
                    perc = arrayfun(@(c)(c.getSensorTime(indices, useGlucoseInterp)), obj);
                end
                return;
            end
            if nargin < 3
                useGlucoseInterp = obj.useGlucoseInterp;
            end
            if useGlucoseInterp
                glucose = obj.glucoseInterp;
            else
                glucose = obj.sensorGlucose;
            end
            if isempty(glucose)
                perc = NaN;
                time = NaN;
                return;
            end
            if nargin < 2
                indices = 1:length(glucose);
            end
            if isempty(indices)
                indices = 1:length(glucose);
            end
            [perc, time] = MAPUtils.sensorTime(glucose(indices), obj.stepTime);
        end
        
        function [perc, time, out] = getTimeIn(obj, fromThresh, toThresh, indices, useGlucoseInterp)
            if numel(obj) > 1
                if nargin < 4
                    perc = arrayfun(@(c)(c.getTimeIn(fromThresh, toThresh)), obj);
                elseif nargin < 5
                    if ~ischar(indices) && numel(indices) == numel(obj)
                        perc = arrayfun(@(c, d)(c.getTimeIn(fromThresh, toThresh, d{1})), obj, indices);
                    else
                        perc = arrayfun(@(c)(c.getTimeIn(fromThresh, toThresh, indices)), obj);
                    end
                else
                    if numel(indices) == numel(obj)
                        perc = arrayfun(@(c, d)(c.getTimeIn(fromThresh, toThresh, d{1}, useGlucoseInterp)), obj, indices);
                    else
                        perc = arrayfun(@(c)(c.getTimeIn(fromThresh, toThresh, indices, useGlucoseInterp)), obj);
                    end
                end
                return;
            elseif isempty(obj)
                perc = NaN;
                time = NaN;
                out = NaN;
                return;
            end
            if nargin < 4
                indices = 1:length(obj.time);
            else
                if ischar(indices)
                    switch indices
                        case {'day', 'daytime'}
                            indices = obj.intervalToIndices(obj.dayInterval);
                        case {'night', 'nighttime'}
                            indices = obj.intervalToIndices(obj.nightInterval);
                        case {'all', 'anytime'}
                            indices = 1:length(obj.time);
                        otherwise
                            error('MAPUtis:getTimeIn unhandled option for indices %s', indices)
                    end
                else
                    indices = intersect(indices, 1:length(obj.time));
                end
            end
            if nargin < 5
                useGlucoseInterp = obj.useGlucoseInterp;
            end
            if useGlucoseInterp
                glucose = obj.glucoseInterp;
            else
                glucose = obj.sensorGlucose;
            end
            if isempty(glucose)
                perc = NaN;
                time = NaN;
                return;
            end
            [perc, time, out] = MAPUtils.timeIn(glucose(indices), fromThresh, toThresh, obj.stepTime);
        end
        
        function [count, out] = getHypoCount(obj, indices, thresh, duration)
            if numel(obj) > 1
                if nargin < 2
                    count = arrayfun(@(c)(c.getHypoCount()), obj);
                elseif nargin < 3
                    count = arrayfun(@(c)(c.getHypoCount(indices)), obj);
                elseif nargin < 4
                    count = arrayfun(@(c)(c.getHypoCount(indices, thresh)), obj);
                else
                    count = arrayfun(@(c)(c.getHypoCount(indices, thresh, duration)), obj);
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.time);
            else
                indices = intersect(indices, 1:length(obj.time));
            end
            if nargin < 3
                thresh = obj.thresholdHypoglycemiaEvent; %mmol/L
            end
            if nargin < 4
                duration = obj.durationHypoglycemiaEvent; %minutes
            end
            if obj.useGlucoseInterp
                glucose = obj.glucoseInterp;
            else
                glucose = obj.sensorGlucose;
            end
            if nansum(obj.bloodGlucose) > 0 && nansum(obj.sensorGlucose) > 0
                minGlucose = min(obj.bloodGlucose, obj.sensorGlucose);
            elseif nansum(obj.bloodGlucose) > 0 && nansum(obj.sensorGlucose) == 0
                minGlucose = obj.bloodGlucose;
            elseif nansum(obj.bloodGlucose) == 0 && nansum(obj.sensorGlucose) > 0
                minGlucose = obj.sensorGlucose;
            else
                minGlucose = nan(size(glucose));
            end
            if isempty(glucose)
                glucose = minGlucose;
            else
                glucose = min(glucose, minGlucose);
            end
            
            if nansum(obj.treats) > 0
                treats = obj.treats;
            else
                treats = zeros(size(glucose));
            end
            [count, out] = MAPUtils.hypoCount(treats(indices), glucose(indices), thresh, duration, obj.stepTime);
        end
        
        function indices_ = intervalToIndices(obj, interval_)
            if numel(obj) > 1
                if nargin < 2
                    indices_ = arrayfun(@(c)(c.intervalToIndices()), obj, 'UniformOutput', false);
                else
                    indices_ = arrayfun(@(c)(c.intervalToIndices(interval_)), obj, 'UniformOutput', false);
                end
                return;
            end
            
            if nargin < 2
                interval_ = [0, 0];
            end
            
            time_ = obj.time + obj.startTime;
            interval_ = mod(interval_, 1440);
            if interval_(2) <= interval_(1)
                indices_ = find(mod(time_, 1440) >= interval_(1) | mod(time_, 1440) < interval_(2));
            else
                indices_ = find(mod(time_, 1440) >= interval_(1) & mod(time_, 1440) < interval_(2));
            end
            
            if ~isnan(obj.outcomeStartTime)
                indices_(obj.time(indices_) < obj.outcomeStartTime-obj.startTime) = [];
            end
            if ~isnan(obj.outcomeEndTime)
                indices_(obj.time(indices_) >= obj.outcomeEndTime-obj.startTime) = [];
            end
            
            indices_ = indices_(:);
        end
        
        function auc = getAUC(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    auc = arrayfun(@(c)(c.getAUC()), obj);
                else
                    auc = arrayfun(@(c)(c.getAUC(indices)), obj);
                end
                return;
            end
            
            if nargin < 2
                indices = 1:length(obj.time);
            else
                indices = intersect(indices, 1:length(obj.time));
            end
            
            if obj.useGlucoseInterp
                glucose = obj.glucoseInterp(indices);
            else
                glucose = obj.sensorGlucose(indices);
            end
            
            auc = MAPUtils.AUC(glucose, obj.stepTime);
        end
        
        function hba1c = getHbA1c(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    hba1c = arrayfun(@(c)(c.getHbA1c()), obj);
                else
                    hba1c = arrayfun(@(c)(c.getHbA1c(indices)), obj);
                end
                return;
            end
            
            if nargin < 2
                indices = 1:length(obj.sensorGlucose);
            else
                indices = intersect(indices, 1:length(obj.sensorGlucose));
            end
            
            if obj.useGlucoseInterp
                glucose = obj.glucoseInterp(indices);
            else
                glucose = obj.sensorGlucose(indices);
            end
            
            meanGlucose = nanmean(glucose);
            hba1c = (meanGlucose * 18.018 + 46.7) / 28.7;
        end
        
        function riskBG = getRiskBG(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    riskBG = arrayfun(@(c)(c.getRiskBG()), obj, 'UniformOutput', false);
                else
                    riskBG = arrayfun(@(c)(c.getRiskBG(indices)), obj, 'UniformOutput', false);
                end
                return;
            end
            
            if nargin < 2
                indices = 1:length(obj.time);
            else
                indices = intersect(indices, 1:length(obj.time));
            end
            
            alpha = 1.026;
            beta = 1.861;
            gamma = 1.794;
            
            if obj.useGlucoseInterp
                glucose = obj.glucoseInterp(indices);
            else
                glucose = obj.sensorGlucose(indices);
            end
            
            riskBG = gamma * (log(glucose).^alpha - beta);
        end
        
        function LBGI = getLBGI(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    LBGI = arrayfun(@(c)(c.getLBGI()), obj);
                else
                    LBGI = arrayfun(@(c)(c.getLBGI(indices)), obj);
                end
                return;
            end
            
            if nargin < 2
                riskBG = obj.getRiskBG();
            else
                riskBG = obj.getRiskBG(indices);
            end
            
            if sum(riskBG < 0) == 0
                LBGI = 0;
            else
                LBGI = -nansum(riskBG(riskBG < 0)) / sum(riskBG < 0);
            end
        end
        
        function HBGI = getHBGI(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    HBGI = arrayfun(@(c)(c.getHBGI()), obj);
                else
                    HBGI = arrayfun(@(c)(c.getHBGI(indices)), obj);
                end
                return;
            end
            
            if nargin < 2
                riskBG = obj.getRiskBG();
            else
                riskBG = obj.getRiskBG(indices);
            end
            
            if sum(riskBG > 0) == 0
                HBGI = 0;
            else
                HBGI = nansum(riskBG(riskBG > 0)) / sum(riskBG > 0);
            end
        end
        
        function [treats, nbrOfTreats] = getTreats(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    [treats, nbrOfTreats] = arrayfun(@(c)(c.getTreats()), obj);
                else
                    [treats, nbrOfTreats] = arrayfun(@(c)(c.getTreats(indices)), obj);
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.treats);
            else
                indices = intersect(indices, 1:length(obj.treats));
            end
            
            if ~isempty(obj.treats) && ~all(isnan(obj.treats))
                treats = nansum(obj.treats(indices));
                nbrOfTreats = nansum(obj.treats(indices) > 0);
            else
                treats = NaN;
                nbrOfTreats = NaN;
            end
        end
        
        function carbs = getCarbs(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    carbs = arrayfun(@(c)(c.getCarbs()), obj);
                else
                    carbs = arrayfun(@(c)(c.getCarbs(indices)), obj);
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.carbs);
            else
                indices = intersect(indices, 1:length(obj.carbs));
            end
            
            if obj.useActualCarbs && nansum(obj.carbsActual) > 0
                carbs = nansum(obj.carbsActual(indices));
            elseif nansum(obj.carbs) > 0
                carbs = nansum(obj.carbs(indices));
            else
                carbs = NaN;
            end
        end
        
        function bolus = getTDCorrBolus(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    bolus = arrayfun(@(c)(c.getTDBolus()), obj);
                else
                    bolus = arrayfun(@(c)(c.getTDBolus(indices)), obj);
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.bolusInsulin);
            else
                indices = intersect(indices, 1:length(obj.bolusInsulin));
            end
            
            if ~isempty(obj.bolusInsulin) && ~all(isnan(obj.bolusInsulin))
                bolusInsulin = obj.bolusInsulin(indices);
                carbs = obj.carbs(indices);
                bolusInsulin(carbs > 0) = 0;
                bolus = nansum(bolusInsulin);
            else
                bolus = NaN;
            end
        end
        
        function [bolus, bolusPerDay] = getTDBolus(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    bolus = arrayfun(@(c)(c.getTDBolus()), obj);
                else
                    if numel(indices) == numel(obj)
                        bolus = arrayfun(@(c, d)(c.getTDBolus(d{1})), obj, indices);
                    else
                        bolus = arrayfun(@(c)(c.getTDBolus(indices)), obj);
                    end
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.bolusInsulin);
            else
                indices = intersect(indices, 1:length(obj.bolusInsulin));
            end
            
            if ~isempty(obj.bolusInsulin) && ~all(isnan(obj.bolusInsulin))
                bolus = nansum(obj.bolusInsulin(indices));
            else
                bolus = NaN;
            end
            
            bolusPerDay = 1440 * bolus / (length(indices) * obj.stepTime);
        end
        
        function [basal, basalPerDay] = getTDBasal(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    basal = arrayfun(@(c)(c.getTDBasal()), obj);
                else
                    if numel(indices) == numel(obj)
                        basal = arrayfun(@(c, d)(c.getTDBasal(d{1})), obj, indices);
                    else
                        basal = arrayfun(@(c)(c.getTDBasal(indices)), obj);
                    end
                end
                return;
            end
            basal = 0;
            basalPerDay = 0;
            basalSource = '';
            if nansum(obj.basalInsulin) > 0 && nansum(obj.basalInjection) > 0
                basalSource = 'basalInjection';
                warning('[MAPUtils][getTDBasal] Basal data from basalInsulin and basalInjection ... will default to basalInsulin.');
            elseif nansum(obj.basalInsulin) == 0 && nansum(obj.basalInjection) > 0
                basalSource = 'basalInjection';
            elseif nansum(obj.basalInsulin) > 0 && nansum(obj.basalInjection) == 0
                basalSource = 'basalInsulin';
            end
            
            switch basalSource
                case 'basalInsulin'
                    if nargin < 2
                        indices = 1:length(obj.basalInsulin);
                    else
                        indices = intersect(indices, 1:length(obj.basalInsulin));
                    end
                    if ~isempty(obj.basalInsulin) && ~all(isnan(obj.basalInsulin))
                        basal = basal + nansum(obj.basalInsulin(indices)) * obj.stepTime / 60;
                        basalPerDay = 1440 * basal / (length(indices) * obj.stepTime);
                    else
                        basal = NaN;
                    end
                case 'basalInjection'
                    if nargin < 2
                        indices = 1:length(obj.basalInjection);
                    else
                        indices = intersect(indices, 1:length(obj.basalInjection));
                    end
                    if ~isempty(obj.basalInjection) && ~all(isnan(obj.basalInjection))
                        basal = basal + nansum(obj.basalInjection(indices));
                        basalPerDay = basal / round((length(indices) * obj.stepTime) / 1440);
                    else
                        basal = NaN;
                    end
            end            
        end
        
        function [tdd, tddPerDay] = getTDD(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    tdd = arrayfun(@(c)(c.getTDD()), obj);
                else
                    if numel(indices) == numel(obj)
                        tdd = arrayfun(@(c, d)(c.getTDD(d{1})), obj, indices);
                    else
                        tdd = arrayfun(@(c)(c.getTDD(indices)), obj);
                    end
                end
                return;
            elseif isempty(obj)
                tdd = NaN;
                tddPerDay = NaN;
                return;
            end
            if nargin < 2
                indices = 1:max(length(obj.basalInsulin), length(obj.basalInjection));
            else
                indices = intersect(indices, 1:max(length(obj.basalInsulin), length(obj.basalInjection)));
            end
            
            tdd = obj.getTDBasal(indices) + obj.getTDBolus(indices);
            tddPerDay = 1440 * tdd / (length(indices) * obj.stepTime);
        end
        
        function meanGlucose = getMeanGlucose(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    meanGlucose = arrayfun(@(c)(c.getMeanGlucose()), obj);
                else
                    meanGlucose = arrayfun(@(c)(c.getMeanGlucose(indices)), obj);
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.glucoseInterp);
            else
                indices = intersect(indices, 1:length(obj.glucoseInterp));
            end
            
            meanGlucose = nanmean(obj.glucoseInterp(indices));
        end
        
        function sdGlucose = getSDGlucose(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    sdGlucose = arrayfun(@(c)(c.getSDGlucose()), obj);
                else
                    sdGlucose = arrayfun(@(c)(c.getSDGlucose(indices)), obj);
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.glucoseInterp);
            else
                indices = intersect(indices, 1:length(obj.glucoseInterp));
            end
            
            sdGlucose = nanstd(obj.glucoseInterp(indices));
        end
        
        function fpg = getFastingGlucose(obj)
            if numel(obj) > 1
                fpg = arrayfun(@(c)(c.getFastingGlucose()), obj);
                return;
            end
            
            computeIOB_ = obj.computeIOB;
            obj.computeIOB = 1;
            indices = intersect(find(obj.iobBolus < 0.05), obj.intervalToIndices(obj.nightInterval));
            
            fpg = nanmean(obj.glucoseInterp(indices));
            obj.computeIOB = computeIOB_;
        end
        
        function [percInCL, timeInCL] = getTimeInCL(obj, indices)
            if numel(obj) > 1
                if nargin < 2
                    percInCL = arrayfun(@(c)(c.getTimeInCL()), obj);
                else
                    percInCL = arrayfun(@(c)(c.getTimeInCL(indices)), obj);
                end
                return;
            end
            if nargin < 2
                indices = 1:length(obj.closedLoopActive);
            else
                indices = intersect(indices, 1:length(obj.closedLoopActive));
            end
            [percInCL, timeInCL] = MAPUtils.timeInCL(obj.closedLoopActive(indices), obj.stepTime);
        end
        
        function out = getMealCarbFactor(obj, meal)
            if numel(obj) > 1
                out = arrayfun(@(c)(c.getMealCarbFactor(meal)), obj);
                return;
            end
            time_ = obj.time + obj.startTime;
            indices = find(mod(time_, 1440) >= obj.([meal, 'Interval'])(1) & mod(time_, 1440) <= obj.([meal, 'Interval'])(2));
            idxMealAll = indices(obj.carbs(indices) > 0);
            idxMealAll = idxMealAll(obj.bolusInsulin(idxMealAll) > 0);
            
            if ~isempty(idxMealAll)
                separateMealsIdx = [(flipud(diff(flipud(idxMealAll))) < -6); true];
                idxMeal = idxMealAll(separateMealsIdx);
                
                out = nanmean(obj.carbFactors(idxMeal));
            else
                out = mode(obj.carbFactors(indices));
            end
        end
        
        function out = getNightBasal(obj)
            time_ = obj.time + 60 * hour(obj.startDate) + minute(obj.startDate);
            indices = mod(time_, 1440) >= obj.nightInterval(1) | mod(time_, 1440) < obj.nightInterval(2);
            
            out = nansum(obj.pumpBasals(indices)) + nansum(obj.basalInjection);
        end
        
        %% Functions to load data
        function fromMAT(obj, filename)
            structname = replace(class(obj), 'MAP', 't');
            load(filename, structname);
            
            eval(['obj.fromStruct(', structname, ');']);
        end
        
        %% Functions to save data
        function export(obj, varargin)
            if mod(numel(varargin), 2) == 0
                obj.toFigure(varargin{:});
            else
                filename = varargin{1};
                [~, ~, ext] = fileparts(filename);
                switch ext
                    case '.mat'
                        obj.toMAT(varargin{:});
                    case '.csv'
                        obj.toCSV(varargin{:});
                    case '.xlsx'
                        obj.toExcel(varargin{:});
                    case '.png'
                        obj.toPNG(varargin{:});
                    otherwise
                        error('[MAPUtilis][export] Unkown format %s', ext)
                end
            end
        end
        
        function toMAT(obj, path, varargin)
            [path_, name_] = obj.ensureFolder(path);
            path = [path_, filesep, name_, '.mat'];
            
            eval([replace(class(obj), 'MAP', 't'), ' = obj.toStruct();']);
            save(path, replace(class(obj), 'MAP', 't'), varargin{:});
        end
        
        function toCSV(obj, filename, varargin)
            if nargin < 2
                filename = obj.name;
            end
            
            sheet = obj.toSheet('individual', true, varargin{:});
            if isstruct(sheet) % multiple sheets
                if exist(filename, 'dir') ~= 7
                    mkdir(filename);
                end
                
                for fn = fieldnames(sheet)'
                    writetable(cell2table(sheet.(fn{1})), [filename, filesep, fn{1}, '.csv'], 'WriteVariableNames', false, 'Delimiter', ';', 'QuoteStrings', true);
                end
            else % One sheet
                writetable(cell2table(sheet), [filename, '.csv'], 'WriteVariableNames', false, 'Delimiter', ';', 'QuoteStrings', true);
            end
        end
        
        function toExcel(obj, filename, varargin)
            if nargin < 2
                filename = obj.name;
            end
            
            [path_, name_] = obj.ensureFolder(filename);
            filename = [path_, filesep, name_, '.xlsx'];
            
            sheet = obj.toSheet(varargin{:}, 'qqplot', [path_, filesep, 'qqplot_', name_, filesep]);
            if isstruct(sheet) % multiple sheets
                for sheetname = fieldnames(sheet)'
                    if ~isempty(sheet.(sheetname{1}))
                        sheetName = sheetname{1};
                        sheetName(32:end) = [];
                        writetable(cell2table(sheet.(sheetname{1})), filename, 'WriteVariableNames', false, 'Sheet', sheetName, 'UseExcel', true);
                    end
                end
            else % One sheet
                writetable(cell2table(sheet), filename, 'WriteVariableNames', false)
            end
            
            MAPUtils.customizeExcel(filename, sheet);
        end
    end
    
    methods (Access = protected)
        function [path_, name_] = ensureFolder(obj, filename)
            if filename(end) == '/' || filename(end) == '\'
                path_ = filename;
                name_ = obj.name;
            else
                [path_, name_] = fileparts(filename);
            end
            
            if ~isempty(path_)
                if exist(path_, 'dir') ~= 7
                    mkdir(path_);
                end
            else
                path_ = '.';
            end
        end
    end
end
