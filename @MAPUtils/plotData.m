function [ax, lgd] = plotData(ax, data, varargin)

% if no ax just use the current one
if mod(nargin, 2) == 1
    if nargin > 1
        varargin(2:end+1) = varargin;
        varargin{1} = data;
    end
    data = ax;
    ax = gca;
end

% Default options
title_ = '';
if isfield(data, 'stepTime') || (isa(data, 'MAPData') && isprop(data, 'stepTime'))
    stepTime_ = data.stepTime;
else
    stepTime_ = 10;
end
if isfield(data, 'startTime')
    startTime_ = data.startTime;
else
    idxTimeStampNonNan = find(~isnan(data.timeStamp), 1);
    startTime_ = mod(data.timeStamp(idxTimeStampNonNan), 1) * 24 * 60 - (idxTimeStampNonNan - 1) * stepTime_;
end
duration_ = length(data.time) * stepTime_;
corrBolus_ = nansum(data.carbFactors) > 0;
plotSecondaryInfusions_ = false;

% Get options
for nVar = 1:2:length(varargin)
    switch lower(varargin{nVar})
        case 'title'
            title_ = varargin{nVar+1};
        case 'corrbolus'
            corrBolus_ = varargin{nVar+1};
        case {'secondary', 'plotsecondary', 'plotsecondaryinfusions'}
            plotSecondaryInfusions_ = varargin{nVar+1};
    end
end

if isempty(title_)
    if isfield(data, 'name') && isfield(data, 'startDate') && isfield(data, 'endDate')
        title(sprintf('%s ~ %s to %s', data.name, data.startDate, data.endDate), 'fontsize', 16, 'Interpreter', 'none');
    end
else
    title(title_, 'fontsize', 16, 'Interpreter', 'none');
end

%colours
blue = [0, 99, 198] / 255;
green = [99, 198, 0] / 255;
lightBlue = [50, 90, 132] / 255;
lightgreen = [90, 132, 50] / 255;
red = [255, 38, 54] / 255;
lightRed = [225, 118, 54] / 255;
orange = [220, 30, 30] / 255;
deeppink = [128, 0, 128] / 255;
magenta = [255, 20, 147] / 255;
grey = [100, 100, 100] / 255;

bolusDetailSize = 10;
mealBolusSize = 12;
mealBolusMarkerSize = 14;

minTextGap = 35 * duration_ / stepTime_ / (24 * 60);

% on the right axis
showInsulinAxis = nansum(data.basalInsulin) > 0 || nansum(data.pumpBasals) > 0 || nansum(data.steps) > 0 || nansum(data.bolusGlucagon) > 0;
if showInsulinAxis
    yyaxis(ax, 'right');
    cla(ax);
    hold(ax, 'on');
    
    if nansum(data.basalInsulin) > 0
        INS_BASAL = stairs(ax, data.time, data.basalInsulin, ...
            'color', blue, ...
            'LineWidth', 3.7, ...
            'LineStyle', '-', ...
            'Marker', 'none');
    end
    
    if plotSecondaryInfusions_ && nansum(data.secondaryBasalInsulin) > 0
        INS_BASAL = stairs(ax, data.time, data.secondaryBasalInsulin, ...
            'color', green, ...
            'LineWidth', 2.7, ...
            'LineStyle', '-', ...
            'Marker', 'none');
    end
    
    if nansum(data.pumpBasals) > 0
        INS_BASAL_OL = stairs(ax, data.time, data.pumpBasals, ...
            '--', ...
            'color', blue, ...
            'LineWidth', 2.7);
    end
    
    try
        if nansum(data.steps) > 0
            PHYSICAL_ACTIVITY = plot(ax, data.time, data.steps/200, ...
                'color', grey, ...
                'LineStyle', '-', ...
                'LineWidth', 3.7, ...
                'Marker', 'none');
        end
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
            throw(ME)
        end
    end
    
    try
        if nansum(data.bolusGlucagon) > 0
            for n = 1:1:length(data.bolusGlucagon)
                if data.bolusGlucagon(n) > 0
                    GLUCAGON = bar(data.time(n), data.bolusGlucagon(n)/2, ...
                        'BarWidth', 7.0, ...
                        'FaceColor', [0, 1, 0], ...
                        'LineStyle', 'none');
                end
            end
        end
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
            throw(ME)
        end
    end
    
    % set up right y axis
    %     if prctile(data.basalInsulin, 90) < 0.8 * 24 / 4
    %         ylim([0, 24 / 4]);
    %         yticks([0:1:6]);
    %     elseif prctile(data.basalInsulin, 90) < 0.8 * 24 / 3
    %         ylim([0, 24 / 3]);
    %         yticks([0:1:8]);
    %     elseif prctile(data.basalInsulin, 90) < 0.8 * 24 / 3
    %         ylim([0, 24 / 2]);
    %         yticks([0:2:12]);
    %     else
    %         ylim([0, 24]);
    %         yticks([0:4:24]);
    %     end
    ylim([0, 24 / 2]);
    yticks([0:2:12]);
    
    ylabelText = '';
    if nansum(data.basalInsulin) > 0
        ylabelText = [ylabelText, 'Insulin (U/h)'];
    end
    
    try
        if nansum(data.bolusGlucagon) > 0
            if isempty(ylabelText)
                ylabelText = [ylabelText, sprintf('Glucagon (0.5xU)')];
            else
                ylabelText = [ylabelText, sprintf(' & Glucagon (0.5xU)')];
            end
        end
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
            throw(ME)
        end
    end
    
    try
        if nansum(data.steps) > 0
            if isempty(ylabelText)
                ylabelText = [ylabelText, sprintf('Step Counts (x200)')];
            else
                ylabelText = [ylabelText, sprintf(' & Step Counts (x200)')];
            end
        end
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
            throw(ME)
        end
    end
    ylabel(ylabelText, 'color', 'k');
    
    % on the left axis
    yyaxis(ax, 'left');
    cla(ax);
end
hold(ax, 'on');

try
    if nansum(data.closedLoopActive) > 0
        CLActive = data.closedLoopActive;
        CLActive(isnan(CLActive)) = 0;
        idxStart = find(diff(CLActive > 0) == -1) + 1;
        if CLActive(1) <= 0
            idxStart = [1; idxStart];
        end
        idxEnd = find(diff(CLActive > 0) == 1);
        if CLActive(end) <= 0
            idxEnd = [idxEnd; length(data.time)];
        end
        
        for k = 1:length(idxStart)
            if CLActive(round((idxStart(k) + idxEnd(k))/2)) == 0
                patch(ax, ...
                    'XData', [data.time(idxStart(k)) - stepTime_ / 2 - 40 * (idxStart(k) == 1), ...
                    data.time(idxEnd(k)) + stepTime_ / 2 + 40 * (idxEnd(k) == length(data.time)), ...
                    data.time(idxEnd(k)) + stepTime_ / 2 + 40 * (idxEnd(k) == length(data.time)), ...
                    data.time(idxStart(k)) - stepTime_ / 2 - 40 * (idxStart(k) == 1)], ...
                    'YData', [25, 25, 0, 0], ...
                    'FaceColor', [0.1, 0.1, 0.1], ...
                    'FaceAlpha', 0.1, ...
                    'EdgeColor', 'none')
            else
                patch(ax, ...
                    'XData', [data.time(idxStart(k)) - stepTime_ / 2 - 40 * (idxStart(k) == 1), ...
                    data.time(idxEnd(k)) + stepTime_ / 2 + 40 * (idxEnd(k) == length(data.time)), ...
                    data.time(idxEnd(k)) + stepTime_ / 2 + 40 * (idxEnd(k) == length(data.time)), ...
                    data.time(idxStart(k)) - stepTime_ / 2 - 40 * (idxStart(k) == 1)], ...
                    'YData', [25, 25, 0, 0], ...
                    'FaceColor', [0.1, 0.1, 0.1], ...
                    'FaceAlpha', 0.1, ...
                    'EdgeColor', 'none')
            end
        end
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
end

% plot glucose indications
timeDense = (-40:1:duration_ + 40);
plot(ax, timeDense, ...
    10.0*ones(size(timeDense)), ...
    '--', 'color', grey, 'linewidth', 2.0, 'Marker', 'none');
plot(ax, timeDense, ...
    8.0*ones(size(timeDense)), ...
    '--', 'color', grey, 'linewidth', 2.0, 'Marker', 'none');
plot(ax, timeDense, 3.9*ones(size(timeDense)), '--', 'color', grey, 'linewidth', 2.0, 'Marker', 'none');

if any(contains(fieldnames(data), 'outcomeStartTime')) && ~isnan(data.outcomeStartTime)
    plot(ax, (data.outcomeStartTime - data.startTime)*[1, 1], ...
        [-100, 100], ...
        '-', 'color', 'k', 'linewidth', 3.0, 'Marker', 'none');
    patch(ax, ...
        'XData', [-40, ...
        data.outcomeStartTime - data.startTime, ...
        data.outcomeStartTime - data.startTime, ...
        -40], ...
        'YData', [25, 25, 0, 0], ...
        'FaceColor', [0.1, 0.1, 0.1], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
end

if any(contains(fieldnames(data), 'outcomeEndTime')) && ~isnan(data.outcomeEndTime)
    plot(ax, (data.outcomeEndTime - data.startTime)*[1, 1], ...
        [-100, 100], ...
        '-', 'color', 'k', 'linewidth', 3.0, 'Marker', 'none');
    patch(ax, ...
        'XData', [data.outcomeEndTime - data.startTime, ...
        duration_ + 40, ...
        duration_ + 40, ...
        data.outcomeEndTime - data.startTime], ...
        'YData', [25, 25, 0, 0], ...
        'FaceColor', [0.1, 0.1, 0.1], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
end

% plot Glucose
try
    if nansum(data.glucoseInterp) > 0
        plot(ax, data.time, data.glucoseInterp, ...
            'color', 'r', ...
            'linestyle', '-', ...
            'linewidth', 1.7, ...
            'Marker', 'none');
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
    
    if nansum(data.sensorGlucose) > 0
        plot(ax, data.time, data.sensorGlucose, ...
            'color', 'r', ...
            'linestyle', '-', ...
            'linewidth', 1.7, ...
            'Marker', 'none');
    end
end

if nansum(data.sensorGlucose) > 0
    GLUCOSE = plot(ax, data.time, data.sensorGlucose, ...
        'color', red, ...
        'LineStyle', 'none', ...
        'LineWidth', 2.7, ...
        'Marker', '.', ...
        'MarkerSize', 27);
end

try
    if nansum(data.sensorCalibration) > 0
        SENSOR_CALIBRATION = plot(ax, data.time, data.sensorCalibration, ...
            'color', red, ...
            'linestyle', 'none', ...
            'linewidth', 2.7, ...
            'Marker', 'o', ...
            'MarkerSize', 13);
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
end

try
    if nansum(data.bloodGlucose) > 0
        BLOOD_GLUCOSE = plot(ax, data.time, data.bloodGlucose, ...
            'color', red, ...
            'linestyle', 'none', ...
            'linewidth', 2.7, ...
            'Marker', '*', ...
            'MarkerSize', 21);
        
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
end

% plot Meals
mealTextPos = 21;
carbsIdx = find(data.carbs > 0);
plotted_carbs = zeros(size(data.time));
for n = carbsIdx(:)'
    delta_n = carbsIdx(abs(carbsIdx-n) < minTextGap);
    delta_time = (mean(data.time(delta_n)) - minTextGap * stepTime_ * mean(0:1:length(delta_n)-1)) + (0:1:length(delta_n) - 1) * minTextGap * stepTime_;
    for k = 1:length(delta_n)
        if plotted_carbs(delta_n(k))
            continue;
        end
        MEALS = plot(ax, data.time(delta_n(k)), mealTextPos+0.9, ...
            '-', ...
            'color', orange, ...
            'Marker', '^', ...
            'MarkerSize', mealBolusMarkerSize, ...
            'MarkerFaceColor', orange, ...
            'MarkerEdgeColor', orange);
        
        text(ax, delta_time(k), mealTextPos, [sprintf('%d', round(data.carbs(delta_n(k)))), 'g'], ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', ...
            'Color', orange, ...
            'FontSize', mealBolusSize, ...
            'FontWeight', 'bold');
    end
    plotted_carbs(delta_n) = true;
end

% plot Actual Meals
try
    if nansum(data.carbsActual) > 0
        actualMealTextPos = mealTextPos + 1.3;
        carbsIdx = find(data.carbsActual > 0 & data.carbsActual ~= data.carbs);
        plotted_actual_carbs = zeros(size(data.time));
        for n = carbsIdx(:)'
            if plotted_actual_carbs(n)
                continue;
            end
            delta_n = carbsIdx(abs(carbsIdx-n) < minTextGap);
            delta_time = (mean(data.time(delta_n)) - minTextGap * stepTime_ * mean(0:1:length(delta_n)-1)) + (0:1:length(delta_n) - 1) * minTextGap * stepTime_;
            for k = 1:length(delta_n)
                plot(ax, data.time(delta_n(k)), actualMealTextPos+0.9, ...
                    '-', ...
                    'color', magenta, ...
                    'Marker', '^', ...
                    'MarkerSize', mealBolusMarkerSize, ...
                    'MarkerFaceColor', magenta, ...
                    'MarkerEdgeColor', magenta);
                
                text(ax, delta_time(k), actualMealTextPos, [sprintf('%d', round(data.carbsActual(delta_n(k)))), 'g'], ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'center', ...
                    'Color', magenta, ...
                    'FontSize', mealBolusSize, ...
                    'FontWeight', 'bold');
            end
            plotted_actual_carbs(delta_n) = true;
        end
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
end

% plot Treatements
try
    treatsTextPos = mealTextPos - 3.6;
    treatsIdx = find(data.treats > 0);
    plotted_treats = zeros(size(data.time));
    for n = treatsIdx(:)'
        delta_n = treatsIdx(abs(treatsIdx-n) < minTextGap);
        delta_time = (mean(data.time(delta_n)) - minTextGap * stepTime_ * mean(0:1:length(delta_n)-1)) + (0:1:length(delta_n) - 1) * minTextGap * stepTime_;
        for k = 1:length(delta_n)
            if plotted_treats(delta_n(k))
                continue;
            end
            TREATS = plot(ax, data.time(delta_n(k)), treatsTextPos+0.9, ...
                '-', ...
                'color', deeppink, ...
                'Marker', '^', ...
                'MarkerSize', mealBolusMarkerSize, ...
                'MarkerFaceColor', deeppink, ...
                'MarkerEdgeColor', deeppink);
            
            text(ax, delta_time(k), treatsTextPos, [sprintf('%d', round(data.treats(delta_n(k)))), 'g'], ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', ...
                'Color', deeppink, ...
                'FontSize', mealBolusSize, ...
                'FontWeight', 'bold');
        end
        plotted_carbs(delta_n) = true;
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
end

% plot Insulin Boluses
try
    bolusTextPos = mealTextPos - 1.7;
    bolusIdx = find(data.bolusInsulin > 0);
    plotted_i_boluses = zeros(size(data.time));
    for n = bolusIdx(:)'
        delta_n = bolusIdx(abs(bolusIdx-n) < minTextGap);
        delta_time = (mean(data.time(delta_n)) - minTextGap * stepTime_ * mean(0:1:length(delta_n)-1)) + (0:1:length(delta_n) - 1) * minTextGap * stepTime_;
        for k = 1:length(delta_n)
            if plotted_i_boluses(delta_n(k))
                continue;
            end
            if ~isempty(data.bolusOverride) && data.bolusOverride(delta_n(k)) > 0
                INS_BOLUS = plot(ax, data.time(delta_n(k)), bolusTextPos-0.3, ...
                    '-', ...
                    'color', blue, ...
                    'Marker', 'v', ...
                    'MarkerSize', mealBolusMarkerSize, ...
                    'MarkerEdgeColor', orange, ...
                    'MarkerFaceColor', blue);
            else
                INS_BOLUS = plot(ax, data.time(delta_n(k)), bolusTextPos-0.3, ...
                    '-', ...
                    'color', blue, ...
                    'Marker', 'v', ...
                    'MarkerSize', mealBolusMarkerSize, ...
                    'MarkerEdgeColor', blue, ...
                    'MarkerFaceColor', blue);
            end
            
            totBolus = data.bolusInsulin(delta_n(k));
            if corrBolus_ && ~isempty(data.carbFactors) && ~isnan(data.carbFactors(delta_n(k))) && ~isempty(data.carbs) && ~isnan(data.carbs(delta_n(k)))
                %             if ~isnan(data.carbsActual(delta_n(k)))
                %                 mealBolus = data.carbsActual(delta_n(k)) / data.carbFactors(delta_n(k));
                %             else
                mealBolus = data.carbs(delta_n(k)) / data.carbFactors(delta_n(k));
                %             end
                if mealBolus > 0
                    text(ax, delta_time(k)-6.5, bolusTextPos+1.2, sprintf('\\div%3.1f', data.carbFactors(delta_n(k))), ...
                        'VerticalAlignment', 'bottom', ...
                        'HorizontalAlignment', 'center', ...
                        'Color', lightRed, ...
                        'FontSize', bolusDetailSize, ...
                        'FontWeight', 'bold');
                    
                    if round(totBolus, 1) ~= round(mealBolus, 1)
                        text(ax, delta_time(k)-6.5, bolusTextPos+0.7, sprintf('%+3.1f', round(totBolus, 1)-round(mealBolus, 1)), ...
                            'VerticalAlignment', 'bottom', ...
                            'HorizontalAlignment', 'center', ...
                            'Color', lightBlue, ...
                            'FontSize', bolusDetailSize, ...
                            'FontWeight', 'bold');
                    end
                    text(ax, delta_time(k)-1.5, bolusTextPos+0.65, '\_\_\_\_\_', ...
                        'VerticalAlignment', 'bottom', ...
                        'HorizontalAlignment', 'center', ...
                        'Color', lightBlue, ...
                        'FontSize', 10, ...
                        'FontWeight', 'bold');
                end
            end
            text(ax, delta_time(k)-3, bolusTextPos, sprintf('%4.1fU', totBolus), ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', ...
                'Color', blue, ...
                'FontSize', 0.95*mealBolusSize, ...
                'FontWeight', 'bold');
            
        end
        plotted_i_boluses(delta_n) = true;
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
end


% plot Secondary Insulin Boluses
if plotSecondaryInfusions_
    bolusTextPos = mealTextPos - 1.7;
    bolusIdx = find(data.secondaryBolusInsulin > 0);
    plotted_i_boluses = zeros(size(data.time));
    for n = bolusIdx(:)'
        delta_n = bolusIdx(abs(bolusIdx-n) < minTextGap);
        delta_time = (mean(data.time(delta_n)) - minTextGap * stepTime_ * mean(0:1:length(delta_n)-1)) + (0:1:length(delta_n) - 1) * minTextGap * stepTime_;
        for k = 1:length(delta_n)
            if plotted_i_boluses(delta_n(k))
                continue;
            end
            if ~isempty(data.bolusOverride) && data.bolusOverride(delta_n(k)) > 0
                INS_BOLUS_SEC = plot(ax, data.time(delta_n(k)), bolusTextPos-0.3, ...
                    '-', ...
                    'color', green, ...
                    'Marker', 'v', ...
                    'MarkerSize', mealBolusMarkerSize, ...
                    'MarkerEdgeColor', orange, ...
                    'MarkerFaceColor', green);
            else
                INS_BOLUS_SEC = plot(ax, data.time(delta_n(k)), bolusTextPos-0.3, ...
                    '-', ...
                    'color', green, ...
                    'Marker', 'v', ...
                    'MarkerSize', mealBolusMarkerSize, ...
                    'MarkerEdgeColor', green, ...
                    'MarkerFaceColor', green);
            end
            
            totBolus = data.secondaryBolusInsulin(delta_n(k));
            if corrBolus_ && ~isnan(data.carbFactors(delta_n(k)))
                %             if ~isnan(data.carbsActual(delta_n(k)))
                %                 mealBolus = data.carbsActual(delta_n(k)) / data.carbFactors(delta_n(k));
                %             else
                mealBolus = data.carbs(delta_n(k)) / data.carbFactors(delta_n(k));
                %             end
                if mealBolus > 0
                    text(ax, delta_time(k)-6.5, bolusTextPos+1.2, sprintf('\\div%3.1f', data.carbFactors(delta_n(k))), ...
                        'VerticalAlignment', 'bottom', ...
                        'HorizontalAlignment', 'center', ...
                        'Color', lightRed, ...
                        'FontSize', bolusDetailSize, ...
                        'FontWeight', 'bold');
                    
                    if round(totBolus, 1) ~= round(mealBolus, 1)
                        text(ax, delta_time(k)-6.5, bolusTextPos+0.7, sprintf('%+3.1f', round(totBolus, 1)-round(mealBolus, 1)), ...
                            'VerticalAlignment', 'bottom', ...
                            'HorizontalAlignment', 'center', ...
                            'Color', lightgreen, ...
                            'FontSize', bolusDetailSize, ...
                            'FontWeight', 'bold');
                    end
                    text(ax, delta_time(k)-1.5, bolusTextPos+0.65, '\_\_\_\_\_', ...
                        'VerticalAlignment', 'bottom', ...
                        'HorizontalAlignment', 'center', ...
                        'Color', lightgreen, ...
                        'FontSize', 10, ...
                        'FontWeight', 'bold');
                end
            end
            text(ax, delta_time(k)-3, bolusTextPos, sprintf('%4.1fU', totBolus), ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', ...
                'Color', green, ...
                'FontSize', 0.95*mealBolusSize, ...
                'FontWeight', 'bold');
            
        end
        plotted_i_boluses(delta_n) = true;
    end
end

% plot Insulin Basal Injections
try
    if nansum(data.basalInjection) > 0
        bolusTextPos = mealTextPos - 1.7;
        plotted_i_basal = zeros(size(data.time));
        for n = 1:1:length(data.time)
            if data.basalInjection(n) > 0
                INS_BASAL_INJEC = plot(ax, data.time(n), bolusTextPos-0.5, ...
                    '-', ...
                    'color', blue, ...
                    'Marker', 'o', ...
                    'MarkerSize', mealBolusMarkerSize, ...
                    'MarkerEdgeColor', blue, ...
                    'MarkerFaceColor', blue);
                
                delta_n = (n - round(minTextGap):n + round(minTextGap));
                delta_n(delta_n <= 0 | delta_n > length(data.time)) = [];
                if ~plotted_i_basal(n)
                    text(ax, data.time(n), bolusTextPos, [sprintf('%4.1f', nansum(data.basalInjection(delta_n))), 'U'], ...
                        'VerticalAlignment', 'bottom', ...
                        'HorizontalAlignment', 'center', ...
                        'Color', blue, ...
                        'FontSize', mealBolusSize, ...
                        'FontWeight', 'bold');
                    
                    plotted_i_boluses(delta_n) = ones(size(delta_n));
                end
            end
        end
    end
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
        throw(ME)
    end
end


% set up axis
ax = gca;
if duration_ <= 2 * 60
    sTick = 0.5 * 60;
elseif duration_ <= 8 * 60
    sTick = 1 * 60;
elseif duration_ <= 24 * 60
    sTick = 2 * 60;
elseif duration_ <= 48 * 60
    sTick = 4 * 60;
else
    sTick = 12 * 60;
end
ax.XTick = (sTick * floor((data.time(1) + startTime_)/(sTick)):sTick:sTick * ceil((data.time(end) + startTime_)/(sTick))) - startTime_;
ax.XTickLabel = [num2str(mod((sTick / 60 * floor((data.time(1) + startTime_)/(sTick)):sTick / 60:sTick / 60 * ceil((data.time(end) + startTime_)/(sTick))), 24)'), repmat(':00', length(sTick/60*floor((data.time(1) + startTime_)/(sTick)):sTick/60:sTick/60*ceil((data.time(end) + startTime_)/(sTick))), 1)];

set(gca, 'FontWeight', 'bold', 'LineWidth', 2.0);
ax.TickDir = 'out';
ax.XAxis.Color = 'k';
ax.YAxis(1).Color = 'k';
ax.YAxis(1).FontSize = 14;
if showInsulinAxis
    ax.YAxis(2).Color = 'k';
    ax.YAxis(2).FontSize = 14;
end
ax.XAxis.FontSize = 14;
xlim(data.time(1)+[- 40, duration_ + 40]);
ylim([0, 24]);
yticks((0:2:24))
xlabel('Time (HH:MM)');
ylabel('Sensor Glucose (mmol/L)');

legendHandlers = [];
legendTitles = {};
if exist('GLUCOSE', 'var') > 0
    legendHandlers(end+1) = GLUCOSE;
    legendTitles{end+1} = 'Sensor Glucose';
end
if exist('BLOOD_GLUCOSE', 'var') > 0
    legendHandlers(end+1) = BLOOD_GLUCOSE;
    legendTitles{end+1} = 'BG Meter';
end
if exist('SENSOR_CALIBRATION', 'var') > 0
    legendHandlers(end+1) = SENSOR_CALIBRATION;
    legendTitles{end+1} = 'Sensor Calibration';
end
if exist('INS_BASAL', 'var') > 0
    legendHandlers(end+1) = INS_BASAL;
    legendTitles{end+1} = 'Insulin Basal';
end
if exist('INS_BASAL_OL', 'var') > 0
    legendHandlers(end+1) = INS_BASAL_OL;
    legendTitles{end+1} = 'Insulin Basal (OL)';
end
if exist('INS_BASAL_INJEC', 'var') > 0
    legendHandlers(end+1) = INS_BASAL_INJEC;
    legendTitles{end+1} = 'Insulin Basal';
end
if exist('INS_BOLUS', 'var') > 0
    legendHandlers(end+1) = INS_BOLUS;
    legendTitles{end+1} = 'Insulin Bolus';
end
if exist('INS_BOLUS_SEC', 'var') > 0
    legendHandlers(end+1) = INS_BOLUS_SEC;
    legendTitles{end+1} = 'Insulin Bolus';
end
if exist('MEALS', 'var') > 0
    legendHandlers(end+1) = MEALS;
    legendTitles{end+1} = 'Meals';
end
if exist('TREATS', 'var') > 0
    legendHandlers(end+1) = TREATS;
    legendTitles{end+1} = 'Treatements';
end
if exist('PHYSICAL_ACTIVITY', 'var') > 0
    legendHandlers(end+1) = PHYSICAL_ACTIVITY;
    legendTitles{end+1} = '(500x) Steps';
end

lgd = legend(legendHandlers, legendTitles, 'location', 'northeastoutside', 'FontSize', 14);
legend boxoff;
drawnow

end
