function [ax, lgd] = plotSummary(ax, data, varargin)
% data is a struct array or a cell of struct arrays

% if no ax just use the current one
if mod(nargin, 2) == 1
    if nargin > 1
        varargin(2:end+1) = varargin;
        varargin{1} = data;
    end
    data = ax;
    ax = gca;
end

% support multiple datasets
if iscell(data)
    M = numel(data);
else
    M = 1;
    data = {data};
end

% Default options
title_ = '';
startTime_ = mode([data{1}.startTime]);
duration_ = minutes(data{1}(1).duration);
averageType_ = 'median';
variabilityType_ = 'iqr';
variabilityStyle_ = 'patch';
if numel(data{1}) > 1
    smooth_ = 3;
else
    smooth_ = 1;
end
showInsulin_ = (M == 1 | numel(data{1}) == 1) & nansum(data{1}(1).basalInsulin) > 0;
hyposType_ = 'treats';
showhypos_ = false;
showCorrections_ = false;
showBolus_ = numel(data{1}) == 1;
showtitle_ = true;
showlegends_ = true;
grayscale_ = false;
armsName_ = cell(1, M);
for m = 1:M
    armsName_{m} = data{m}(1).name;
end
colorLevel_ = [];
incremental_ = 0;
markerGlucose_ = cellstr(repmat('none', M, 1));
defualtMatlabColor = [; ...
    0, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.9290, 0.6940, 0.1250; ...
    0.4940, 0.1840, 0.5560; ...
    0.4660, 0.6740, 0.1880; ...
    0.3010, 0.7450, 0.9330; ...
    0.6350, 0.0780, 0.1840];
colors_ = [];

% Get options
for nVar = 1:2:length(varargin)
    switch lower(varargin{nVar})
        case 'title'
            title_ = varargin{nVar+1};
        case 'starttime'
            startTime_ = varargin{nVar+1};
        case 'duration'
            duration_ = floor(varargin{nVar+1});
        case 'averagetype'
            averageType_ = varargin{nVar+1};
        case 'variabilitytype'
            variabilityType_ = varargin{nVar+1};
        case 'variabilitystyle'
            variabilityStyle_ = varargin{nVar+1};
        case 'smooth'
            smooth_ = varargin{nVar+1};
        case 'showinsulin'
            showInsulin_ = varargin{nVar+1};
        case 'showhypos'
            showhypos_ = varargin{nVar+1};
        case 'hypostype'
            hyposType_ = varargin{nVar+1};
        case 'showbolus'
            showBolus_ = varargin{nVar+1};
        case 'showcorrections'
            showCorrections_ = varargin{nVar+1};
        case 'showlegends'
            showlegends_ = varargin{nVar+1};
        case 'showtitle'
            showtitle_ = varargin{nVar+1};
        case 'grayscale'
            grayscale_ = varargin{nVar+1};
        case 'armsname'
            armsName_ = varargin{nVar+1};
        case 'colorlevel'
            colorLevel_ = varargin{nVar+1};
        case 'incremental'
            incremental_ = varargin{nVar+1};
        case 'markerglucose'
            markerGlucose_ = varargin{nVar+1};
        case 'colors'
            colors_ = varargin{nVar+1};
    end
end

if showtitle_
    switch averageType_
        case 'median'
            suffixTitle = 'Median';
        case 'mean'
            suffixTitle = 'Mean';
        otherwise
            error('Unkown average type');
    end
    switch variabilityType_
        case 'se'
            suffixTitle = [suffixTitle, ' (SE)'];
        case 'sd'
            suffixTitle = [suffixTitle, ' (SD)'];
        case 'iqr'
            suffixTitle = [suffixTitle, ' (IQR)'];
        otherwise
            error('Unkown average type');
    end
    if isempty(title_)
        title(suffixTitle, 'fontsize', 18);
    else
        title(sprintf('%s ~ %s', suffixTitle, title_), 'fontsize', 18, 'Interpreter', 'none');
    end
end

%colours
grey = [100, 100, 100] / 255;
if isempty(colors_)
    if M > 7 || grayscale_
        gColor = cell(M, 1);
        iColor = cell(M, 1);
        if grayscale_
            for m = 1:M
                gColor{m} = (5 / 7 - (4 / 7) * (m - 1) / (M - 1)) * ones(1, 3);
                iColor{m} = (5 / 7 - (4 / 7) * (m - 1) / (M - 1)) * ones(1, 3);
            end
        else
            parulaColor = parula(M);
            gColor = mat2cell(parulaColor(1:M, :), ones(1, M));
            iColor = mat2cell(parulaColor(1:M, :), ones(1, M));
        end
        
        if ~isempty(colorLevel_) && M == 1
            if ~grayscale_
                gColor{1} = [255 - 215 * colorLevel_, 85, 55] / 255;
                iColor{1} = [0, 0, 255 - 215 * colorLevel_] / 255;
            else
                warning('TODO');
                gColor{1} = (5 / 7 - (4 / 7) * colorLevel_) * ones(1, 3);
                iColor{1} = (5 / 7 - (4 / 7) * colorLevel_) * ones(1, 3);
            end
        end
    else
        if M == 1 && showInsulin_
            gColor = {[255, 0, 0] / 255};
            iColor = {[0, 0, 255] / 255};
        else
            gColor = mat2cell(defualtMatlabColor(1:M, :), ones(1, M));
            iColor = mat2cell(defualtMatlabColor(1:M, :), ones(1, M));
        end
    end
else
    gColor = mat2cell(colors_, ones(1, M));
    iColor = mat2cell(colors_, ones(1, M));
end

if incremental_ > 0
    showInsulin_ = false;
    showhypos_ = false;
    showBolus_ = false;
end

if showInsulin_
    yyaxis(ax, 'right');
    cla(ax);
    
    yyaxis(ax, 'left');
end
cla(ax);

mealTextPos = 20;
bolusTextPos = mealTextPos - 1.0;

for m = 1:M
    % get mean values
    startTime = data{m}(1).startTime;
%     nightInterval = data{m}(1).nightInterval;
    offsetTime = (floor(M/2) - m + 1) * 2.5;
    time = data{m}(1).time + offsetTime;
    name = armsName_{m};
    if contains(name, '#')
        name = name(1:strfind(name, '#')-1);
    end
    
    if showInsulin_
        % on the right axis
        yyaxis(ax, 'right');
        hold on;
        
        bolusAsBasal_ = false;
        if bolusAsBasal_
            basalInsulinAll = [data{m}.basalInsulin] + [data{m}.bolusInsulin]*6;
        else
            basalInsulinAll = [data{m}.basalInsulin];
        end
        
        if strcmp(averageType_, 'mean')
            basalInsulin = nanmean(basalInsulinAll, 2);
        else
            basalInsulin = nanmedian(basalInsulinAll, 2);
        end
        
        basalInsulinUp = smooth(prctile(basalInsulinAll, 75, 2), smooth_);
        
        basalInsulinDown = smooth(prctile(basalInsulinAll, 25, 2), smooth_);
        
        time_grade = kron(time(2:end), [1; 1]);
        time_grade(2:end+1) = time_grade;
        time_grade(1) = time(1);
        basal_grade = kron(basalInsulin(1:end-1), [1; 1]);
        basal_grade(end+1) = basalInsulin(end);
        basalInsulinUp_grade = kron(basalInsulinUp(1:end-1), [1; 1]);
        basalInsulinUp_grade(end+1) = basalInsulinUp(end);
        basalInsulinUp_grade(isnan(basalInsulinUp_grade)) = 0;
        basalInsulinDown_grade = kron(basalInsulinDown(1:end-1), [1; 1]);
        basalInsulinDown_grade(end+1) = basalInsulinDown(end);
        basalInsulinDown_grade(isnan(basalInsulinDown_grade)) = 0;
        patch([time_grade; flipud(time_grade)], [basalInsulinUp_grade; flipud(basalInsulinDown_grade)], ...
            iColor{m}, ...
            'FaceAlpha', 0.2, ...
            'LineStyle', 'none');
        INS_BASAL{m} = plot(time_grade, basal_grade, ...
            'color', iColor{m}, ...
            'LineStyle', '-', ...
            'linewidth', 1.6, ...
            'Marker', 'none');
        
        % set up right y axis
        if prctile(basalInsulin, 90) < 0.8 * 24 / 4
            ylim([0, 24 / 4]);
            yticks([0:1:6]);
        elseif prctile(basalInsulin, 90) < 0.8 * 24 / 3
            ylim([0, 24 / 3]);
            yticks([0:1:8]);
        elseif prctile(basalInsulin, 90) < 0.8 * 24 / 3
            ylim([0, 24 / 2]);
            yticks([0:2:12]);
        else
            ylim([0, 24]);
            yticks([0:4:24]);
        end
        ylim([0, 24 / 2]);
        yticks([0:1:12]);

        ylabelText = 'Insulin (u/h)';
        ylabel(ylabelText, 'color', 'k');
        % on the left axis
        yyaxis(ax, 'left');
    end
    
    hold on;
    if m == 1
        timeDense = (-40:1:duration_ + 40);
        if incremental_ > 0
            plot(ax, timeDense, zeros(size(timeDense)), '--', 'color', grey, 'linewidth', 2.0, 'Marker', 'none');
        else
            plot(ax, timeDense, ...
                10.0*ones(size(timeDense)), ...
                '--', 'color', grey, 'linewidth', 2.0, 'Marker', 'none');
%             plot(ax, timeDense, ...
%                 7.8*ones(size(timeDense)), ...
%                 '--', 'color', grey, 'linewidth', 2.0, 'Marker', 'none');
            plot(ax, timeDense, 3.9*ones(size(timeDense)), '--', 'color', grey, 'linewidth', 2.0, 'Marker', 'none');
        end
    end
    
    % plot Glucose
    glucose = [data{m}.glucoseInterp];
    if incremental_ > 0
        glucose = glucose - glucose(incremental_, :);
    end
    
    if strcmp(averageType_, 'mean')
        sensorGlucose = smooth(nanmean(glucose, 2), smooth_);
    else
        sensorGlucose = smooth(nanmedian(glucose, 2), smooth_);
    end
    if strcmpi(variabilityType_, 'sd')
        sensorGlucoseUp = nanmean(glucose, 2) + smooth(nanstd(glucose, [], 2), smooth_);
        sensorGlucoseDown = [];
    elseif strcmpi(variabilityType_, 'se')
        sensorGlucoseUp = nanmean(glucose, 2) + smooth(nanstd(glucose, [], 2)/sqrt(size(glucose, 2)), smooth_);
        sensorGlucoseDown = [];
    else
        sensorGlucoseUp = smooth(prctile(glucose, 75, 2), smooth_);
        sensorGlucoseDown = smooth(prctile(glucose, 25, 2), smooth_);
        sensorGlucoseUp(isnan(sensorGlucoseUp)) = 0.0;
        sensorGlucoseDown(isnan(sensorGlucoseDown)) = 0.0;
    end
    
    resample = 3;
    if isempty(sensorGlucoseDown)
        if strcmpi(variabilityStyle_, 'patch')
            patch(ax, [time; flipud(time)], [sensorGlucoseUp; flipud(sensorGlucoseUp)], ...
                gColor{m}, ...
                'FaceAlpha', 0.2, ...
                'LineStyle', 'none');
        elseif strcmpi(variabilityStyle_, 'bar')
            errorbar(ax, time(1:resample:end), sensorGlucose(1:resample:end), [], sensorGlucoseUp(1:resample:end)-sensorGlucose(1:resample:end), ...
                'color', gColor{m}, ...
                'linewidth', 3.0, ...
                'LineStyle', 'none');
        end
    else
        if strcmpi(variabilityStyle_, 'patch')
            patch(ax, [time; flipud(time)], [sensorGlucoseUp; flipud(sensorGlucoseDown)], ...
                gColor{m}, ...
                'FaceAlpha', 0.2, ...
                'LineStyle', 'none');
        elseif strcmpi(variabilityStyle_, 'bar')
            errorbar(ax, time, sensorGlucose, sensorGlucoseUp-sensorGlucose, sensorGlucoseDown-sensorGlucose, ...
                'color', gColor{m}, ...
                'LineStyle', 'none');
        end
    end
    
    GLUCOSE{m} = plot(ax, time, sensorGlucose, ...
        'color', gColor{m}, ...
        'LineStyle', '-', ...
        'LineWidth', 2.7, ...
        'Marker', '.', ...
        'MarkerSize', 22);
    
    plot(ax, time(1:resample:end), sensorGlucose(1:resample:end), ...
        'color', gColor{m}, ...
        'LineStyle', 'none', ...
        'LineWidth', 2.7, ...
        'Marker', markerGlucose_{m}, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', gColor{m});
    
    offsetBolusCarbs = (floor(M/2) - m + 1) * 0.35;
    if showBolus_
        meanCarbs = nanmean([data{m}.carbs], 2);
        numCarbs = sum([data{m}.carbs] > 0, 2);
        
        % plot Meals
        for n = find(meanCarbs(:) > 10)'
            plot(ax, time(n), mealTextPos+0.7+offsetBolusCarbs, ...
                '-', ...
                'color', gColor{m}, ...
                'Marker', '^', ...
                'MarkerSize', 6*(1 + log10(numCarbs(n))), ...
                'MarkerFaceColor', gColor{m}, ...
                'MarkerEdgeColor', gColor{m});
            text(ax, time(n), mealTextPos+offsetBolusCarbs, [sprintf('%d', round(nansum(meanCarbs(n)))), 'g'], ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', ...
                'Color', gColor{m}, ...
                'FontSize', 8, ...
                'FontWeight', 'bold');
        end
        
        meanBolus = nanmean([data{m}.bolusInsulin], 2);
        numBolus = sum([data{m}.bolusInsulin] > 0, 2);
        
        % plot Insulin Boluses
        for n = find(meanBolus(:) > 0.1)'
            plot(ax, time(n), bolusTextPos-0.3+offsetBolusCarbs, ...
                '-', ...
                'color', iColor{m}, ...
                'Marker', 'v', ...
                'MarkerSize', 6*(1 + log10(numBolus(n))), ...
                'MarkerEdgeColor', iColor{m}, ...
                'MarkerFaceColor', iColor{m});
            
            text(ax, time(n)+offsetBolusCarbs, bolusTextPos, [sprintf('%4.1f', meanBolus(n)), 'U'], ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', ...
                'Color', iColor{m}, ...
                'FontSize', 8, ...
                'FontWeight', 'bold');
        end
    elseif showCorrections_
        boluses = [data{m}.bolusInsulin];
        meals = [data{m}.carbs];
        corrBoluses = boluses;
        corrBoluses(meals > 0) = 0;
        meanBolus = nanmean(corrBoluses, 2);
        numBolus = sum(corrBoluses > 0, 2);
        % plot Insulin Boluses
        for n = find(meanBolus(:) > 0.1)'
            plot(ax, time(n), bolusTextPos-0.3+offsetBolusCarbs, ...
                '-', ...
                'color', iColor{m}, ...
                'Marker', 'v', ...
                'MarkerSize', 6*(1 + log10(numBolus(n))), ...
                'MarkerEdgeColor', iColor{m}, ...
                'MarkerFaceColor', iColor{m});
%             
%             text(ax, time(n)+offsetBolusCarbs, bolusTextPos, [sprintf('%4.1f', meanBolus(n)), 'U'], ...
%                 'VerticalAlignment', 'bottom', ...
%                 'HorizontalAlignment', 'center', ...
%                 'Color', iColor{m}, ...
%                 'FontSize', 8, ...
%                 'FontWeight', 'bold');
        end
    end
    
    if showhypos_
        if numel(data{m}) > 1
            hypoEventsPos = 5.5;
        else
            hypoEventsPos = 2.5;
        end
        offsetHypoEvents = (floor(M/2) - m + 1) * 0.35;
        switch hyposType_
            case {'treats', 'real', 'true'}
                hypoEvents = nanmean([data{m}.treats] > 0, 2);
            case {'calculated', 'computed', 'augmented', 'thresh'}
                hypoEvents = nanmean([data{m}.hypoEvents], 2);
        end
        for n = find(hypoEvents(:) > 0.25)'
            HYPO_EVENTS{m} = plot(ax, time(n), hypoEventsPos+offsetHypoEvents, ...
                'linestyle', 'none', ...
                'color', gColor{m}, ...
                'Marker', 'o', ...
                'MarkerSize', 15*(1 + log10(hypoEvents(n))), ...
                'MarkerEdgeColor', gColor{m}, ...
                'MarkerFaceColor', gColor{m});
        end
    end
    
    % set up axis
    ax = gca;
    if duration_ <= 2 * 60
        sTick = 30;
    elseif duration_ <= 8 * 60
        sTick = 1 * 60;
    elseif duration_ <= 24 * 60
        sTick = 2 * 60;
    elseif duration_ <= 48 * 60
        sTick = 4 * 60;
    else
        sTick = 12 * 60;
    end
    ax.XTick = (sTick * floor((time(1) + startTime_)/(sTick)):sTick:sTick * ceil((time(end) + startTime_)/(sTick))) - startTime_;
    if duration_ <= 2 * 60
        ax.XTickLabel = [num2str(floor(mod((sTick / 60 * floor((time(1) + startTime_)/(sTick)):sTick / 60:sTick / 60 * ceil((time(end) + startTime_)/(sTick))), 24)')),...
            repmat(':00', length(sTick/60*floor((time(1) + startTime_)/(sTick)):sTick/60:sTick/60*ceil((time(end) + startTime_)/(sTick))), 1)];
        halfHourLabels = [num2str(floor(mod((sTick / 60 * floor((time(1) + startTime_)/(sTick)):sTick / 60:sTick / 60 * ceil((time(end) + startTime_)/(sTick))), 24)')),...
            repmat(':30', length(sTick/60*floor((time(1) + startTime_)/(sTick)):sTick/60:sTick/60*ceil((time(end) + startTime_)/(sTick))), 1)];
        ax.XTickLabel(mod(sTick/60*floor((time(1) + startTime_)/(sTick)):sTick/60:sTick/60*ceil((time(end) + startTime_)/(sTick)), 1) ~= 0, :) = halfHourLabels(mod(sTick/60*floor((time(1) + startTime_)/(sTick)):sTick/60:sTick/60*ceil((time(end) + startTime_)/(sTick)), 1) ~= 0, :);
    else
        ax.XTickLabel = [num2str(mod((sTick / 60 * floor((time(1) + startTime_)/(sTick)):sTick / 60:sTick / 60 * ceil((time(end) + startTime_)/(sTick))), 24)'), repmat(':00', length(sTick/60*floor((time(1) + startTime_)/(sTick)):sTick/60:sTick/60*ceil((time(end) + startTime_)/(sTick))), 1)];
    end
    
    set(gca, 'FontWeight', 'bold', 'LineWidth', 2.0);
    ax.TickDir = 'out';
    ax.XAxis.Color = 'k';
    if showInsulin_
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        ax.YAxis(1).FontSize = 18;
        ax.YAxis(2).FontSize = 18;
    else
        ax.YAxis.Color = 'k';
        ax.YAxis.FontSize = 18;
    end
    ax.XAxis.FontSize = 18;
    if duration_ <= 4 * 60
        xlim(time(1)+[- 10, duration_ + 10]);
    else
        xlim(time(1)+[- 30, duration_ + 30]);
    end
    xlabel('Time (HH:MM)');
    
    if incremental_
        ylim([-4, 10]);
        yticks((-10:2:12));
        ylabel('Incremental sensor glucose (mmol/L)');
    else
        ylim([0, 24]);
        yticks((0:2:24));
        ylabel('Sensor glucose (mmol/L)');
    end
end

lgd = [];
if showlegends_
    legendHandlers = [];
    legendTitles = {};
    if exist('GLUCOSE', 'var') > 0
        for m = 1:length(GLUCOSE)
            if ~isempty(GLUCOSE{m})
                legendHandlers(end+1) = GLUCOSE{m};
                name = armsName_{m};
                if contains(name, '#')
                    name = name(1:strfind(name, '#')-1);
                end
                if ~showInsulin_ && ~showhypos_
                    legendTitles{end+1} = name;
                else
                    legendTitles{end+1} = ['Sensor glucose: ', name];
                end
            end
        end
    end
    if exist('INS_BASAL', 'var') > 0
        for m = 1:length(INS_BASAL)
            if ~isempty(INS_BASAL{m})
                name = armsName_{m};
                if contains(name, '#')
                    name = name(1:strfind(name, '#')-1);
                end
                legendHandlers(end+1) = INS_BASAL{m};
                legendTitles{end+1} = ['Basal insulin: ', name];
            end
        end
    end
    if exist('HYPO_EVENTS', 'var') > 0
        for m = 1:length(HYPO_EVENTS)
            if ~isempty(HYPO_EVENTS{m})
                name = armsName_{m};
                if contains(name, '#')
                    name = name(1:strfind(name, '#')-1);
                end
                legendHandlers(end+1) = HYPO_EVENTS{m};
                legendTitles{end+1} = ['Hypoglycemia events: ', name];
            end
        end
    end
    
    lgd = legend(legendHandlers, legendTitles, 'location', 'northeast', 'FontSize', 18);
    legend boxoff;
    drawnow
end
end
