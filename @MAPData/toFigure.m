function [fig, ax, lgd] = toFigure(obj, varargin)
if isempty(obj)
    return;
end

fig = [];
ax = [];
lgd = [];
if numel(obj) > 1
    arrayfun(@(c)(c.toFigure(varargin{:})), obj);
    return;
end

% Default options
inScreen_ = true;
startTime_ = obj.startTime;
duration_ = inf;
daily_ = false; % plot day per day, each day in figure
days_ = [];
title_ = '';
stats_ = true;
average_ = false; %
plotFun = @obj.plotData;

% Get options
for nVar = 1:2:length(varargin)
    switch lower(varargin{nVar})
        case 'inscreen'
            inScreen_ = varargin{nVar+1};
        case 'starttime'
            startTime_ = varargin{nVar+1};
        case 'duration'
            duration_ = varargin{nVar+1};
        case 'daily'
            daily_ = varargin{nVar+1};
        case 'days'
            days_ = varargin{nVar+1};
        case 'title'
            title_ = varargin{nVar+1};
        case 'stats'
            stats_ = varargin{nVar+1};
        case 'average'
            average_ = varargin{nVar+1};
        case {'plotfun', 'plotfct', 'plotfuncttion'}
            plotFun = varargin{nVar+1};
    end
end

% constuct default title
defaultTitle = '';
if ~isempty(obj.studyName)
    defaultTitle = [defaultTitle, sprintf('%s ~ ', obj.studyName)];
end
if ~isempty(obj.patientName)
    defaultTitle = [defaultTitle, sprintf('%s ~ ', obj.patientName)];
end
if ~isempty(obj.name)
    defaultTitle = [defaultTitle, sprintf('%s ~ ', obj.name)];
end
if ~isempty(defaultTitle)
    defaultTitle(end-1:end) = [];
end

% days argument overrides duration and daily
if ~isempty(days_)
    daily_ = true;
    duration_ = inf;
end

% average format overrides duration and daily
if average_
    daily_ = false;
    duration_ = inf;
end

% assume minutes if number
if ~isinf(duration_) && ~isduration(duration_)
    duration_ = minutes(duration_);
end

if daily_ && obj.duration > days(1)
    if isempty(days_)
        relTime = obj.time + obj.startTime - startTime_;
        days_ = 1:ceil(relTime(end)/(24 * 60));
    end
    
    fig = gobjects(1, length(days_));
    for k = 1:length(days_)
        if inScreen_
            fig(k) = figure(days_(k));
            clf;
        else
            fig(k) = figure('Visible', 'Off');
            clf;
        end
        if length(days_) > 1
            fig(k).UserData = days_(k);
        end
        
        plot_(fig(k), obj.getDay(days_(k), 'startTime', startTime_));
    end
else
    if inScreen_
        fig = figure(mod(prod(uint16(sprintf('%s%s%s', obj.studyName, obj.patientName, obj.name))), 7919));
        clf;
    else
        fig = figure('Visible', 'Off');
        clf;
    end
    
    if isinf(duration_)
        if average_ && ~isnan(obj.outcomeStartTime)
            plot_(fig, obj.resizeData(duration_, 'startTime', obj.outcomeStartTime));
        else
            plot_(fig, obj);
        end
    else
        plot_(fig, obj.resizeData(duration_, 'startTime', startTime_));
    end
end

    function plot_(fig, data)
        if isempty(data.time)
            close(fig);
            return;
        end
        
        set(fig, 'name', sprintf('MAPUtils::toFigure::%s', data.startDate), ...
            'numbertitle', 'off', ...
            'units', 'normalized', ...
            'outerposition', [0, 0, 1, 1], ...
            'defaultAxesColorOrder', [[1, 0, 0]; [0, 0, 1]]);
        
        if average_
            plotSummaryAndStats(fig, data);
        else
            plotDataAndStats(fig, data);
        end
    end

    function plotSummaryAndStats(fig, data)
        if stats_
            axisPositions = [0.04, 0.15, 0.75, 0.75];
        else
            axisPositions = [0.04, 0.15, 0.81, 0.75];
        end
        
        ax = subplot('Position', axisPositions);
        
        if isempty(title_)
            plotTitle = sprintf('%s ~ %s to %s', defaultTitle, data.startDate, data.endDate);
        else
            plotTitle = title_;
        end
        
        if isinf(duration_)
            data = data.toStruct('format', 'struct', 'starttime', startTime_, 'days', days_);
        else
            data = data.toStruct('format', 'struct', 'starttime', mod(startTime_, 24*60), 'duration', duration_, 'days', days_);
        end
        
        [ax, lgd] = MAPUtils.plotSummary(ax, ...
            data, ...
            varargin{:}, ...
            'title', plotTitle);
        ax.Position = axisPositions;
        if stats_
            MAPUtils.plotStatsSummary(fig, data)
        end
    end

    function plotDataAndStats(fig, data)
        if stats_
            axisPositions = [0.05, 0.15, 0.75, 0.75];
        else
            axisPositions = [0.05, 0.15, 0.75, 0.75];
        end
        
        ax = subplot('Position', axisPositions);
        
        if isempty(title_)
            plotTitle = sprintf('%s ~ %s to %s', defaultTitle, data.startDate, data.endDate);
        else
            plotTitle = title_;
        end
        
        if isinf(duration_)
            if startTime_ == data.startTime
                [ax, lgd] = plotFun(ax, ...
                    data, ...
                    varargin{:}, ...
                    'title', plotTitle);
            else
                [ax, lgd] = plotFun(ax, ...
                    data.toStruct('format', 'array', 'starttime', startTime_), ...
                    varargin{:}, ...
                    'title', plotTitle);
            end
        else
            [ax, lgd] = plotFun(ax, ...
                data.toStruct('format', 'array', 'starttime', mod(startTime_, 24*60), 'duration', duration_), ...
                varargin{:}, ...
                'title', plotTitle);
        end
        
        ax.Position = axisPositions;
        
        if stats_
            MAPUtils.plotStats(fig, data)
        end
    end
end