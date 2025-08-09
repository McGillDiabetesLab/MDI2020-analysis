function s = toStruct(obj, varargin)
if isempty(obj)
    s = struct();
    return;
end

% Default options
format_ = 'array';
startTime_ = obj.startTime;
duration_ = inf;
days_ = [];
outcomes_ = false;

% Get options
for nVar = 1:2:length(varargin)
    switch lower(varargin{nVar})
        case 'format'
            format_ = varargin{nVar+1};
        case 'starttime'
            startTime_ = varargin{nVar+1};
        case 'duration'
            duration_ = varargin{nVar+1};
        case 'days'
            days_ = varargin{nVar+1};
        case 'outcomes'
            outcomes_ = varargin{nVar+1};
        otherwise
            error('[MAPData][toStruct] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
    end
end

if numel(obj) > 1 && (strcmp(format_, 'array') || strcmp(format_, 'average'))
    s = arrayfun(@(c)(c.toStruct(varargin{:})), obj, 'UniformOutput', false);
    s = [s{:}];
    return;
end

if outcomes_
    if ~isnan(obj.outcomeStartTime)
        startTime_ = obj.outcomeStartTime;
    end
    if ~isnan(obj.outcomeEndTime)
        duration_ = obj.outcomeEndTime - obj.outcomeStartTime + obj.stepTime;
    end
end

if isempty(startTime_)
    startTime_ = obj.startTime;
end

if isduration(duration_)
    duration_ = minutes(duration_);
end

if strcmp(format_, 'array')
    if ~isempty(days_)
        startTime_ = startTime_ + (days_(1) - 1) * 24 * 60;
        duration_ = (days_(end) - days_(1) + 1) * 24 * 60;
    end
    
    relTime = round((obj.time + obj.startTime - startTime_)/obj.stepTime) * obj.stepTime;
    
    s = struct();
    for fn = obj.getStaticFields
        if any(strcmp(fn{1}, {'name', 'startDate', 'endDate', 'duration'}))
            continue;
        end
        s.(fn{1}) = obj.(fn{1});
    end
    
    idx = 0 <= relTime & relTime < duration_;
    
    for fn = obj.getTimeFields
        if isempty(obj.(fn{1}))
            s.(fn{1}) = nan(sum(idx), 1);
        elseif isstruct(obj.(fn{1}))
            s.(fn{1}) = [];
        else
            s.(fn{1}) = obj.(fn{1})(idx);
        end
    end
    
    % pad or shrink data because of starttime
    if ~isempty(s.timeStamp) && obj.stepTime * round(mod(s.timeStamp(1), 1)*24*60/obj.stepTime) > mod(startTime_, 24*60)
        % add nan in other fields
        nanVectLen = (obj.stepTime * round(mod(s.timeStamp(1), 1)*24*60/obj.stepTime) - startTime_) / obj.stepTime;
        for fn = obj.getTimeFields
            if any(strcmp(fn{1}, {'timeStamp'}))
                continue;
            end
            if isdatetime(obj.(fn{1}))
                nanVect = NaT(nanVectLen, 1);
            else
                nanVect = nan(nanVectLen, 1);
            end
            s.(fn{1}) = [nanVect; s.(fn{1})];
        end
        % add fictional timeStamp
        actualStartTime = obj.stepTime * round(mod(s.timeStamp(1), 1)*24*60/obj.stepTime);
        if ~isempty(s.timeStamp)
            s.timeStamp = [s.timeStamp(1)-(actualStartTime-startTime_:-obj.stepTime:obj.stepTime)'/(1440); s.timeStamp];
        else
            s.timeStamp = datenum(obj.startDate) - (actualStartTime-startTime_-obj.stepTime:-obj.stepTime:0)'/(1440);
        end
    end
    
    % pad or shrink data because of duration
    if ~isinf(duration_) && obj.stepTime * length(s.timeStamp) < duration_
        nanVectLen = (duration_ - obj.stepTime * length(s.timeStamp)) / obj.stepTime;
        for fn = obj.getTimeFields
            if any(strcmp(fn{1}, {'timeStamp'}))
                continue;
            end
            if isdatetime(obj.(fn{1}))
                nanVect = NaT(nanVectLen, 1);
            else
                nanVect = nan(nanVectLen, 1);
            end
            s.(fn{1}) = [s.(fn{1}); nanVect];
        end
        % add fictional timeStamp
        actualDuration = obj.stepTime * length(s.timeStamp);
        if ~isempty(s.timeStamp)
            s.timeStamp = [s.timeStamp; s.timeStamp(end)+(obj.stepTime:obj.stepTime:duration_-actualDuration)'/(1440)];
        else
            s.timeStamp = datenum(obj.startDate) + (0:obj.stepTime:duration_-actualDuration-obj.stepTime)'/(1440);
        end
    end
    
    s.time = obj.stepTime * (0:1:(length(s.timeStamp) - 1))';
    s.name = obj.name;
    if s.outcomeStartTime < startTime_
        s.outcomeStartTime = NaN;
    end
    s.startTime = mod(startTime_, 1440);
    s.startDate = datetime(s.timeStamp(find(~isnan(s.timeStamp), 1))+1e-9, 'ConvertFrom', 'datenum');
    s.endDate = datetime(s.timeStamp(find(~isnan(s.timeStamp), 1, 'last'))+1e-9, 'ConvertFrom', 'datenum');
    s.duration = MAPUtils.roundToXmin(s.endDate-s.startDate+minutes(s.stepTime), s.stepTime);
else % daily or average
    % Extract data name
    if ~isempty(strfind(obj.name, '#'))
        warning('[MAPData][toStruct] Data is renamed from %s to %s.', obj.name, obj.name(1:strfind(obj.name, '#')-1));
        name_ = obj.name(1:strfind(obj.name, '#')-1);
    else
        name_ = obj.name;
    end
    
    % Construct structure
    s = struct([]);
    relTime = round((obj.time + obj.startTime - startTime_)/obj.stepTime) * obj.stepTime;
    if isempty(days_)
        nbrDays = ceil(relTime(end)/(24 * 60));
        days_ = 1:nbrDays;
    end
    for d = 1:length(days_)
        for fn = obj.getStaticFields
            if any(strcmp(fn{1}, {'name', 'startTime', 'startDate', 'endDate', 'duration'}))
                continue;
            end
            s(d).(fn{1}) = obj.(fn{1});
        end
        s(d).startTime = startTime_;
        s(d).time = (0:obj.stepTime:(24 * 60 - obj.stepTime))';
        for fn = obj.getTimeFields
            if strcmp(fn{1}, 'time')
                continue;
            end
            if isdatetime(obj.(fn{1}))
                s(d).(fn{1}) = NaT(24*60/obj.stepTime, 1);
            else
                s(d).(fn{1}) = nan(24*60/obj.stepTime, 1);
            end
            if ~isempty(obj.(fn{1}))
                idxR = any(relTime == (s(d).time + (days_(d) - 1) * 24 * 60)', 2);
                idxL = any((s(d).time + (days_(d) - 1) * 24 * 60) == relTime', 2);
                s(d).(fn{1})(idxL) = obj.(fn{1})(idxR);
            end
        end
        
        s(d).name = [name_, '#', num2str(days_(d), '%02d')];
        s(d).startDate = datetime(s(d).timeStamp(find(~isnan(s(d).timeStamp), 1))+1e-9, 'ConvertFrom', 'datenum');
        s(d).endDate = datetime(s(d).timeStamp(find(~isnan(s(d).timeStamp), 1, 'last'))+1e-9, 'ConvertFrom', 'datenum');
        s(d).duration = MAPUtils.roundToXmin(s(d).endDate-s(d).startDate+minutes(s(d).stepTime), s(d).stepTime);
    end
    
    % return average data if requested
    if strcmp(format_, 'average')
        sAll = s;
        s = struct();
        for fn = obj.getStaticFields
            if any(strcmp(fn{1}, {'name', 'startDate', 'endDate', 'duration'}))
                continue;
            end
            s.(fn{1}) = sAll(1).(fn{1});
        end
        s.time = (0:obj.stepTime:(24 * 60 - obj.stepTime))';
        s.timeStamp = sAll(1).timeStamp;
        for fn = obj.getTimeFields
            if strcmp(fn{1}, 'time')
                continue;
            end
            if strcmp(fn{1}, 'timeStamp')
                continue;
            end
            s.(fn{1}) = nanmean([sAll.(fn{1})], 2);
        end
        s.name = [name_, 'Summary'];
        s.startTime = startTime_;
        s.startDate = min([sAll.startDate]);
        s.endDate = max([sAll.endDate]);
        s.duration = max([sAll.duration]);
    end
end
end
