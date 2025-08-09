function fromStruct(obj, s, softMerge_, reset_)
if nargin < 4
    reset_ = false;
end
if nargin < 3
    softMerge_ = obj.softMerge;
end

% convert from struct array to array
if numel(s) > 1
    s = MAPUtils.structMerge(s, obj.stepTime);
end

% Load time information
if ~isfield(s, 'timeStamp')
    if isfield(s, 'startDate') && isfield(s, 'time')
        s.timeStamp = datenum(s.startDate) + s.time / 60 / 24;
    elseif isfield(s, 'startDate') && isfield(s, 'duration')
        s.timeStamp = datenum(s.startDate):obj.stepTime:datenum(s.startDate+s.duration);
    elseif isfield(s, 'startDate') && isfield(s, 'endDate')
        s.timeStamp = datenum(s.startDate):obj.stepTime:datenum(s.endDate);
    elseif isfield(s, 'time') && ~isnat(obj.startDate)
        s.timeStamp = datenum(obj.startDate) + s.time / 60 / 24;
    else
        warning('[MAPData][fromStruct] Unable to load data from struct since timeStamp info could not be inferred.');
        return;
    end
end

% round seconds
s.timeStamp = round(s.timeStamp*(24*60))/(24*60);

% round to closet 10 min
if nansum(obj.timeStamp) > 0
    s.timeStamp = MAPUtils.roundToXmin([obj.timeStamp(1); s.timeStamp(:)], obj.stepTime);
    s.timeStamp(1) = [];
else
    s.timeStamp = MAPUtils.roundToXmin(s.timeStamp(:), obj.stepTime);
end

% Load Options
for fn = setdiff(obj.getStaticFields, obj.getUnchangedFields)
    if isfield(s, fn{1})
        obj.(fn{1}) = s.(fn{1});
    end
end

% Load fields
if isempty(obj.timeStamp) || reset_
    timeStamp_ = unique(s.timeStamp);
else
    timeStamp_ = unique([obj.timeStamp; s.timeStamp]);
end

timeStamp_ = MAPUtils.roundToXmin( ...
    (max(datenum(obj.startDateLimit)-1e-6, timeStamp_(find(~isnan(timeStamp_), 1))): ...
    obj.stepTime / (60 * 24): ...
    min(datenum(obj.endDateLimit)+1e-6, timeStamp_(find(~isnan(timeStamp_), 1, 'last'))))', min(1, obj.stepTime/10), 0);

for fn = setdiff(obj.getTimeFields, obj.getUnchangedFields)
    if strcmp(fn{1}, 'timeStamp')
        continue;
    end
    if reset_ && ~isstruct(s.(fn{1}))
        obj.(fn{1}) = [];
    end
    resizeArrayFlag = false;
    if isfield(s, fn{1})
        if isstruct(s.(fn{1})) % for now we don't have a way to merge struct fields
            if isstruct(obj.(fn{1})) % only do this of the obj is also struct
                obj.(fn{1}) = s.(fn{1});
            end
        elseif ~all(isnan(s.(fn{1}))) && ~isempty(s.(fn{1}))
            if (isempty(obj.(fn{1})) || all(isnan(obj.(fn{1})))) && length(timeStamp_) == length(s.(fn{1})) % shortcut
                temp_ = s.(fn{1});
            else % need to match indices
                temp_ = nan(size(timeStamp_));
                idxIn = ~isnan(s.(fn{1})(:)) & any(abs(s.timeStamp - timeStamp_') < obj.stepTime/(24 * 60)/2, 2);
                idxOutFromIn = any(abs(timeStamp_-s.timeStamp(idxIn)') < obj.stepTime/(24 * 60)/2, 2);
                if isempty(obj.(fn{1})) || all(isnan(obj.(fn{1}))) % Just copy all
                    temp_(idxOutFromIn) = s.(fn{1})(idxIn);
                else % if fn{1} is not empty in obj we need to merge
                    IdxCurrent = ~isnan(obj.(fn{1}));
                    idxOutFromCurrent = any(abs(timeStamp_-obj.timeStamp(IdxCurrent)') < obj.stepTime/(24 * 60)/2, 2);
                    if softMerge_ % only copies non NaN values over non NaN values
                        temp_(idxOutFromIn) = s.(fn{1})(idxIn);
                        temp_(idxOutFromCurrent) = obj.(fn{1})(IdxCurrent);
                    else % copies all non NaN values
                        temp_(idxOutFromCurrent) = obj.(fn{1})(IdxCurrent);
                        temp_(idxOutFromIn) = s.(fn{1})(idxIn);
                    end
                end
            end
            obj.(fn{1}) = temp_(:);
        else
            resizeArrayFlag = true;
        end
    else
        resizeArrayFlag = true;
    end
    
    if resizeArrayFlag 
        if ~isstruct(obj.(fn{1})) && ~isempty(obj.(fn{1}))
            temp_ = nan(size(timeStamp_));
            idxOutFromCurrent = any(abs(timeStamp_-obj.timeStamp') < obj.stepTime/(24 * 60)/2, 2);
            temp_(idxOutFromCurrent) = obj.(fn{1});
            obj.(fn{1}) = temp_(:);
        end
    end
    
end

obj.timeStamp = timeStamp_;
end
