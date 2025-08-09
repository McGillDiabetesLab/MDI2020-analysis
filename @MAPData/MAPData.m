classdef MAPData < MAPUtils
    %%MAPDATA Handle data collected by McGill Artificial Pancreas Lab
    
    properties (GetAccess = public, SetAccess = public)
        name = 'data' % Data set name
        units = 'uk'; % Choose between 'uk' or 'us'
        nested = 0 % if 1 a folder with the data 'name' will be created.
        outcome = true % if true, data will be used as part of outcomes.
        outcomeStartTime = NaN % if set outcome calculation start at this time, should be bigger than startTime
        outcomeEndTime = NaN % if set outcome calculation end at this time, should be smaller than startTime
        softMerge = true % if true, data is only added if was NaN or empty previously.
        useActualCarbs = true % report actual carbs if applicable
        useGlucoseInterp = true % if false don't use interpolated glucose for outcomes
        useBloodGlucoseInGlucoseInterp = false % if false don't use blood glucose entries in interpolated glucose for outcomes
        glucoseInterpHolesSize = 90; % maximum hole in the glucose data to interpolate
        glucoseInterpHolesSizeAfterMealsOrBoluses = 90 % maximum hole in the glucose data to interpolate when there is a meal or a bolus
        thresholdHypoglycemiaEvent = 3.9 % Threshold used to flag a hypoglycemia event
        durationHypoglycemiaEvent = 20 % Duration of time less than thresholdHypoglycemiaEvent to flag a hypoglycemia event
        nightInterval = [0, 6] * 60
        dayInterval = [6, 0] * 60
        breakfastInterval = [7, 10] * 60 % interval is inclusive in both extermities
        lunchInterval = [11.5, 14] * 60 % interval is inclusive in both extermities
        dinnerInterval = [16, 21] * 60 % interval is inclusive in both extermities
        snacksThreshold = 15 % threshold to consider a meal as a snack (inclusive)
        onlyCountBolusedMeals = true % in getMealInfo only meals accompagnied with bolus are considred
        onlyCountBiggerMeal = false % in getMealInfo only get the bigger meal in the specified interval
        mealInfoSnapshot = [-4 * 60, 6 * 60] % interval used for meal snapshots
        computeIOB = false % Flag to compute IOB, only set if needed IOB calculation slows down computations
        IMAPLogsReproductionMode = false; % In this mode, data are loaded for reproduction purposes
        plotSecondaryInfusions = false; % show or not show 
        listIgnoreFolders = {'obsolete', 'deprecated', 'ignore'};
        userData = struct(); % This is a free user data holder
    end
    
    properties (GetAccess = public, SetAccess = ?MAPPatient)
        patientHandle % The patient who had this data
        studyHandle % The study where this data was collected
    end
    
    properties (GetAccess = public, SetAccess = ?MAPPatient)
        patientName % The patient who had this data
        studyName % The study where this data was collected
    end
    
    properties (GetAccess = public, SetAccess = private)
        timeStamp % array for datetime
        sensorGlucose % array of sensor glucose (mmol/L)
        sensorScan % array recording when sensor is scanned (logical)
        sensorCalibration % array recording when sensor calibration happens (logical)
        sensor2Glucose % array of secondary sensor glucose (mmol/L)
        sensor2Calibration % array recording when secondary sensor calibration happens (logical)
        bloodGlucose % array of blood glucose (mmol/L)
        basalInsulin % infused basal insulin
        bolusInsulin % infused bolused insulin
        usualBolus % For fixed dose patients their usual bolus is saved here
        mealBolus % part of bolus considered to be meal carbohydrate
        corrBolus % part of bolus considered to be correction bolus
        bolusRecommend % algorithm bolus recommendation
        bolusExternal % manual bolus not recommended by algorithm
        bolusUser % bolus amount chosen by the user
        bolusStatus % SUCCESS - bolus delivered, FAILURE - bolus undelivered
        bolusOverride % bolus was overriden
        basalGlucagon % infused basal glucagon
        bolusGlucagon % infused basal glucagon
        basalPramlintide % infused basal pramlintide
        bolusPramlintide % infused basal pramlintide
        carbsActual % Actual Carbohydrate intake (g)
        carbs % carbohydrate intake (g)
        carbsFree % free carbohydrate intake (g)
        treats % treatment for hypoglycemia
        mealType % a meal parameter set by user, i.e. 1: breakfast, 2: lunch, ...
        tdd % used total daily dose
        exercise %1-exercise on, 0-exercise off
        targetGlucose % used target glucose
        targetAlgo % target that changes every 10 minutes within the algorithm
        carbFactors % programmed carb factors
        pumpBasals % programmed pump basals
        insulinSensitivityFactor % programmed insulin sensitivity factor
        basalInjection % injected long-acting insulin
        basalOverride % basal was overriden
        basalExternal % basal was external
        steps % array for number of steps per stepTime
        heartrate % heart rate
        appVersion
        algoVersion
        closedLoopActive % true if closed-loop is active
        gainR
        gainQ
        MOBfiltered % meal on board
        MOBpredicted
        IOB
        modelError
        modelIndex
        corr10mins %correction bolus calculated every 10 minutes
        secondaryBasalInsulin % infused basal insulin (secondary controller)
        secondaryBolusInsulin % infused bolused insulin (secondary controller)
        prevCarbs % Previously consumed carbs
        prevCarbsTime % Time when previous carbs was consumed
        prevBolus % Previously delivred insulin bolus
        prevBolusTime % Time when previous bolus was delivred
        internalState1 % This variable hold component 1 of the interenal state
        internalState2 % This variable hold component 2 of the interenal state
        internalState3 % This variable hold component 3 of the interenal state
        internalState4 % This variable hold component 4 of the interenal state
        internalState5 % This variable hold component 5 of the interenal state
        internalState6 % This variable hold component 6 of the interenal state
        internalState7 % This variable hold component 7 of the interenal state
        internalState8 % This variable hold component 8 of the interenal state
        tddBasal % tdd used for basal calculations
        tddBolus % tdd used for boluse calculations
    end
    
    properties (Dependent)
        timestampString
        time % array for time (min), starts at 0 and increments with stepTime
        iobBolus % Computed IOB Bolus using an expenentiol decay
        iobBasal % Computed IOB Basal using an expenentiol decay
        glucoseInterp % interpolated glucose (mmol/L)
        hypoEvents % a hypoglycemia event
        startTime % time start of data (min)
        duration % duration of the study (min)
        startDate % datetime object specifying the start of the study
        endDate % datetime object specifying the end of the study
        isEmpty % is Data empty ? (see how function is defined below)
        durationOutcome % Duration measured using outcome start and end
    end
    
    properties (GetAccess = public, SetAccess = immutable)
        stepTime % sampling time of data in min
        startDateLimit = NaT % Limit data to start with startDateLimit
        endDateLimit = NaT % Limit data to end with endDateLimit
    end
    
    methods (Static)
        % I had to explictly write this list to accelerate the code.
        % This list should be updated anytime a new variable is added.
        
        % getUnchangedFields: Fields that are computed online or never
        % changes
        function im = getUnchangedFields()
            im = { ...
                'timestampString', ...
                'time', ...
                'iobBolus', ...
                'iobBasal', ...
                'glucoseInterp', ...
                'hypoEvents', ...
                'stepTime', ...
                'startDateLimit', ...
                'endDateLimit', ...
                'startTime', ...
                'duration', ...
                'startDate', ...
                'endDate', ...
                'isEmpty', ...
                };
        end
        
        % getStaticFields: Fields that are time invarient
        function cst = getStaticFields()
            cst = { ...
                'appVersion', ...
                'algoVersion', ...
                'name', ...
                'nested', ...
                'units', ...
                'outcome', ...
                'outcomeStartTime', ...
                'outcomeEndTime', ...
                'softMerge', ...
                'useGlucoseInterp', ...
                'useBloodGlucoseInGlucoseInterp', ...
                'glucoseInterpHolesSize', ...
                'glucoseInterpHolesSizeAfterMealsOrBoluses', ...
                'thresholdHypoglycemiaEvent', ...
                'durationHypoglycemiaEvent', ...
                'nightInterval', ...
                'dayInterval', ...
                'breakfastInterval', ...
                'lunchInterval', ...
                'dinnerInterval', ...
                'snacksThreshold', ...
                'onlyCountBolusedMeals', ...
                'onlyCountBiggerMeal', ...
                'mealInfoSnapshot', ...
                'computeIOB', ...
                'startTime', ...
                'duration', ...
                'startDate', ...
                'endDate', ...
                'isEmpty', ...
                'stepTime', ...
                'startDateLimit', ...
                'endDateLimit', ...
                'studyName', ...
                'patientName', ...
                'userData', ...
                };
        end
        
        % getTimeFields: Fields that are time varient
        function signals = getTimeFields()
            signals = { ...
                'basalGlucagon', ...
                'basalInjection', ...
                'basalOverride', ...
                'basalExternal', ...
                'basalInsulin', ...
                'basalPramlintide', ...
                'bloodGlucose', ...
                'bolusExternal', ...
                'bolusGlucagon', ...
                'bolusInsulin', ...
                'bolusPramlintide', ...
                'bolusRecommend', ...
                'bolusStatus', ...
                'bolusUser', ...
                'bolusOverride', ...
                'usualBolus', ...
                'mealBolus', ...
                'corrBolus', ...
                'carbFactors', ...
                'carbs', ...
                'carbsFree', ...
                'carbsActual', ...
                'closedLoopActive', ...
                'exercise', ...
                'gainQ', ...
                'gainR', ...
                'glucoseInterp', ...
                'heartrate', ...
                'hypoEvents', ...
                'IOB', ...
                'MOBfiltered', ...
                'MOBpredicted', ...
                'modelError', ...
                'modelIndex', ...
                'corr10mins', ...
                'pumpBasals', ...
                'insulinSensitivityFactor', ...
                'sensor2Calibration', ...
                'sensor2Glucose', ...
                'sensorCalibration', ...
                'sensorGlucose', ...
                'sensorScan', ...
                'steps', ...
                'targetGlucose', ...
                'targetAlgo', ...
                'tdd', ...
                'timestampString', ...
                'time', ...
                'timeStamp', ...
                'treats', ...
                'mealType', ...
                'secondaryBasalInsulin', ...
                'secondaryBolusInsulin', ...
                'prevCarbs', ...
                'prevBolus', ...
                'prevCarbsTime', ...
                'prevBolusTime', ...
                'internalState1', ...
                'internalState2', ...
                'internalState3', ...
                'internalState4', ...
                'internalState5', ...
                'internalState6', ...
                'internalState7', ...
                'internalState8', ...
                'tddBasal',...
                'tddBolus',...
                };
        end
    end
    
    methods
        function out = get.outcomeStartTime(obj)
            if isnan(obj.outcomeStartTime) || ...
                    isempty(obj.outcomeStartTime) || ...
                    obj.outcomeStartTime < obj.startTime
                out = NaN;
            else
                out = obj.outcomeStartTime;
            end
        end
        
        function out = get.outcomeEndTime(obj)
            if isnan(obj.outcomeEndTime) || ...
                    isempty(obj.outcomeEndTime) || ...
                    obj.outcomeEndTime > obj.startTime + minutes(obj.duration)
                out = NaN;
            else
                out = obj.outcomeEndTime;
            end
        end
        
        function out = get.durationOutcome(obj)
            out = obj.duration;
            if ~isnan(obj.outcomeStartTime)
                out = out - minutes(obj.outcomeStartTime-obj.startTime);
            end
            if ~isnan(obj.outcomeEndTime)
                out = out + minutes(obj.outcomeEndTime-obj.startTime) - obj.duration;
            end
        end
        
        function out = get.timestampString(obj)
            out = [];
            if ~isempty(obj.timeStamp)
                out = datetime(obj.timeStamp, 'ConvertFrom', 'datenum');
            end
        end
        
        function t = get.time(obj)
            t = [];
            if ~isempty(obj.timeStamp)
                t = obj.stepTime * round((obj.timeStamp - obj.timeStamp(1))*24*60/obj.stepTime);
            end
        end
        
        function t = get.startTime(obj)
            t = NaN;
            if ~isnat(obj.startDate)
                t = obj.stepTime * round(mod(obj.timeStamp(1), 1)*24*60/obj.stepTime);
            end
        end
        
        function t = get.duration(obj)
            t = MAPUtils.roundToXmin(obj.endDate-obj.startDate+minutes(obj.stepTime), obj.stepTime);
        end
        
        function t = get.startDate(obj)
            t = NaT;
            if ~isempty(obj.timeStamp)
                t = datetime(obj.timeStamp(1)+1e-9, 'ConvertFrom', 'datenum');
            end
        end
        
        function t = get.endDate(obj)
            t = NaT;
            if ~isempty(obj.timeStamp)
                t = datetime(obj.timeStamp(end)+1e-9, 'ConvertFrom', 'datenum');
            end
        end
        
        function g = get.glucoseInterp(obj)
            g = obj.sensorGlucose;
            if ~isempty(g) && obj.useGlucoseInterp
                if obj.useBloodGlucoseInGlucoseInterp && nansum(obj.bloodGlucose) > 0
                    g(isnan(g)) = obj.bloodGlucose(isnan(g));
                end
                action = false(size(g));
                if nansum(obj.carbs) > 0
                    action = action | (obj.carbs > 0);
                end
                if nansum(obj.treats) > 0
                    action = action | (obj.treats > 0);
                end
                if nansum(obj.bolusInsulin) > 0
                    action = action | (obj.bolusInsulin > 0);
                end
                g(g == 0) = NaN;
                g = MAPUtils.smooth(g, action, obj.glucoseInterpHolesSize, obj.glucoseInterpHolesSizeAfterMealsOrBoluses, obj.stepTime);
            end
        end
        
        function e = get.hypoEvents(obj)
            if obj.useGlucoseInterp
                glucose_ = obj.glucoseInterp;
            else
                glucose_ = obj.sensorGlucose;
            end
            if nansum(obj.bloodGlucose) > 0 && nansum(obj.sensorGlucose) > 0
                minGlucose = min(obj.bloodGlucose, obj.sensorGlucose);
            elseif nansum(obj.bloodGlucose) > 0 && nansum(obj.sensorGlucose) == 0
                minGlucose = obj.bloodGlucose;
            elseif nansum(obj.bloodGlucose) == 0 && nansum(obj.sensorGlucose) > 0
                minGlucose = obj.sensorGlucose;
            else
                minGlucose = nan(size(glucose_));
            end
            if isempty(glucose_)
                glucose_ = minGlucose;
            else
                glucose_ = min(glucose_, minGlucose);
            end
            
            if nansum(obj.treats) > 0
                treats_ = obj.treats;
            else
                treats_ = zeros(size(glucose_));
            end
            [~, e] = MAPUtils.hypoCount(treats_, glucose_, obj.thresholdHypoglycemiaEvent, obj.durationHypoglycemiaEvent, obj.stepTime);
        end
        
        function iob = get.iobBasal(obj)
            iob = zeros(size(obj.time));
            if obj.computeIOB
                td = 4 * 60;
                tp = 75;
                
                for n = 1:length(iob)
                    for idx = find(obj.time <= obj.time(n) & obj.time(n) < obj.time+td)'
                        dt = obj.time(n) - obj.time(idx);
                        
                        tau = tp * (1 - tp / td) / (1 - 2 * tp / td);
                        a = 2 * tau / td;
                        S = 1 / (1 - a + (1 + a) * exp(-td/tau));
                        iob(n) = iob(n) + ...
                            (obj.stepTime / 60) * (obj.basalInsulin(idx) - obj.pumpBasals(idx)) * (1 - S * (1 - a) * ((dt.^2 / (tau * td * (1 - a)) - dt / tau - 1) .* exp(-dt/tau) + 1));
                    end
                end
            end
        end
        
        function iob = get.iobBolus(obj)
            iob = zeros(size(obj.time));
            if obj.computeIOB && nansum(obj.bolusInsulin) > 0
                td = 4 * 60;
                tp = 75;
                
                idxBoluses = find(obj.bolusInsulin > 0);
                
                for n = 1:length(iob)
                    for idx = idxBoluses(obj.time(idxBoluses) <= obj.time(n) & obj.time(n) < obj.time(idxBoluses)+td)'
                        dt = obj.time(n) - obj.time(idx);
                        
                        tau = tp * (1 - tp / td) / (1 - 2 * tp / td);
                        a = 2 * tau / td;
                        S = 1 / (1 - a + (1 + a) * exp(-td/tau));
                        iob(n) = iob(n) + ...
                            obj.bolusInsulin(idx) * (1 - S * (1 - a) * ((dt.^2 / (tau * td * (1 - a)) - dt / tau - 1) .* exp(-dt/tau) + 1));
                    end
                end
            end
        end
        
        function flag = get.isEmpty(obj)
            flag = true;
            if nansum(obj.sensorGlucose) > 0 || ...
                    nansum(obj.basalInsulin) > 0 || ...
                    nansum(obj.bolusInsulin) > 0 || ...
                    nansum(obj.basalInjection) > 0
                flag = false;
            end
        end
    end
    
    methods (Access = public)
        
        %% Constructor
        function obj = MAPData(varargin)
            obj.stepTime = 10;
            data = [];
            durationLimit = [];
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'name'
                        obj.name = varargin{nVar+1};
                    case 'steptime'
                        obj.stepTime = varargin{nVar+1};
                    case 'format'
                        obj.format = varargin{nVar+1};
                    case 'startdate'
                        obj.startDateLimit = varargin{nVar+1};
                    case 'enddate'
                        obj.endDateLimit = varargin{nVar+1};
                    case 'outcome'
                        obj.outcome = varargin{nVar+1};
                    case 'duration'
                        durationLimit = varargin{nVar+1};
                    case 'nested'
                        obj.nested = varargin{nVar+1};
                    case {'data', 'struct'}
                        data = varargin{nVar+1};
                    otherwise
                        error('[MAPData][publish] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
                end
            end
            
            if ~isempty(durationLimit) && ~isnat(obj.startDateLimit) && isnat(obj.endDateLimit)
                obj.endDateLimit = obj.startDateLimit + durationLimit - minutes(obj.stepTime);
            end
            
            if ~isempty(durationLimit) && ~isnat(obj.endDateLimit) && isnat(obj.startDateLimit)
                obj.startDateLimit = obj.endDateLimit - durationLimit + minutes(obj.stepTime);
            end
            
            if ~isempty(obj.endDateLimit)
                obj.endDateLimit = MAPUtils.floorToXmin(obj.endDateLimit, 1);
            end
            
            if ~isempty(obj.startDateLimit)
                obj.startDateLimit = MAPUtils.floorToXmin(obj.startDateLimit, 1);
            end
            
            if ~isempty(data)
                if isfield(data(1), 'stepTime')
                    obj.stepTime = data(1).stepTime;
                end
                obj.fromStruct(data);
            end
        end
        
        function clip(obj, varargin)
            fromIdx = 0;
            toIdx = inf;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'startdate'
                        fromIdx = find(abs(obj.timeStamp - datenum(varargin{nVar+1})) < obj.stepTime/60/24/2);
                    case 'enddate'
                        toIdx = find(abs(obj.timeStamp - datenum(varargin{nVar+1})) < obj.stepTime/60/24/2);
                    case 'fromidx'
                        fromIdx = varargin{nVar+1};
                    case 'toidx'
                        toIdx = varargin{nVar+1};
                    otherwise
                        error('[MAPData][clip] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
                end
            end
            
            if ~isempty(fromIdx) && fromIdx > 0
                fromTimeStamp_ = obj.timeStamp(fromIdx);
            else
                fromTimeStamp_ = obj.timeStamp(1);
            end
            
            if ~isempty(toIdx) && toIdx < inf
                toTimeStamp_ = obj.timeStamp(toIdx);
            else
                toTimeStamp_ = obj.timeStamp(end);
            end
            
            idx = obj.timeStamp < MAPUtils.roundToXmin(fromTimeStamp_, obj.stepTime, false) ...
                | obj.timeStamp > MAPUtils.roundToXmin(toTimeStamp_, obj.stepTime, false);
            
            for fn = setdiff(obj.getTimeFields, obj.getUnchangedFields)
                if isempty(obj.(fn{1}))
                    continue;
                end
                obj.(fn{1})(idx) = [];
            end
        end
        
        function dataOut = clear(obj, varargin)
            if nargout == 1
                dataOut = copy(obj);
                dataOut.clear(varargin{:});
                return
            end
            
            from_ = obj.startDate;
            to_ = obj.endDate;
            fields_ = setdiff(obj.getTimeFields, obj.getUnchangedFields);
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'from'
                        from_ = varargin{nVar+1};
                    case 'to'
                        to_ = varargin{nVar+1};
                    case 'fields'
                        fields_ = varargin{nVar+1};
                    otherwise
                        error('[MAPData][clear] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
                end
            end
            
            if ischar(fields_)
                fields_ = {fields_};
            end
            
            fromTimeStamp_ = datenum(from_);
            toTimeStamp_ = datenum(to_);
            
            idx = obj.timeStamp >= MAPUtils.roundToXmin(fromTimeStamp_, obj.stepTime, false) ...
                & obj.timeStamp <= MAPUtils.roundToXmin(toTimeStamp_, obj.stepTime, false);
            
            for fn = fields_(:)'
                if any(strcmp(fn{1}, {'time', 'timeStamp'}))
                    continue;
                end
                
                obj.(fn{1})(idx) = NaN;
            end
        end
        
        function dataOut = resize(obj, duration_, startTime_)
            if nargout == 1
                if numel(obj) > 1
                    if nargin < 3
                        dataOut = arrayfun(@(c)(c.resize(duration_)), obj, 'UniformOutput', false);
                    else
                        if numel(startTime_) == numel(obj)
                            dataOut = arrayfun(@(c, s)(c.resize(duration_, s)), obj, startTime_, 'UniformOutput', false);
                        else
                            dataOut = arrayfun(@(c)(c.resize(duration_, startTime_)), obj, 'UniformOutput', false);
                        end
                    end
                    dataOut = [dataOut{:}];
                    return;
                end
                dataOut = copy(obj);
                dataOut.resize(duration_, startTime_);
                return
            else
                if numel(obj) > 1
                    if nargin < 3
                        arrayfun(@(c)(c.resize(duration_, startTime_)), obj, 'UniformOutput', false);
                    else
                        if numel(startTime_) == numel(obj)
                            arrayfun(@(c, s)(c.resize(duration_, s)), obj, startTime_, 'UniformOutput', false);
                        else
                            arrayfun(@(c)(c.resize(duration_, startTime_)), obj, 'UniformOutput', false);
                        end
                    end
                    return;
                end
            end
            
            if nargin < 3
                startTime_ = obj.startTime;
            end

            obj.fromStruct(obj.toStruct('duration', duration_, 'startTime', startTime_, 'format', 'array'), false, true);
        end
        
        function data = getDays(obj, varargin)
            if numel(obj) > 1
                data = arrayfun(@(c)(c.getDays(varargin{:})), obj, 'UniformOutput', false);
                data = [data{:}];
                return;
            end
            
            startTime_ = obj.startTime;
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'starttime'
                        startTime_ = varargin{nVar+1};
                    otherwise
                        error('[MAPData][getDay] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
                end
            end
            
            for d = ceil(days(obj.duration)):-1:1
                data(d) = obj.getDay(d, 'starttime', startTime_);
            end
        end
        
        function data = getDay(obj, day_, varargin)
            startTime_ = obj.startTime;
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'starttime'
                        startTime_ = varargin{nVar+1};
                    otherwise
                        error('[MAPData][getDay] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
                end
            end
            if numel(obj) > 1
                data = arrayfun(@(c)(c.getDay(day_, varargin{:})), obj);
                return;
            end
            if isdatetime(day_)
                day_.Format = 'dd-MMM-yyyy';
                idx = find(obj.timestampString >= day_, 1);
                day_ = floor((obj.time(idx) + obj.startTime)/1440) + 1;
            end
            if any(day_ > 0)
                data = MAPData('data', obj.toStruct('days', day_, 'startTime', startTime_, 'format', 'array'));
            elseif any(day_ < 0)
                day_ = round(days(obj.duration)) + day_ + 1;
                data = MAPData('data', obj.toStruct('days', day_, 'startTime', startTime_, 'format', 'array'));
            else
                warning('[MAPData][getDay] Unrecognized input')
                return
            end
        end
        
        function data = getSnapshot(obj, from, to)            
            data = MAPData('data', obj.toStruct('startTime', from, 'duration', to-from));
        end
        
        function toPNG(obj, path, varargin)
            if mod(nargin, 2) == 1
                error('[MAPData][toPNG] You need to specify one argument which is the path of the png file.');
            end
            
            nested_ = obj.nested;
            inScreen_ = false;
            average_ = false;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'average'
                        average_ = varargin{nVar+1};
                    case 'nested'
                        nested_ = varargin{nVar+1};
                    case 'inscreen'
                        inScreen_ = varargin{nVar+1};
                end
            end
            
            figs = obj.toFigure( ...
                varargin{:}, ...
                'inScreen', inScreen_, ...
                'nested', nested_);
            
            for fig = figs
                if nested_ > 0
                    dPath = [path, obj.name, '/Day'];
                else
                    if isempty(obj.patientName) || contains([path, obj.name], obj.patientName)
                        dPath = [path, obj.name];
                    else
                        dPath = [path, obj.patientName, '_', obj.name];
                    end
                    if average_
                        dPath = [dPath, '_Summary'];
                    end
                end
                obj.ensureFolder(dPath);
                
                if isempty(fig.UserData)
                    print(fig, [dPath, '.png'], '-dpng');
                else
                    print(fig, [dPath, num2str(fig.UserData, '_%02d'), '.png'], '-dpng');
                end
                
                if ~inScreen_
                    close(fig);
                end
            end
        end
        
        
        function detailsTab = details(obj)
            row = 0;
            row = row + 1;
            detailsName{row, 1} = 'Start Date (dd-mmm-yyyy HH:MM)';
            if ~isnat(obj.startDate)
                if ~isnan(obj.outcomeStartTime)
                    detailsValue{row, 1} = datestr(obj.startDate+minutes(obj.outcomeStartTime-obj.startTime), 'dd-mmm-yyyy HH:MM');
                else
                    detailsValue{row, 1} = datestr(obj.startDate, 'dd-mmm-yyyy HH:MM');
                end
            else
                detailsValue{row, 1} = '';
            end
            
            row = row + 1;
            detailsName{row, 1} = 'End Date (dd-mmm-yyyy HH:MM)';
            if ~isnat(obj.endDate)
                if ~isnan(obj.outcomeEndTime)
                    detailsValue{row, 1} = datestr(obj.endDate+minutes(obj.outcomeEndTime-obj.startTime-minutes(obj.duration)), 'dd-mmm-yyyy HH:MM');
                else
                    detailsValue{row, 1} = datestr(obj.endDate, 'dd-mmm-yyyy HH:MM');
                end
            else
                detailsValue{row, 1} = '';
            end
            
            row = row + 1;
            detailsName{row, 1} = 'Duration (days)';
            detailsValue{row, 1} = days(obj.durationOutcome);
            
            indices = obj.intervalToIndices([0, 0]);
            
            row = row + 1;
            detailsName{row, 1} = 'Raw Sensor time (%)';
            detailsValue{row, 1} = obj.getSensorTime(indices, false);
            
            row = row + 1;
            detailsName{row, 1} = 'Interpolated Sensor time (%)';
            detailsValue{row, 1} = obj.getSensorTime(indices);
            
            detailsTab = table(detailsValue, 'RowNames', detailsName);
            detailsTab.Properties.VariableNames = {obj.name};
            detailsTab.Properties.VariableDescriptions = {'Data Name'};
        end
        
        function dataOut = excludeData(obj, timeStart, timeEnd)
            idx = obj.timestampString > timeStart & obj.timestampString < timeEnd;
            
            s = obj.toStruct('format', 'array');
            for fn = obj.getTimeFields
                s.(fn{1})(idx) = [];
            end
            
            dataOut = MAPData('data', s);
        end
        
        function outcomeTab = outcomes(obj, varargin)
            interval_ = [];
            glucoseOutcomes_ = true;
            insulinOutcomes_ = true;
            mealsOutcomes_ = days(obj.duration) > 0.99;
            otherOutcomes_ = [];
            dayNo_ = [];
            startTime_ = 0;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'interval'
                        interval_ = varargin{nVar+1};
                    case 'dayno'
                        dayNo_ = varargin{nVar+1};
                    case 'starttime'
                        startTime_ = varargin{nVar+1};
                    case 'glucoseoutcomes'
                        glucoseOutcomes_ = varargin{nVar+1};
                    case 'insulinoutcomes'
                        insulinOutcomes_ = varargin{nVar+1};
                    case 'mealsoutcomes'
                        mealsOutcomes_ = varargin{nVar+1};
                    case 'otheroutcomes'
                        otherOutcomes_ = varargin{nVar+1};
                end
            end
            
            if ~isempty(dayNo_)
                obj = obj.getDay(dayNo_, 'startTime', startTime_);
            end
            
            if obj.isEmpty == 1
                outcomeTab = [];
                return
            end
            
            if ~isempty(interval_) && isempty(obj.intervalToIndices(interval_))
                outcomeTab = [];
                return;
            end
            
            if glucoseOutcomes_
                glucose_ = obj.glucoseOutcomes(interval_);
            else
                glucose_ = [];
            end
            
            if insulinOutcomes_
                insulin_ = obj.insulinOutcomes(interval_);
            else
                insulin_ = [];
            end
            
            if mealsOutcomes_
                meals_ = obj.mealsOutcomes(interval_);
            else
                meals_ = [];
            end
            
            if ~isempty(otherOutcomes_) && isa(otherOutcomes_, 'function_handle')
                other_ = otherOutcomes_(obj, interval_);
            else
                other_ = [];
            end
            
            outcomeTab = [glucose_; insulin_; meals_; other_];
        end
        
        function obj_ = regroupMeals(obj, diffInMinutes, repeat)
            if nargin < 2
                diffInMinutes = 60;
            end
            if nargin < 3
                repeat = 10;
            end

            if repeat > 1
                for rr = 1:repeat
                    obj = obj.regroupMeals(diffInMinutes, 1);
                end
                obj_ = obj.copy();
                return;
            end

            obj_ = obj.copy();

            allMealTypes = unique(obj_.mealType(~isnan(obj_.mealType)));
            for k = 1:length(allMealTypes)
                idxMealAll = find(obj_.mealType == allMealTypes(k) | (isnan(obj_.mealType) & obj_.bolusInsulin > 0));
%                 idxMealAll = find(obj_.mealType == allMealTypes(k));
                idxClose = find(diff(idxMealAll) <= diffInMinutes / obj_.stepTime, 1);
                if isempty(idxClose)
                    continue;
                end
                
                n = idxMealAll(idxClose);
                m = idxMealAll(idxClose + 1);

                if ~isempty(obj_.carbs)
                    if isnan(obj_.carbs(n))
                        obj_.carbs(n) = obj_.carbs(m);
                    elseif ~isnan(obj_.carbs(m))
                        obj_.carbs(n) = obj_.carbs(m) + obj_.carbs(n);
                    end
                    obj_.carbs(m) = 0;
                end

                if isnan(obj_.bolusInsulin(n))
                    obj_.bolusInsulin(n) = obj_.bolusInsulin(m);
                elseif ~isnan(obj_.bolusInsulin(m))
                    obj_.bolusInsulin(n) = obj_.bolusInsulin(m) + obj_.bolusInsulin(n);
                end
                obj_.bolusInsulin(m) = 0;

                if isnan(obj_.mealBolus(n))
                    obj_.mealBolus(n) = obj_.mealBolus(m);
                elseif ~isnan(obj_.mealBolus(m))
                    obj_.mealBolus(n) = obj_.mealBolus(m) + obj_.mealBolus(n);
                end
                obj_.mealBolus(m) = 0;

                if isnan(obj_.corrBolus(n))
                    obj_.corrBolus(n) = obj_.corrBolus(m);
                elseif ~isnan(obj_.corrBolus(m))
                    obj_.corrBolus(n) = obj_.corrBolus(m) + obj_.corrBolus(n);
                end
                obj_.corrBolus(m) = 0;

                if isnan(obj_.bolusRecommend(n))
                    obj_.bolusRecommend(n) = obj_.bolusRecommend(m);
                elseif ~isnan(obj_.bolusRecommend(m))
                    obj_.bolusRecommend(n) = obj_.bolusRecommend(m) + obj_.bolusRecommend(n);
                end
                obj_.bolusRecommend(m) = 0;

                if isnan(obj_.bolusOverride(n))
                    obj_.bolusOverride(n) = obj_.bolusOverride(m);
                elseif ~isnan(obj_.bolusOverride(m))
                    obj_.bolusOverride(n) = max([obj_.bolusOverride(n),  obj_.bolusOverride(m)]);
                end
                obj_.bolusOverride(m) = NaN;

                if ~isempty(obj_.carbFactors)
                    if isnan(obj_.carbFactors(n))
                        obj_.carbFactors(n) = obj_.carbFactors(m);
                    end
                end

                if isnan(obj_.mealType(n))
                    obj_.mealType(n) = obj_.mealType(m);
                end
                obj_.mealType(m) = NaN;
            end
        end
    end
    
    methods (Access = private)
        function data = syncBoluses(obj, data)
            if nansum(obj.bolusInsulin) > 0 ...
                    && isfield(data, 'bolusInsulin') ...
                    && nansum(data.bolusInsulin) > 0
                for n = find(data.bolusInsulin > 0)'
                    if obj.bolusInsulin(n) == data.bolusInsulin(n)
                        continue;
                    end
                    if isfield(data, 'timeStamp')
                        m = find(abs(obj.timeStamp-data.timeStamp(n))*60*24 <= obj.stepTime/2);
                    else
                        m = find(abs(obj.time-data.time(n)) <= obj.stepTime/2);
                    end
                    if ~isempty(m)
                        interval = (-20:obj.stepTime:20)/obj.stepTime;
                        interval(interval == 0) = [];
                        mInterval = m + interval;
                        mInterval(mInterval <= 0) = [];
                        mInterval(mInterval > length(obj.timeStamp)) = [];
                        mIdx1 = mInterval(find(abs(obj.bolusInsulin(mInterval)-data.bolusInsulin(n)) < 1e-2, 1));
                        if ~isempty(mIdx1)
                            if ~isempty(data.carbs) && ~isnan(data.carbs(n)) % change obj to match data
                                if ~isempty(obj.carbs)
                                    mIdx2 = mInterval(find(abs(obj.carbs(mInterval)-data.carbs(n)) < 1e-2, 1));
                                    if ~isempty(mIdx2)
                                        obj.carbs(mIdx2) = 0.0;
                                        obj.carbs(m) = data.carbs(n);
                                    end
                                end
                                obj.bolusInsulin(mIdx1) = 0.0;
                                obj.bolusInsulin(m) = data.bolusInsulin(n);
                            else % change data to match obj
                                if isfield(data, 'timeStamp')
                                    nIdx = find(abs(data.timeStamp-obj.timeStamp(mIdx1))*60*24 <= obj.stepTime/2, 1);
                                else
                                    nIdx = find(abs(data.time-obj.time(mIdx1)) <= data.stepTime/2, 1);
                                end
                                data.bolusInsulin(n) = 0.0;
                                data.bolusInsulin(nIdx) = obj.bolusInsulin(mIdx1);
                                if obj.carbs(mIdx1) > 0 && data.carbs(nIdx) == 0 % TODO think about this
                                    data.carbs(nIdx) = obj.carbs(mIdx1);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
