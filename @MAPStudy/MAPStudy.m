classdef MAPStudy < MAPUtils
    %MAPStudy Handle studies conducted by McGill Artificial Pancreas Lab
    
    properties (GetAccess = public, SetAccess = public)
        name = 'study' % Study name
        type = 'matched' % support crossover, matched (this is crossover without period effect) or parallel studies
        nested = 0 % if 1 a folder with the data 'name' will be created.
    end
    
    properties (GetAccess = public, SetAccess = private)
        patients = MAPPatient.empty % List of patients in the atudy
    end
    
    properties (Dependent)
        units
        isEmpty
        stepTime
        duration
        arms % Name of Arms
        useGlucoseInterp % if false don't use interpolated glucose for outcomes
        useBloodGlucoseInGlucoseInterp % if false don't use blood glucose entries in interpolated glucose for outcomes
        glucoseInterpHolesSize % maximum hole in the glucose data to interpolate
        glucoseInterpHolesSizeAfterMealsOrBoluses % maximum hole in the glucose data to interpolate when there is a meal or a bolus
        thresholdHypoglycemiaEvent % Threshold used to flag a hypoglycemia event
        durationHypoglycemiaEvent % Duration of time less than thresholdHypoglycemiaEvent to flag a hypoglycemia event
        nightInterval
        dayInterval
        breakfastInterval % interval is inclusive in both extermities
        lunchInterval % interval is inclusive in both extermities
        dinnerInterval % interval is inclusive in both extermities
        snacksThreshold % threshold to consider a meal as a snack (inclusive)
        onlyCountBolusedMeals
        onlyCountBiggerMeal
        mealInfoSnapshot % interval used for meal snapshots
        computeIOB
        outcomeStartTime % if set outcome calculation start at this time, should be bigger than startTime
        outcomeEndTime % if set outcome calculation end at this time, should be bigger than startTime
        softMerge % if true, data is only added if was NaN or empty previously.
    end
    
    methods
        function flag = get.isEmpty(obj)
            flag = isempty(obj.patients);
        end
        
        function set.name(obj, name_)
            obj.name = name_;
            if isprop(obj, 'patients')
                for pIdx = 1:length(obj.patients)
                    obj.patients(pIdx).studyName = obj.name;
                end
            end
        end
        
        function val = get.arms(obj)
            val = '';
            if ~isempty(obj.patients)
                val = [obj.patients.arms];
                val = unique(val);
            end
        end
        
        function dt = get.stepTime(obj)
            dt = [];
            if ~isempty(obj.patients) && ~isempty([obj.patients.stepTime])
                dt = mode([obj.patients.stepTime]);
                if any([obj.patients.stepTime] ~= dt)
                    error('Patients must have the same step time!');
                end
            end
        end
        
        function t = get.duration(obj)
            t = median([obj.patients.duration]);
        end
        
        function int = get.useGlucoseInterp(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).useGlucoseInterp;
            end
        end
        
        function set.useGlucoseInterp(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).useGlucoseInterp = int_;
            end
        end
        
        function int = get.useBloodGlucoseInGlucoseInterp(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).useBloodGlucoseInGlucoseInterp;
            end
        end
        
        function set.useBloodGlucoseInGlucoseInterp(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).useBloodGlucoseInGlucoseInterp = int_;
            end
        end
        
        function int = get.glucoseInterpHolesSize(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).glucoseInterpHolesSize;
            end
        end
        
        function set.glucoseInterpHolesSize(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).glucoseInterpHolesSize = int_;
            end
        end
        
        function int = get.glucoseInterpHolesSizeAfterMealsOrBoluses(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).glucoseInterpHolesSizeAfterMealsOrBoluses;
            end
        end
        
        function set.glucoseInterpHolesSizeAfterMealsOrBoluses(obj, val_)
            for p = 1:numel(obj.patients)
                obj.patients(p).glucoseInterpHolesSizeAfterMealsOrBoluses = val_;
            end
        end
        
        function int = get.thresholdHypoglycemiaEvent(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).thresholdHypoglycemiaEvent;
            end
        end
        
        function set.thresholdHypoglycemiaEvent(obj, val_)
            for p = 1:numel(obj.patients)
                obj.patients(p).thresholdHypoglycemiaEvent = val_;
            end
        end
        
        function int = get.durationHypoglycemiaEvent(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).durationHypoglycemiaEvent;
            end
        end
        
        function set.durationHypoglycemiaEvent(obj, val_)
            for p = 1:numel(obj.patients)
                obj.patients(p).durationHypoglycemiaEvent = val_;
            end
        end
        
        function int = get.nightInterval(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).nightInterval;
            end
        end
        
        function set.nightInterval(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).nightInterval = int_;
            end
        end
        
        function int = get.dayInterval(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).dayInterval;
            end
        end
        
        function set.dayInterval(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).dayInterval = int_;
            end
        end
        
        function int = get.breakfastInterval(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).breakfastInterval;
            end
        end
        
        function set.breakfastInterval(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).breakfastInterval = int_;
            end
        end
        
        function int = get.lunchInterval(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).lunchInterval;
            end
        end
        
        function set.lunchInterval(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).lunchInterval = int_;
            end
        end
        
        function int = get.dinnerInterval(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).dinnerInterval;
            end
        end
        
        function set.dinnerInterval(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).dinnerInterval = int_;
            end
        end
        
        function int = get.snacksThreshold(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).snacksThreshold;
            end
        end
        
        function set.snacksThreshold(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).snacksThreshold = int_;
            end
        end
        
        function int = get.onlyCountBolusedMeals(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).onlyCountBolusedMeals;
            end
        end
        
        function set.onlyCountBolusedMeals(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).onlyCountBolusedMeals = int_;
            end
        end
        
        function int = get.onlyCountBiggerMeal(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).onlyCountBiggerMeal;
            end
        end
        
        function set.onlyCountBiggerMeal(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).onlyCountBiggerMeal = int_;
            end
        end
        
        function int = get.mealInfoSnapshot(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).mealInfoSnapshot;
            end
        end
        
        function set.mealInfoSnapshot(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).mealInfoSnapshot = int_;
            end
        end
        
        function int = get.computeIOB(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).computeIOB;
            end
        end
        
        function set.computeIOB(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).computeIOB = int_;
            end
        end
        
        function int = get.outcomeStartTime(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).outcomeStartTime;
            end
        end
        
        function set.outcomeStartTime(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).outcomeStartTime = int_;
            end
        end
        
        function int = get.outcomeEndTime(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).outcomeEndTime;
            end
        end
        
        function set.outcomeEndTime(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).outcomeEndTime = int_;
            end
        end
        
        function int = get.units(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).units;
            end
        end
        
        function set.units(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).units = int_;
            end
        end
        
        function int = get.softMerge(obj)
            int = [];
            if ~isempty(obj.patients)
                int = obj.patients(1).softMerge;
            end
        end
        
        function set.softMerge(obj, int_)
            for p = 1:numel(obj.patients)
                obj.patients(p).softMerge = int_;
            end
        end
    end
    
    methods (Static)
        function cst = getUnchangedFields()
            cst = { ...
                'stepTime', ...
                'duration', ...
                'arms', ...
                };
        end
        function cst = getStaticFields()
            cst = { ...
                'name', ...
                'type', ...
                'nested', ...
                'stepTime', ...
                'duration', ...
                'arms', ...
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
                'computeIOB', ...
                'outcomeStartTime', ...
                'outcomeEndTime', ...
                'softMerge', ...
                };
        end
    end
    
    methods (Access = public)
        
        %% Constructor
        function obj = MAPStudy(varargin)
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'name'
                        obj.name = varargin{nVar+1};
                    case 'type'
                        obj.type = varargin{nVar+1};
                    case 'mat'
                        obj.fromMAT(varargin{nVar+1});
                    otherwise
                        error('[MAPData][MAPStudy] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
                end
            end
        end
        
        function add(obj, mapPatient)
            if ~isa(mapPatient, 'MAPPatient')
                error('[MAPStudy][add] Only supports adding a MAPStudy object.');
            end
            
            idxPatient = find([obj.patients.ID] == mapPatient.ID, 1);
            
            if isempty(idxPatient)
                obj.patients(end+1) = mapPatient.copy;
                obj.patients(end).studyName = obj.name;
            else
                obj.patients(idxPatient) = mapPatient.copy;
                obj.patients(idxPatient).studyName = obj.name;
            end
        end
        
        function remove(obj, idx, id)
            if isnumeric(idx) || islogical(idx)
                obj.patients(idx) = [];
            elseif strcmpi(idx, 'id')
                obj.patients(any([obj.patients.ID] == id(:), 1)) = [];
            else
                error('Unknown input');
            end
        end
        
        function res = hasPatient(obj, id)
            res = ~isempty(obj.getPatient(id));
        end
        
        function mapPatient = getPatient(obj, id)
            if nargin < 2
                disp(strcat(num2str((1:length(obj.patients))', '%d -> '), ...
                    {obj.patients.name}', ...
                    cellfun(@(c)(strjoin([{' ->'}, c], ' ')), {obj.patients.arms}, 'UniformOutput', false)'));
                return;
            end
            
            mapPatient = MAPPatient.empty;
            if isempty([obj.patients.ID])
                return;
            end
            
            if ischar(id)
                idNum = regexp(id, '\d+.?', 'match');
                if isempty(idNum)
                    return;
                end
                idNum = str2double(idNum(end));
                mapPatient = obj.patients([obj.patients.ID] == idNum);
            else
                for idx = numel(id):-1:1
                    if iscell(id(idx))
                        idNum = regexp(id{idx}, '\d+.?', 'match');
                        if isempty(idNum)
                            continue;
                        end
                        idNum = str2double(idNum(end));
                    elseif ischar(id(idx))
                        idNum = regexp(id, '\d+.?', 'match');
                        if isempty(idNum)
                            continue;
                        end
                        idNum = str2double(idNum(end));
                    else
                        idNum = id(idx);
                    end
                    if any([obj.patients.ID] == idNum)
                        mapPatient(idx) = obj.patients([obj.patients.ID] == idNum);
                    end
                end
            end
        end
        
        function mapData = getData(obj, name)
            if nargout == 1
                mapData = MAPData.empty;
            end
            if nargin < 2
                disp(obj.arms);
                return;
            end
            mapData = obj.patients.getData(name);
        end
        
        %% Functions to load data
        function fromStruct(obj, s)
            for fn = setdiff(obj.getStaticFields, obj.getUnchangedFields)
                if isfield(s, fn{1})
                    obj.(fn{1}) = s.(fn{1});
                end
            end
            
            if isfield(s, 'patients')
                nbrPatients = numel(s.patients);
                if nbrPatients > 0
                    for p = nbrPatients:-1:1
                        map = MAPPatient('patient', s.patients(p));
                        map.studyName = obj.name;
                        if ~map.isEmpty
                            obj.add(map);
                        end
                    end
                end
            end
        end
        
        function fromSimulator(obj, s, varargin)
            dataName_ = cell(1, length(s.patients));
            startDate_ = NaT;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'dataname'
                        dataName_ = varargin{nVar+1};
                    case 'startdate'
                        startDate_ = varargin{nVar+1};
                    otherwise
                        error('[MAPData][fromSimulator] Unkown option %s, use doc MAPStudy for more information.', varargin{nVar});
                end
            end
            
            if numel(dataName_) == 1 && ~isempty(dataName_{1})
                dataName_ = cellstr(repmat(dataName_, [length(s.patients), 1]));
            end
            
            % loop through patients
            for idx = 1:length(s.patients)
                sData = struct();
                sData.stepTime = s.options.simulationStepSize;
                sData.time = 0:s.options.simulationStepSize:(s.options.simulationDuration - s.options.simulationStepSize);
                if isnat(startDate_)
                    startDate_ = datetime(datestr(floor(now)+(s.options.simulationStartTime)/(24 * 60)));
                end
                sData.startDate = startDate_;
                utime = cell2mat(s.resultsManagers{idx}.primaryInfusionTimes);
                sData.sensorGlucose = cell2mat(s.resultsManagers{idx}.glucoseMeasurements);
                sData.sensorGlucose(end) = [];
                uval = cell2mat(s.resultsManagers{idx}.primaryInfusions);
                if ~isempty(uval)
                    sData.basalInsulin = full([uval.basalInsulin]);
                    sData.bolusInsulin = full([uval.bolusInsulin]);
                    if isfield(uval, 'bolusGlucagon')
                        sData.bolusGlucagon = [uval.bolusGlucagon];
                    end
                end
                sData.carbsActual = full(s.resultsManagers{idx}.patient.mealPlan.getMeal(utime).value);
                sData.carbs = full(s.resultsManagers{idx}.patient.getMeal(utime).value);
                sData.treats = full(s.resultsManagers{idx}.patient.getTreatment(utime));
                if isprop(s.resultsManagers{idx}, 'secondaryInfusions')
                    uval = cell2mat(s.resultsManagers{idx}.secondaryInfusions);
                    if ~isempty(uval)
                        sData.secondaryBasalInsulin = full([uval.basalInsulin]);
                        sData.secondaryBolusInsulin = full([uval.bolusInsulin]);
                    end
                end
                prop = s.patients{idx}.getProperties();
                if length(sData.time) == length(prop.pumpBasals.time) && all(sData.time(1)+s.options.simulationStartTime == prop.pumpBasals.time(1))
                    sData.pumpBasals = prop.pumpBasals.value;
                else
                    sData.pumpBasals = obj.struct2array(prop.pumpBasals, sData.time(:)+s.options.simulationStartTime);
                end
                if length(sData.time) == length(prop.carbFactors.time) && all(sData.time(1)+s.options.simulationStartTime == prop.carbFactors.time(1))
                    sData.carbFactors = prop.carbFactors.value;
                else
                    sData.carbFactors = obj.struct2array(prop.carbFactors, sData.time(:)+s.options.simulationStartTime);
                end
                
                if isempty(dataName_{idx})
                    dataName_{idx} = s.primaryControllers{idx}.name;
                end
                pName_ = s.patients{idx}.name;
                if isempty(obj.getPatient(pName_))
                    pMAP = MAPPatient('name', pName_);
                    if isfield(s.patients{idx}.param, 'TDD')
                        pMAP.TDD = s.patients{idx}.param.TDD;
                    end
                    if isfield(s.patients{idx}.param, 'w')
                        pMAP.weight = s.patients{idx}.param.w;
                    end
                    mData = MAPData('name', dataName_{idx}, 'data', sData);
                    pMAP.add(mData);
                    obj.add(pMAP);
                else
                    if isempty(obj.getPatient(pName_).getData(dataName_{idx}))
                        mData = MAPData('name', dataName_{idx}, 'data', sData);
                        obj.getPatient(pName_).add(mData);
                    else
                        obj.getPatient(pName_).getData(dataName_{idx}).fromStruct(sData, false);
                    end
                end
            end
        end
        
        function makeData(obj, name_, varargin)
            for p = 1:numel(obj.patients)
                obj.patients(p).makeData(name_, varargin{:});
            end
        end
        
        function removeData(obj, name_)
            for p = 1:numel(obj.patients)
                obj.patients(p).removeData(name_);
            end
        end
        
        %% Functions to save data
        function s = toStruct(obj, varargin)
            patients_ = [obj.patients.ID];
            average_ = false;
            arms_ = obj.arms;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'average'
                        average_ = varargin{nVar+1};
                    case 'patients'
                        patients_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = varargin{nVar+1};
                end
            end
            
            s = struct();
            
            for fn = obj.getStaticFields()
                s.(fn{1}) = obj.(fn{1});
            end
            
            for p = 1:numel(patients_)
                if ~isempty(intersect({obj.getPatient(patients_(p)).data.name}, arms_))
                    s.patients(p) = obj.getPatient(patients_(p)).toStruct( ...
                        varargin{:}, ...
                        'average', average_);
                end
            end
        end
        
        function demoSheet = getDemographSheet(obj, varargin)
            patients_ = [obj.patients.ID];
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'patients'
                        patients_ = varargin{nVar+1};
                end
            end
            
            % Reorder indices
            [~, idx] = sort({obj.patients.name});
            idxValid = any([obj.patients.ID] == patients_(:));
            idx = idx(idxValid);
            
            % Demographics sheet
            demoSheet = cell(0);
            
            row = 1;
            demoSheet{row, 1} = 'Study';
            demoSheet{row, 2} = obj.name;
            
            row = row + 1;
            col = 0;
            len = numel(idx);
            
            data_ = [obj.patients(idx).ID]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Patient ID';
                demoSheet{row+1, col} = 'ID';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%4d', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = {'Mean'; 'SD'; 'Median'; '25_IQR'; '75_IQR'; 'Max'; 'Min'};
            end
            
            data_ = [obj.patients(idx).admissionDate]';
            if ~all(isnat(data_))
                col = col + 1;
                demoSheet{row, col} = 'Admission date';
                demoSheet{row+1, col} = 'admissionDate';
                strdata_ = cell(size(data_));
                strdata_(~isnat(data_)) = cellstr(datestr(data_(~isnat(data_))));
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%s', c), strdata_, 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = {mean(data_); std(data_); median(data_); NaT; NaT; max(data_); min(data_)};
            end
            
            data_ = [obj.patients(idx).age]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Age (years)';
                demoSheet{row+1, col} = 'age';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.2f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).durationDiabetes]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Diabetes Duration (years)';
                demoSheet{row+1, col} = 'diabetesDuration';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.2f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).weight]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Weight (Kg)';
                demoSheet{row+1, col} = 'weight';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.2f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).height]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Height (cm)';
                demoSheet{row+1, col} = 'height';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.2f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).bmi]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'BMI (Kg/m2)';
                demoSheet{row+1, col} = 'bmi';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.2f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).female]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Female';
                demoSheet{row+1, col} = 'female';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%d', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).TDD]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Total Daily Dose (U)';
                demoSheet{row+1, col} = 'TDD';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.3f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).tddPerWeight]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'TDD per Weight (U/Kg)';
                demoSheet{row+1, col} = 'tddPerWeight';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.3f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).TDB]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'Total Daily Basal (U)';
                demoSheet{row+1, col} = 'TDB';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.3f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
            
            data_ = [obj.patients(idx).hba1c]';
            if nansum(data_) > 0
                col = col + 1;
                demoSheet{row, col} = 'HbA1c (%)';
                demoSheet{row+1, col} = 'hba1c';
                demoSheet(row+2:row+len+1, col) = cellfun(@(c) sprintf('%6.3f', c), num2cell(data_), 'UniformOutput', false);
                demoSheet(row+len+3:row+len+9, col) = cellstr(num2str([nanmean(data_); nanstd(data_); nanmedian(data_); prctile(data_, 25); prctile(data_, 75); max(data_); min(data_)], '%6.3f'));
            end
        end
        
        function summarySheets = getSummarySheet(obj, varargin)
            summarySheets = struct();
            
            dayOutcomes_ = true;
            nightOutcomes_ = true;
            
            arms_ = obj.arms;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'dayoutcomes'
                        dayOutcomes_ = varargin{nVar+1};
                    case 'nightoutcomes'
                        nightOutcomes_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = varargin{nVar+1};
                end
            end
            
            detailsTab = obj.details(varargin{:}, 'arms', arms_);
            outcomeTab = obj.outcomes(varargin{:}, 'arms', arms_);
            
            for idx = 1:numel(arms_)
                sheet = cell(0);
                
                % Adding details sheet
                gRow = 1;
                sheet{gRow, 1} = obj.arms{idx};
                for row = 1:size(detailsTab{idx}, 1)
                    sheet{row+gRow, 1} = detailsTab{idx}.Row{row};
                    for col = 1:size(detailsTab{idx}, 2)
                        sheet{gRow, col+1} = detailsTab{idx}.Properties.VariableNames{col};
                        if iscell(detailsTab{idx}{row, col})
                            if isnumeric(detailsTab{idx}{row, col}{1})
                                sheet{row+gRow, col+1} = sprintf('%6.3f', detailsTab{idx}{row, col}{1});
                            else
                                sheet{row+gRow, col+1} = detailsTab{idx}{row, col}{1};
                            end
                        else
                            sheet{row+gRow, col+1} = sprintf('%6.3f', detailsTab{idx}{row, col});
                        end
                    end
                end
                
                % Adding outcomes sheet
                gRow = row + gRow + 2;
                sheet{gRow, 1} = 'Overall';
                for row = 1:size(outcomeTab{idx}, 1)
                    sheet{row+gRow, 1} = outcomeTab{idx}.Row{row};
                    for col = 1:size(outcomeTab{idx}, 2)
                        sheet{gRow, col+1} = outcomeTab{idx}.Properties.VariableNames{col};
                        sheet{row+gRow, col+1} = sprintf('%6.3f', outcomeTab{idx}{row, col});
                    end
                end
                
                if dayOutcomes_ && days(obj.duration) >= 0.99
                    outcomeTabDay = obj.outcomes('interval', obj.dayInterval, varargin{:});
                    
                    gRow = row + gRow + 2;
                    sheet{gRow, 1} = sprintf('Daytime [%02d-%02d]', round(obj.dayInterval/60));
                    for row = 1:size(outcomeTabDay{idx}, 1)
                        sheet{row+gRow, 1} = outcomeTabDay{idx}.Row{row};
                        for col = 1:size(outcomeTabDay{idx}, 2)
                            sheet{gRow, col+1} = outcomeTabDay{idx}.Properties.VariableNames{col};
                            sheet{row+gRow, col+1} = sprintf('%6.3f', outcomeTabDay{idx}{row, col});
                        end
                    end
                end
                
                if nightOutcomes_ && days(obj.duration) >= 0.99
                    outcomeTabNight = obj.outcomes('interval', obj.nightInterval, varargin{:});
                    
                    gRow = row + gRow + 2;
                    sheet{gRow, 1} = sprintf('Nighttime [%02d-%02d]', round(obj.nightInterval/60));
                    for row = 1:size(outcomeTabNight{idx}, 1)
                        sheet{row+gRow, 1} = outcomeTabNight{idx}.Row{row};
                        for col = 1:size(outcomeTabNight{idx}, 2)
                            sheet{gRow, col+1} = outcomeTabNight{idx}.Properties.VariableNames{col};
                            sheet{row+gRow, col+1} = sprintf('%6.3f', outcomeTabNight{idx}{row, col});
                        end
                    end
                end
                
                summarySheets.(['Summary', arms_{idx}]) = sheet;
            end
        end
        
        function outSheet = getSummaryOutcomes(obj, varargin)
            arms_ = obj.arms;
            baseline_ = [];
            interval_ = [];
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'baseline'
                        baseline_ = varargin{nVar+1};
                    case 'interval'
                        interval_ = varargin{nVar+1};
                end
            end
            
            if ~isempty(baseline_)
                arms_ = [baseline_, arms_];
            end
            
            outSheet = cell(0);
            if ~isempty(interval_)
                outSheet{1, 1} = sprintf('[%02d-%02d]', round(interval_/60));
            end
            rowOffset = 2;
            outTab = obj.outcomes(varargin{:}, 'interval', interval_, 'arms', arms_, 'tablewithstats', false);
            if isempty(interval_)
                detailTab = obj.details(varargin{:}, 'arms', arms_);
            else
                detailTab = [];
            end
            if isempty(outTab)
                return;
            end
            rowNames = [];
            for n = 1:numel(outTab)
                rowNames = [rowNames, outTab{n}.Properties.VariableNames(contains(outTab{n}.Properties.VariableNames, 'P'))];
            end
            rowNames = unique(rowNames);
            rowNames = [rowNames, {'', 'Mean', 'SD', 'Median', 'IQR_25', 'IQR_75'}];
            
            colOffset = 2;
            if ~isempty(detailTab)
                for m = 1:numel(detailTab{1}.Row) % iterate through outcomes
                    for n = 1:numel(detailTab) % iterate through interventions
                        loadColumn(detailTab, n, m, colOffset);
                    end
                end
                colOffset = colOffset + numel(detailTab{1}.Row) * (numel(arms_) + 1);
            end
            for m = 1:numel(outTab{1}.Row) % iterate through outcomes
                for n = 1:numel(outTab) % iterate through interventions
                    loadColumn(outTab, n, m, colOffset);
                end
            end
            
            function loadColumn(tab, n, m, colOffset)
                if n == 1
                    outSheet{rowOffset-1, colOffset+(m - 1)*(numel(arms_) + 1)+n} = tab{n}.Properties.RowNames{m};
                else
                    outSheet{rowOffset-1, colOffset+(m - 1)*(numel(arms_) + 1)+n} = '##merge{-1}{0}';
                end
                outSheet{rowOffset, colOffset+(m - 1)*(numel(arms_) + 1)+n} = arms_{n};
                for p = 1:length(rowNames) % iterate through patients
                    if n == 1 && m == 1
                        outSheet{rowOffset+p, 1} = rowNames{p};
                    end
                    
                    switch rowNames{p}
                        case 'Mean'
                            outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = ...
                                sprintf('=IFERROR(ROUND(AVERAGE(%s%d:%s%d), 2), "")', ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+1, ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+length(rowNames)-6);
                        case 'SD'
                            outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = ...
                                sprintf('=IFERROR(ROUND(STDEV(%s%d:%s%d), 2), "")', ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+1, ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+length(rowNames)-6);
                        case 'Median'
                            outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = ...
                                sprintf('=IFERROR(ROUND(MEDIAN(%s%d:%s%d), 2), "")', ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+1, ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+length(rowNames)-6);
                        case 'IQR_25'
                            outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = ...
                                sprintf('=IFERROR(ROUND(QUARTILE(%s%d:%s%d, 1), 2), "")', ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+1, ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+length(rowNames)-6);
                        case 'IQR_75'
                            outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = ...
                                sprintf('=IFERROR(ROUND(QUARTILE(%s%d:%s%d, 3), 2), "")', ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+1, ...
                                obj.xlscol(colOffset+(m - 1)*(numel(arms_) + 1)+n), ...
                                rowOffset+length(rowNames)-6);
                        otherwise
                            if any(strcmp(tab{n}.Properties.VariableNames, rowNames{p}))
                                if isnumeric(tab{n}.(rowNames{p})(m))
                                    outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = round(tab{n}.(rowNames{p})(m), 2);
                                elseif isnumeric(tab{n}.(rowNames{p}){m})
                                    outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = round(tab{n}.(rowNames{p}){m}, 2);
                                else
                                    outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = tab{n}.(rowNames{p}){m};
                                end
                            else
                                outSheet{rowOffset+p, colOffset+(m - 1)*(numel(arms_) + 1)+n} = '';
                            end
                    end
                end
            end
        end
        
        function sheet = getCompareSheet(obj, varargin)
            arms_ = obj.arms;
            interval_ = [];
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'interval'
                        interval_ = varargin{nVar+1};
                end
            end
            
            prefix = join(arms_, '_');
            prefix = prefix{1};
            
            if isempty(interval_)
                suffix = '';
            elseif all(interval_ == obj.dayInterval)
                suffix = '_Day';
            elseif all(interval_ == obj.nightInterval)
                suffix = '_Night';
            end
            
            data = obj.getData('*');
            varNames = data(1).outcomes(varargin{:}, 'interval', interval_, 'arms', arms_).Properties.RowNames;
            dataCount = min(cellfun(@(c)sum(strcmp({data.name}, c)), arms_));
            
            sheet = cell(0);
            gRow = 1;
            for row = 1:length(varNames)
                gCol = 1;
                col = 0;
                
                if ~isempty(interval_)
                    sheet{gRow, gCol} = sprintf('[%02d-%02d]', round(interval_/60));
                end
                
                sheet{gRow+row, gCol} = varNames{row};
                for n = 1:length(arms_)
                    % generate mean
                    gCol = gCol + 1;
                    header = {'Mean', 'SD'};
                    for k = 1:length(header)
                        col = col + 1;
                        if row == 1
                            sheet{gRow, gCol+col} = [arms_{n}, '_', header{k}];
                        end
                        sheet{gRow+row, gCol+col} = sprintf('=OFFSET(INDEX(Outcomes%s!$%d:$%d,MATCH(A%d,Outcomes%s!$1:$1,0)), 0, %d)', ...
                            suffix, ...
                            dataCount+3+k, ...
                            dataCount+3+k, ...
                            gRow+row, ...
                            suffix, ...
                            n-1);
                    end
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = [arms_{n}, '_Mean(SD)'];
                    end
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2) & " (" & FIXED(%s%d,2) & ")"', ...
                        obj.xlscol(gCol+col-2), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-1), ...
                        gRow+row);
                    
                    % generate median
                    gCol = gCol + 1;
                    header = {'Median', 'IQR_25', 'IQR_75'};
                    for k = 1:length(header)
                        col = col + 1;
                        if row == 1
                            sheet{gRow, gCol+col} = [arms_{n}, '_', header{k}];
                        end
                        sheet{gRow+row, gCol+col} = sprintf('=OFFSET(INDEX(Outcomes%s!$%d:$%d,MATCH(A%d,Outcomes%s!$1:$1,0)), 0, %d)', ...
                            suffix, ...
                            dataCount+5+k, ...
                            dataCount+5+k, ...
                            gRow+row, ...
                            suffix, ...
                            n-1);
                    end
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = [arms_{n}, '_Median[IQR]'];
                    end
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2) & " [" & FIXED(%s%d,2) & "–" & FIXED(%s%d,2) & "]"', ...
                        obj.xlscol(gCol+col-3), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-2), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-1), ...
                        gRow+row);
                end
                
                % generate diff
                gCol = gCol + 1;
                header = {[prefix, '_Diff'], [prefix, '_CI_2_5'], [prefix, '_CI_97_5']};
                for k = 1:length(header)
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = header{k};
                    end
                    sheet{gRow+row, gCol+col} = sprintf('=INDEX(LME%s!$%d:$%d,MATCH(A%d,LME%s!$1:$1,0))', ...
                        suffix, ...
                        numel(data)+21+k, ...
                        numel(data)+21+k, ...
                        gRow+row, ...
                        suffix);
                end
                col = col + 1;
                if row == 1
                    sheet{gRow, gCol+col} = [prefix, '_Diff(CI)'];
                end
                sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2) & " (" & FIXED(%s%d,2) & " to " & FIXED(%s%d,2) & ")"', ...
                    obj.xlscol(gCol+col-3), ...
                    gRow+row, ...
                    obj.xlscol(gCol+col-2), ...
                    gRow+row, ...
                    obj.xlscol(gCol+col-1), ...
                    gRow+row);
                
                % generate p-value
                gCol = gCol + 1;
                header = {'Normality', 'P_Value'};
                for k = 1:length(header)
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = header{k};
                    end
                    sheet{gRow+row, gCol+col} = sprintf('=INDEX(LME%s!$%d:$%d,MATCH(A%d,LME%s!$1:$1,0))', ...
                        suffix, ...
                        numel(data)+19+k, ...
                        numel(data)+19+k, ...
                        gRow+row, ...
                        suffix);
                end
            end
        end
        
        function sheet = getOutcomeSheets(obj, varargin)
            arms_ = obj.arms;
            interval_ = [];
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'interval'
                        interval_ = varargin{nVar+1};
                end
            end
            
            prefix = join(arms_, '_');
            prefix = prefix{1};
            
            outcomeMean = obj.compare(varargin{:}, 'type', 'mean', 'arms', arms_);
            outcomeMedian = obj.compare(varargin{:}, 'type', 'median', 'arms', arms_);
            
            varNames = outcomeMean.Properties.RowNames;
            sheet = cell(0);
            gRow = 1;
            for row = 1:length(varNames)
                gCol = 1;
                col = 0;
                
                if ~isempty(interval_)
                    sheet{gRow, gCol} = sprintf('[%02d-%02d]', round(interval_/60));
                end
                
                sheet{gRow+row, gCol} = varNames{row};
                for n = 1:length(arms_)
                    % generate mean
                    gCol = gCol + 1;
                    header = {'Mean', 'SD'};
                    for k = 1:length(header)
                        col = col + 1;
                        if row == 1
                            sheet{gRow, gCol+col} = [arms_{n}, '_', header{k}];
                        end
                        sheet{gRow+row, gCol+col} = outcomeMean.([arms_{n}, '_', header{k}])(row);
                    end
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = [arms_{n}, '_Mean(SD)'];
                    end
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2) & " (" & FIXED(%s%d,2) & ")"', ...
                        obj.xlscol(gCol+col-2), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-1), ...
                        gRow+row);
                end
                
                % Normality p value
                gCol = gCol + 1;
                col = col + 1;
                if row == 1
                    sheet{gRow, gCol+col} = [prefix, '_Normality_P_value'];
                end
                sheet{gRow+row, gCol+col} = outcomeMean.([prefix, '_Normality'])(row);
                
                for n = 1:length(arms_)
                    % generate median
                    gCol = gCol + 1;
                    header = {'Median', 'IQR_25', 'IQR_75'};
                    for k = 1:length(header)
                        col = col + 1;
                        if row == 1
                            sheet{gRow, gCol+col} = [arms_{n}, '_', header{k}];
                        end
                        sheet{gRow+row, gCol+col} = outcomeMedian.([arms_{n}, '_', header{k}])(row);
                    end
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = [arms_{n}, '_Median[IQR]'];
                    end
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2) & " [" & FIXED(%s%d,2) & "–" & FIXED(%s%d,2) & "]"', ...
                        obj.xlscol(gCol+col-3), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-2), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-1), ...
                        gRow+row);
                end
                
                % Mean Diff
                gCol = gCol + 1;
                header = {[prefix, '_Diff'], [prefix, '_CI_2_5'], [prefix, '_CI_97_5'], 'P_Value'};
                for k = 1:length(header)
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = ['Mean_', header{k}];
                    end
                    sheet{gRow+row, gCol+col} = outcomeMean.(header{k})(row);
                end
                
                % Median Diff
                gCol = gCol + 1;
                header = {[prefix, '_Diff'], [prefix, '_CI_2_5'], [prefix, '_CI_97_5'], 'P_Value'};
                for k = 1:length(header)
                    col = col + 1;
                    if row == 1
                        sheet{gRow, gCol+col} = ['Median_', header{k}];
                    end
                    sheet{gRow+row, gCol+col} = outcomeMedian.(header{k})(row);
                end
                
                % final CI/p-value
                gCol = gCol + 1;
                col = col + 1;
                if row == 1
                    sheet{gRow, gCol+col} = [prefix, '_Normality'];
                end
                isNormal = outcomeMean.([prefix, '_Normality'])(row) > 0.05;
                sheet{gRow+row, gCol+col} = sprintf('=IF(ISBLANK(%s%d), "", IF(%s%d>0.05,"NORMAL", "SKEWED"))', ...
                    obj.xlscol(gCol+col-22), ...
                    gRow+row, ...
                    obj.xlscol(gCol+col-22), ...
                    gRow+row);
                col = col + 1;
                if row == 1
                    sheet{gRow, gCol+col} = [prefix, '_Diff(CI)'];
                end
                if isNormal
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2) & " (" & FIXED(%s%d,2) & " to " & FIXED(%s%d,2) & ")"', ...
                        obj.xlscol(gCol+col-11), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-10), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-9), ...
                        gRow+row);
                else
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2) & " (" & FIXED(%s%d,2) & " to " & FIXED(%s%d,2) & ")"', ...
                        obj.xlscol(gCol+col-6), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-5), ...
                        gRow+row, ...
                        obj.xlscol(gCol+col-4), ...
                        gRow+row);
                end
                col = col + 1;
                if row == 1
                    sheet{gRow, gCol+col} = [prefix, '_P_Value'];
                end
                if isNormal
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2)', ...
                        obj.xlscol(gCol+col-9), ...
                        gRow+row);
                else
                    sheet{gRow+row, gCol+col} = sprintf('=FIXED(%s%d,2)', ...
                        obj.xlscol(gCol+col-4), ...
                        gRow+row);
                end
            end
        end
        
        function sheet = getLMESheet(obj, varargin)
            arms_ = obj.arms;
            interval_ = [];
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'interval'
                        interval_ = varargin{nVar+1};
                end
            end
            
            tab = obj.lme(varargin{:}, 'interval', interval_, 'arms', arms_);
            headerLen = 4;
            
            data = obj.getData('*');
            varNames = data(1).outcomes(varargin{:}, 'interval', interval_, 'arms', arms_).Properties.RowNames;
            
            sheet = cell(0);
            for row = 1:size(tab, 1)
                if contains(tab.Row{row}, 'Row')
                    gRow = 1;
                else
                    gRow = 2;
                end
                for col = 1:size(tab, 2)
                    if row == 1
                        if col <= headerLen
                            sheet{1, col} = tab.Properties.VariableNames{col};
                        else
                            sheet{1, col} = varNames{col-headerLen};
                        end
                    end
                    if iscell(tab{row, col})
                        sheet{row+1, col} = tab{row, col}{1};
                    elseif iscategorical(tab{row, col})
                        if isundefined(tab{row, col})
                            if col == 1
                                sheet{row+gRow, col} = tab.Row{row};
                            else
                                sheet{row+gRow, col} = '##merge{-1}{0}';
                            end
                        else
                            sheet{row+gRow, col} = char(tab{row, col});
                        end
                    else
                        sheet{row+gRow, col} = tab{row, col};
                    end
                end
            end
            row = size(tab, 1) + 1;
            
            offsetNormality = find(strcmp(tab.Properties.RowNames, 'Shapiro Normality Test'));
            offsetSquence = find(strcmp(tab.Properties.RowNames, 'P-Value Sequence Effects'));
            offsetPeriod = find(strcmp(tab.Properties.RowNames, 'P-Value Period Effects'));
            offsetPValueNormal = find(strcmp(tab.Properties.RowNames, 'P-Value Intervention Effects'));
            offsetPValueSkewed = find(strcmp(tab.Properties.RowNames, 'P-Value Wilcoxon Test'));
            offsetDiffNormal = find(strcmp(tab.Properties.RowNames, 'Fixed Effects'));
            offsetDiffSkewed = find(strcmp(tab.Properties.RowNames, 'Fixed Effects (Hodges-Lehmann)'));
            offsetCI1Normal = find(strcmp(tab.Properties.RowNames, '2.5% Quantile'));
            offsetCI1Skewed = find(strcmp(tab.Properties.RowNames, '2.5% Quantile (Hodges-Lehmann)'));
            offsetCI2Normal = find(strcmp(tab.Properties.RowNames, '97.5% Quantile'));
            offsetCI2Skewed = find(strcmp(tab.Properties.RowNames, '97.5% Quantile (Hodges-Lehmann)'));
            
            gRow = gRow + 1;
            for col = 1:size(tab, 2)
                if col == 1
                    sheet{row+gRow, col} = '__Sequence order effect__';
                elseif col <= headerLen
                    sheet{row+gRow, col} = '##merge{-1}{0}';
                else
                    sheet{row+gRow, col} = sprintf('=IF(OR(ISBLANK(%s%d), ISBLANK(%s%d)), "", IF(AND(%s%d>0.05, %s%d>0.05),"NO", "YES"))', ...
                        obj.xlscol(col), ...
                        row+offsetSquence-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetPeriod-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetSquence-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetPeriod-height(tab)+1);
                end
            end
            
            gRow = gRow + 1;
            for col = 1:size(tab, 2)
                if col == 1
                    sheet{row+gRow, col} = '__Normality__';
                elseif col <= headerLen
                    sheet{row+gRow, col} = '##merge{-1}{0}';
                else
                    sheet{row+gRow, col} = sprintf('=IF(ISBLANK(%s%d), "", IF(%s%d>0.05,"NORMAL", "SKEWED"))', ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1);
                end
            end
            
            gRow = gRow + 1;
            for col = 1:size(tab, 2)
                if col == 1
                    sheet{row+gRow, col} = '__P-Value__';
                elseif col <= headerLen
                    sheet{row+gRow, col} = '##merge{-1}{0}';
                else
                    sheet{row+gRow, col} = sprintf('=IF(ISBLANK(%s%d), "", IF(%s%d>0.05,%s%d,%s%d))', ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetPValueNormal-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetPValueSkewed-height(tab)+1);
                end
            end
            
            gRow = gRow + 1;
            for col = 1:size(tab, 2)
                if col == 1
                    sheet{row+gRow, col} = '__Fixed Effects__';
                elseif col <= headerLen
                    sheet{row+gRow, col} = '##merge{-1}{0}';
                else
                    sheet{row+gRow, col} = sprintf('=IF(ISBLANK(%s%d), "", IF(%s%d>0.05,%s%d,%s%d))', ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetDiffNormal-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetDiffSkewed-height(tab)+1);
                end
            end
            
            gRow = gRow + 1;
            for col = 1:size(tab, 2)
                if col == 1
                    sheet{row+gRow, col} = '__2.5% Quantile__';
                elseif col <= headerLen
                    sheet{row+gRow, col} = '##merge{-1}{0}';
                else
                    sheet{row+gRow, col} = sprintf('=IF(ISBLANK(%s%d), "", IF(%s%d>0.05,%s%d,%s%d))', ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetCI1Normal-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetCI1Skewed-height(tab)+1);
                end
            end
            
            gRow = gRow + 1;
            for col = 1:size(tab, 2)
                if col == 1
                    sheet{row+gRow, col} = '__97.5% Quantile__';
                elseif col <= headerLen
                    sheet{row+gRow, col} = '##merge{-1}{0}';
                else
                    sheet{row+gRow, col} = sprintf('=IF(ISBLANK(%s%d), "", IF(%s%d>0.05,%s%d,%s%d))', ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetNormality-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetCI2Normal-height(tab)+1, ...
                        obj.xlscol(col), ...
                        row+offsetCI2Skewed-height(tab)+1);
                end
            end
        end
        
        function sheets = toSheet(obj, varargin)
            demographics_ = true;
            individual_ = false;
            arms_ = obj.arms;
            type_ = obj.type;
            baseline_ = '';
            lmeSheet_ = false;
            compareSheet_ = false;
            outcomesSheet_ = true;
            dayNightSheets_ = true;
            summarySheet_ = false;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'individual'
                        individual_ = varargin{nVar+1};
                    case {'outcomes', 'outcome', 'outcomessheet', 'outcomesheet'}
                        outcomesSheet_ = varargin{nVar+1};
                    case {'summary', 'summarysheet'}
                        summarySheet_ = varargin{nVar+1};
                    case {'demo', 'demographic', 'demographics'}
                        demographics_ = varargin{nVar+1};
                    case {'compare', 'comparesheet'}
                        compareSheet_ = varargin{nVar+1};
                    case {'lme', 'lmesheet', 'linearmixedeffects', 'linear-mixed-effects'}
                        lmeSheet_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'type'
                        type_ = varargin{nVar+1};
                    case 'baseline'
                        baseline_ = varargin{nVar+1};
                    case {'day', 'night', 'daynight', 'daysheet', 'nightsheet', 'daynightsheet'}
                        dayNightSheets_ = varargin{nVar+1};
                end
            end
            
            if compareSheet_
                lmeSheet_ = true;
                outcomesSheet_ = true;
            end
            
            sheets = struct();
            
            % Demographics sheet
            if demographics_
                sheets.Demographics = obj.getDemographSheet;
            end
            
            % summary of outcomes
            if outcomesSheet_
                sheets.Outcomes = obj.getSummaryOutcomes(varargin{:});
                if dayNightSheets_ && days(obj.duration) >= 0.99
                    sheets.Outcomes_Day = obj.getSummaryOutcomes(varargin{:}, 'interval', obj.dayInterval);
                end
                if dayNightSheets_ && days(obj.duration) >= 0.99
                    sheets.Outcomes_Night = obj.getSummaryOutcomes(varargin{:}, 'interval', obj.nightInterval);
                end
            end
            
            if strcmp(type_, 'matched')
                % Compare sheets
                if compareSheet_ && numel(arms_) >= 2
                    couples = nchoosek(1:numel(arms_), 2);
                    for i = 1:size(couples, 1)
                        sheets.(['Compare_', arms_{couples(i, 1)}, '_', arms_{couples(i, 2)}]) = obj.getOutcomeSheets( ...
                            varargin{:}, ...
                            'arms', arms_(couples(i, :)));
                        if dayNightSheets_ && days(obj.duration) >= 0.99
                            sheets.(['Compare_', arms_{couples(i, 1)}, '_', arms_{couples(i, 2)}, '_Day']) = obj.getOutcomeSheets( ...
                                varargin{:}, ...
                                'interval', obj.dayInterval, ...
                                'arms', arms_(couples(i, :)));
                        end
                        if dayNightSheets_ && days(obj.duration) >= 0.99
                            sheets.(['Compare_', arms_{couples(i, 1)}, '_', arms_{couples(i, 2)}, '_Night']) = obj.getOutcomeSheets( ...
                                varargin{:}, ...
                                'interval', obj.nightInterval, ...
                                'arms', arms_(couples(i, :)));
                        end
                    end
                end
            elseif strcmp(type_, 'crossover')
                % LME sheets
                if lmeSheet_ && numel(arms_) == 2
                    sheets.LME = obj.getLMESheet(varargin{:});
                    if dayNightSheets_ && days(obj.duration) >= 0.99
                        sheets.LME_Day = obj.getLMESheet(varargin{:}, 'interval', obj.dayInterval);
                    end
                    if dayNightSheets_ && days(obj.duration) >= 0.99
                        sheets.LME_Night = obj.getLMESheet(varargin{:}, 'interval', obj.nightInterval);
                    end
                end
                
                % Compare sheets
                if compareSheet_ % TODO multiple arms ?
                    sheets.Compare = obj.getCompareSheet(varargin{:});
                    if dayNightSheets_ && days(obj.duration) >= 0.99
                        sheets.Compare_Day = obj.getCompareSheet(varargin{:}, 'interval', obj.dayInterval);
                    end
                    if dayNightSheets_ && days(obj.duration) >= 0.99
                        sheets.Compare_Night = obj.getCompareSheet(varargin{:}, 'interval', obj.nightInterval);
                    end
                end
            else % TODO parallel ?
            end
            
            % Summary sheets (optional)
            if summarySheet_
                if ~isempty(obj.patients)
                    summarySheets = obj.getSummarySheet(varargin{:});
                    for fn = fieldnames(summarySheets)'
                        sheets.(fn{1}) = summarySheets.(fn{1});
                    end
                    if ~isempty(baseline_)
                        summarySheets = obj.getSummarySheet(varargin{:}, 'arms', baseline_);
                        for fn = fieldnames(summarySheets)'
                            sheets.(fn{1}) = summarySheets.(fn{1});
                        end
                    end
                end
            end
            
            % Patients sheet
            if individual_
                for p = 1:numel(obj.patients)
                    sheets.(obj.patients(p).name) = obj.patients(p).toSheet(varargin{:});
                end
            end
        end
        
        function [fig_, ax_, lgd_, data_] = toFigure(obj, varargin)
            inScreen_ = true;
            average_ = false;
            arms_ = obj.arms;
            title_ = '';
            days_ = [];
            patients_ = [obj.patients.ID];
            stats_ = true;
            axis_ = [];
            summaryFct_ = @MAPUtils.plotSummary;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'inscreen'
                        inScreen_ = varargin{nVar+1};
                    case 'average'
                        average_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'title'
                        title_ = varargin{nVar+1};
                    case 'days'
                        days_ = varargin{nVar+1};
                    case 'patients'
                        patients_ = varargin{nVar+1};
                    case 'stats'
                        stats_ = varargin{nVar+1};
                    case 'axis'
                        axis_ = varargin{nVar+1};
                    case {'summaryfun', 'summaryfct', 'summaryfunction'}
                        summaryFct_ = varargin{nVar+1};
                end
            end
            
            data_ = cell(numel(arms_), 1);
            for n = 1:numel(arms_)
                data_{n} = [];
                for p = 1:length(obj.patients)
                    if ~any(patients_ == obj.patients(p).ID)
                        continue;
                    end
                    if ~any(strcmp({obj.patients(p).data.name}, arms_{n}))
                        continue;
                    end
                    pStruct = obj.patients(p).toStruct( ...
                        'outcomes', true, ...
                        'days', days_, ...
                        'average', average_, ...
                        'arms', arms_{n});
                    data_{n} = [data_{n}, pStruct.data];
                end
                if any(mode([data_{n}.duration])-[data_{n}.duration] > minutes(1))
                    warning('[MAPStudy][toFigure] Data duration is different between interventions in arm %s ... can''t show summary.', arms_{n});
                    fig_ = gobjects(0);
                    return;
                end
            end
            
            if numel(unique(round(minutes([data_{1}.duration])))) > 1
                warning('[MAPStudy][toFigure] Data duration is different between patients ... can''t show summary.');
                fig_ = gobjects(0);
                return;
            end
            duration_ = minutes(data_{1}(1).duration);
            
            if ~isempty(axis_)
                fig_ = [];
                ax_ = axis_;
                stats_ = false;
            else
                doNotRefresh = false;
                if inScreen_
                    if ishandle(1111)
                        doNotRefresh = true;
                    end
                    fig_ = figure(1111);
                else
                    fig_ = figure('Visible', 'Off');
                end
                clf(fig_);
                if ~doNotRefresh
                    set(fig_, 'name', sprintf('MAPStudy::toFigure::%s', obj.name), ...
                        'numbertitle', 'off', ...
                        'units', 'normalized', ...
                        'outerposition', [0, 0, 1, 1]);
                end
                
                if stats_
                    axisPositions = [0.04, 0.15, 0.75, 0.75];
                else
                    axisPositions = [0.09, 0.13, 0.8, 0.8];
                end
                
                ax_ = subplot('Position', axisPositions);
                ax_.Position = axisPositions;
            end
            
            if isempty(title_)
                title_ = sprintf(['%s Summary of ', repmat('%s, ', 1, numel(arms_)), '~ %6.3f hours'], obj.name, arms_{:}, hours(minutes(duration_)));
            end
            [ax_, lgd_] = summaryFct_(ax_, data_, ...
                varargin{:}, ...
                'duration', duration_, ...
                'title', title_);
            
            if stats_
                MAPUtils.plotStatsSummary(fig_, data_)
            end
        end
        
        function toPNG(obj, path, varargin)
            if nargin < 2
                error('[MAPStudy][toPNG] Needs at least one argument.');
            end
            
            inScreen_ = false;
            inFolder_ = true;
            daily_ = false;
            arms_ = obj.arms;
            nested_ = obj.nested;
            patients_ = [obj.patients.ID];
            summary_ = true;
            individual_ = true;
            average_ = false;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'inscreen'
                        inScreen_ = varargin{nVar+1};
                    case 'infolder'
                        inFolder_ = varargin{nVar+1};
                    case 'nested'
                        nested_ = varargin{nVar+1};
                    case 'daily'
                        daily_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'patients'
                        patients_ = varargin{nVar+1};
                    case 'summary'
                        summary_ = varargin{nVar+1};
                    case 'individual'
                        individual_ = varargin{nVar+1};
                    case 'average'
                        average_ = varargin{nVar+1};
                end
            end
            
            if individual_
                if nested_ > 0
                    if contains(path, obj.name)
                        pPath = [path, filesep];
                    else
                        pPath = [path, obj.name, filesep];
                    end
                else
                    if contains(path, obj.name)
                        pPath = path;
                    else
                        pPath = [path, obj.name, '_'];
                    end
                end
                for p = 1:numel(obj.patients)
                    if ~any(patients_ == obj.patients(p).ID)
                        continue;
                    end
                    
                    obj.patients(p).toPNG(pPath, ...
                        varargin{:}, ...
                        'inscreen', inScreen_, ...
                        'infolder', inFolder_, ...
                        'daily', daily_, ...
                        'average', average_, ...
                        'summary', summary_, ...
                        'arms', arms_, ...
                        'nested', nested_-1);
                end
            end
            
            if summary_
                fig = obj.toFigure( ...
                    varargin{:}, ...
                    'inscreen', inScreen_);
                
                if ~isempty(fig)
                    if path(end) == '\' || path(end) == '/'
                        if length([arms_{:}]) < 10
                            dPath = [path, obj.name, 'Summary_', [arms_{:}]];
                        else
                            dPath = [path, obj.name, 'Summary'];
                        end
                    else
                        dPath = path;
                    end
                    obj.ensureFolder(dPath);
                    print(fig, [dPath, '.png'], '-dpng');
                    
                    if ~inScreen_
                        close(fig);
                    end
                end
            end
        end
        
        function detailsTab = details(obj, varargin)
            arms_ = obj.arms;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'arms'
                        arms_ = varargin{nVar+1};
                end
            end
            
            rowNames_ = '';
            detailsStruct = struct([]);
            for p = 1:numel(obj.patients)
                out_ = obj.patients(p).details(varargin{:});
                if isempty(rowNames_) && ~isempty(out_.Properties.RowNames)
                    rowNames_ = out_.Properties.RowNames;
                end
                for idx = 1:length(arms_)
                    if any(strcmpi(out_.Properties.VariableNames, arms_{idx}))
                        detailsStruct(idx).(obj.patients(p).name) = out_.(arms_{idx});
                    end
                end
            end
            
            detailsTab = cell(1, numel(detailsStruct));
            for idx = 1:numel(detailsStruct)
                fn_ = fieldnames(detailsStruct(idx));
                detailsTab{idx} = struct2table(rmfield(detailsStruct(idx), fn_(structfun(@isempty, detailsStruct(idx)))));
                
                detailsTab{idx} = detailsTab{idx}(:, sort(detailsTab{idx}.Properties.VariableNames));
                detailsTab{idx}.Properties.RowNames = rowNames_;
            end
        end
        
        function outcomeTab = outcomes(obj, varargin)
            arms_ = obj.arms;
            patientID_ = [obj.patients.ID];
            tableWithStats_ = true;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'patients'
                        patientID_ = varargin{nVar+1};
                    case 'tablewithstats'
                        tableWithStats_ = varargin{nVar+1};
                end
            end
            
            RowNames = [];
            outcomeStruct = struct([]);
            for p = 1:numel(obj.patients)
                if all(obj.patients(p).ID ~= patientID_)
                    continue;
                end
                out_ = obj.patients(p).outcomes(varargin{:});
                if isempty(out_)
                    continue;
                end
                if isempty(RowNames)
                    RowNames = out_.Properties.RowNames;
                end
                for idx = 1:numel(out_.Properties.VariableNames)
                    outcomeStruct(strcmpi(arms_, out_.Properties.VariableNames{idx})).(obj.patients(p).name) = ...
                        out_.(out_.Properties.VariableNames{idx});
                end
            end
            
            outcomeTab = cell(1, numel(outcomeStruct));
            for idx = 1:numel(outcomeStruct)
                fn_ = fieldnames(outcomeStruct(idx));
                outcomeTab{idx} = struct2table(rmfield(outcomeStruct(idx), fn_(structfun(@isempty, outcomeStruct(idx)))));
                outcomeTab{idx} = outcomeTab{idx}(:, sort(outcomeTab{idx}.Properties.VariableNames));
                outcomeTab{idx}.Properties.RowNames = RowNames;
                if tableWithStats_
                    n_ = size(outcomeTab{idx}.Variables, 2);
                    outcomeTab{idx}.Mean = nanmean(outcomeTab{idx}{:, 1:n_}, 2);
                    outcomeTab{idx}.SD = nanstd(outcomeTab{idx}{:, 1:n_}, 0, 2);
                    outcomeTab{idx}.Median = nanmedian(outcomeTab{idx}{:, 1:n_}, 2);
                    outcomeTab{idx}.IQR_25 = prctile(outcomeTab{idx}{:, 1:n_}, 25, 2);
                    outcomeTab{idx}.IQR_75 = prctile(outcomeTab{idx}{:, 1:n_}, 75, 2);
                end
            end
        end
        
        function outcomeTab = compare(obj, varargin)
            type_ = 'mean';
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'type'
                        type_ = varargin{nVar+1};
                end
            end
            
            if strcmpi(type_, 'mean')
                outcomeTab = compareMean(obj, varargin{:});
            else
                outcomeTab = compareMedian(obj, varargin{:});
            end
        end
    end
    
    %     methods (Access = private)
    %         function out_ = getField(obj, fn)
    %             out_ = [];
    %             timeStamp_ = [];
    %             for p = 1:numel(obj.patients)
    %                 timeStamp_ = [timeStamp_; obj.patients(p).timeStamp];
    %                 if nansum(obj.patients(p).(fn)) > 0
    %                     out_ = [out_; obj.patients(p).(fn)];
    %                 else
    %                     out_ = [out_; nan(size(obj.patients(p).timeStamp))];
    %                 end
    %             end
    %             [sortedTimeStamp, iSort] = sort(timeStamp_);
    %             [~, iUnique] = unique(round(sortedTimeStamp*60*24));
    %             n = length(timeStamp_);
    %             out_ = out_(iSort(any(iUnique == 1:n)));
    %         end
    %     end
end
