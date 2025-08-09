classdef MAPPatient < MAPUtils
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = public)
        nested = 0 % if 1 a folder with the data 'name' will be created.
    end
    
    properties (GetAccess = public, SetAccess = private)
        data = MAPData.empty % List of data for the patient
    end
    
    properties (GetAccess = public, SetAccess = ?MAPStudy)
        studyName % The study where this patient participated
    end
    
    properties (GetAccess = public, SetAccess = immutable)
        ID % Patient ID
    end
    
    properties (GetAccess = public, SetAccess = public)
        TDD = NaN % Total daily dose
        TDB = NaN
        weight = NaN % Weight in (Kg)
        height = NaN % in (cm)
        hba1c = NaN
        gender = '' % 'M' or 'F'
        diagnosticYear = NaN % a double for the year of diagnotistics (not a duration)
        birthDate = NaT % datetime object
        admissionDate = NaT % datetime object for admission date
        studyStartDate = NaT % datetime object for study start date
        sensorType = '' % Name/type of glucose sensor that patient use
        sensorDuration = NaT - NaT % Duration of glucose sensor usage
        therapyType = '' % this is a string, example: 'pump', 'mdi', 'mdi_carb_counting', 'mdi_fix_dose'
        therapyParam = struct() % patient default parameter
    end
    
    properties (Dependent)
        units
        isEmpty
        stepTime
        name % Patient name
        tddPerWeight % Total daily dose per Weight
        bmi
        female % is this patient a female (True or false)
        durationDiabetes
        age % This is the age of the patient at the admission date
        ageToday
        arms % Name of Arms
        armsSequence
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
        outcomeEndTime % if set outcome calculation end at this time, should be smaller than startTime
        softMerge % if true, data is only added if was NaN or empty previously.
    end
    
    methods
        function set.studyName(obj, study_)
            obj.studyName = study_;
            if isprop(obj, 'data')
                for dIdx = 1:length(obj.data)
                    obj.data(dIdx).studyName = study_;
                end
            end
        end
        
        function flag = get.isEmpty(obj)
            flag = true;
            if nansum(obj.sensorGlucose) > 0 || ...
                    nansum(obj.bloodGlucose) > 0 || ...
                    nansum(obj.basalInsulin) > 0 || ...
                    nansum(obj.bolusInsulin) > 0 || ...
                    nansum(obj.basalInjection) > 0
                flag = false;
            end
        end
        
        function dt = get.stepTime(obj)
            dt = [];
            if ~isempty(obj.data)
                dt = mode([obj.data.stepTime]);
                if any([obj.data.stepTime] ~= dt)
                    error('Data must have the same step time!');
                end
            end
        end
        
        function set.name(obj, name_)
            obj.ID = str2double(regexp(name_, '\d+.?', 'match'));
            if isprop(obj, 'data')
                for dIdx = 1:length(obj.data)
                    obj.data(dIdx).patientName = obj.name;
                    obj.data(dIdx).patientHandle = obj;
                end
            end
        end
        
        function n = get.name(obj)
            n = num2str(obj.ID, 'P%03d');
        end
        
        function val = get.tddPerWeight(obj)
            val = obj.TDD / obj.weight;
        end
        
        function val = get.bmi(obj)
            val = 1e4 * obj.weight / (obj.height^2);
        end
        
        function val = get.female(obj)
            val = NaN;
            if ~isempty(obj.gender)
                val = strcmpi(obj.gender, 'F');
            end
        end
        
        function val = get.durationDiabetes(obj)
            val = year(obj.admissionDate) - obj.diagnosticYear;
        end
        
        function set.age(obj, duration)
            if isduration(duration)
                obj.birthDate = obj.admissionDate - duration;
            else
                obj.birthDate = obj.admissionDate - years(duration);
            end
        end
        
        function val = get.age(obj)
            val = years(obj.admissionDate-obj.birthDate);
        end
        
        function val = get.ageToday(obj)
            val = years(datetime(now, 'ConvertFrom', 'datenum')-obj.birthDate);
        end
        
        function val = get.arms(obj)
            val = {obj.data.name};
            val = val([obj.data.outcome] == 1);
        end
        
        function val = get.armsSequence(obj)
            [~, indices] = sort([obj.data.startDate]);
            indices([obj.data.outcome] == 0) = [];
            val = join(obj.arms(indices), '_');
            val = val{1};
        end
        
        function int = get.useGlucoseInterp(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).useGlucoseInterp;
            end
        end
        
        function set.useGlucoseInterp(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).useGlucoseInterp = int_;
            end
        end
        
        function int = get.useBloodGlucoseInGlucoseInterp(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).useBloodGlucoseInGlucoseInterp;
            end
        end
        
        function set.useBloodGlucoseInGlucoseInterp(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).useBloodGlucoseInGlucoseInterp = int_;
            end
        end
        
        function int = get.glucoseInterpHolesSize(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).glucoseInterpHolesSize;
            end
        end
        
        function set.glucoseInterpHolesSize(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).glucoseInterpHolesSize = int_;
            end
        end
        
        function int = get.glucoseInterpHolesSizeAfterMealsOrBoluses(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).glucoseInterpHolesSizeAfterMealsOrBoluses;
            end
        end
        
        function set.glucoseInterpHolesSizeAfterMealsOrBoluses(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).glucoseInterpHolesSizeAfterMealsOrBoluses = int_;
            end
        end
        
        function int = get.thresholdHypoglycemiaEvent(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).thresholdHypoglycemiaEvent;
            end
        end
        
        function set.thresholdHypoglycemiaEvent(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).thresholdHypoglycemiaEvent = int_;
            end
        end
        
        function int = get.durationHypoglycemiaEvent(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).durationHypoglycemiaEvent;
            end
        end
        
        function set.durationHypoglycemiaEvent(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).durationHypoglycemiaEvent = int_;
            end
        end
        
        function int = get.nightInterval(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).nightInterval;
            end
        end
        
        function set.nightInterval(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).nightInterval = int_;
            end
        end
        
        function int = get.dayInterval(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).dayInterval;
            end
        end
        
        function set.dayInterval(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).dayInterval = int_;
            end
        end
        
        function int = get.breakfastInterval(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).breakfastInterval;
            end
        end
        
        function set.breakfastInterval(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).breakfastInterval = int_;
            end
        end
        
        function int = get.lunchInterval(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).lunchInterval;
            end
        end
        
        function set.lunchInterval(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).lunchInterval = int_;
            end
        end
        
        function int = get.dinnerInterval(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).dinnerInterval;
            end
        end
        
        function set.dinnerInterval(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).dinnerInterval = int_;
            end
        end
        
        function int = get.snacksThreshold(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).snacksThreshold;
            end
        end
        
        function set.snacksThreshold(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).snacksThreshold = int_;
            end
        end
        
        function int = get.onlyCountBolusedMeals(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).onlyCountBolusedMeals;
            end
        end
        
        function set.onlyCountBolusedMeals(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).onlyCountBolusedMeals = int_;
            end
        end
        
        function int = get.onlyCountBiggerMeal(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).onlyCountBiggerMeal;
            end
        end
        
        function set.onlyCountBiggerMeal(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).onlyCountBiggerMeal = int_;
            end
        end
        
        function int = get.mealInfoSnapshot(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).mealInfoSnapshot;
            end
        end
        
        function set.mealInfoSnapshot(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).mealInfoSnapshot = int_;
            end
        end
        
        function int = get.computeIOB(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).computeIOB;
            end
        end
        
        function set.computeIOB(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).computeIOB = int_;
            end
        end
        
        function int = get.outcomeStartTime(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(strcmp({obj.data.name}, obj.arms{1})).outcomeStartTime;
            end
        end
        
        function set.outcomeStartTime(obj, int_)
            for d = 1:numel(obj.data)
                if any(strcmp(obj.arms, obj.data(d).name))
                    obj.data(d).outcomeStartTime = int_;
                end
            end
        end
        
        function int = get.outcomeEndTime(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(strcmp({obj.data.name}, obj.arms{1})).outcomeEndTime;
            end
        end
        
        function set.outcomeEndTime(obj, int_)
            for d = 1:numel(obj.data)
                if any(strcmp(obj.arms, obj.data(d).name))
                    obj.data(d).outcomeEndTime = int_;
                end
            end
        end
        
        function int = get.softMerge(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).softMerge;
            end
        end
        
        function set.softMerge(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).softMerge = int_;
            end
        end
        
        function int = get.units(obj)
            int = [];
            if ~isempty(obj.data)
                int = obj.data(1).units;
            end
        end
        
        function set.units(obj, int_)
            for d = 1:numel(obj.data)
                obj.data(d).units = int_;
            end
        end
    end
    
    % variables dependent on data
    properties (Dependent)
        timeStamp
        sensorGlucose % array of sensor glucose (mmol/L)
        glucoseInterp % array of interpolated glucose (mmol/L)
        sensorCalibration % array recording when sensor calibration happens (logical)
        sensorScan % array recording when sensor is scanned (logical)
        bloodGlucose % array of blood glucose (mmol/L)
        basalInsulin % infused basal insulin
        basalInjection % injected long-acting insulin
        basalOverride % basal was overriden
        basalExternal % basal was external (true/false)
        bolusRecommend % algorithm bolus recommendation
        bolusInsulin % infused bolused insulin
        usualBolus % For fixed dose patients their usual bolus is saved here
        mealBolus % part of bolus considered to be meal carbohydrate
        corrBolus % part of bolus considered to be correction bolus
        bolusGlucagon % infused basal glucagon
        bolusExternal
        bolusUser
        bolusOverride % bolus was overriden
        carbs % carbohydrate intake (g)
        carbsFree % free carbohydrate intake (g)
        treats % treatment for hypoglycemia
        mealType % a meal parameter set by user, i.e. 1: breakfast, 2: lunch, ...
        carbFactors % programmed carb factors
        pumpBasals % programmed pump basals
        closedLoopActive % true if closed-loop is active
        time % array for time (min), starts at 0 and increments with stepTime
        startDate
        endDate
        duration
        startTime
        basalPramlintide
        bolusPramlintide
        steps
        heartrate
        exercise
        tdd
        tddBasal
        tddBolus
    end
    
    methods
        function t = get.timeStamp(obj)
            t = unique(vertcat(obj.data.timeStamp));
        end
        
        function t = get.sensorGlucose(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.glucoseInterp(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.sensorCalibration(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.sensorScan(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.bloodGlucose(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.basalInsulin(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.basalPramlintide(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.bolusPramlintide(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.basalInjection(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.basalOverride(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.basalExternal(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.bolusRecommend(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.bolusInsulin(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.usualBolus(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.mealBolus(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.corrBolus(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end

        function t = get.bolusGlucagon(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.bolusExternal(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.bolusOverride(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.bolusUser(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.carbs(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.treats(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.mealType(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.carbFactors(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.pumpBasals(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.closedLoopActive(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.steps(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.heartrate(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.time(obj)
            t = [];
            if ~isempty(obj.timeStamp)
                t = round((obj.timeStamp - obj.timeStamp(1))*24*60);
            end
        end
        
        function t = get.startDate(obj)
            t = min([obj.data.startDate]);
        end
        
        function t = get.endDate(obj)
            t = max([obj.data.endDate]);
        end
        
        function t = get.duration(obj)
            t = max([obj.data.duration]);
        end
        
        function t = get.startTime(obj)
            t = [];
            if ~isempty(obj.data)
                t = mode([obj.data.startTime]);
            end
        end
        
        function t = get.exercise(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.tdd(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.tddBasal(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
        
        function t = get.tddBolus(obj)
            st = dbstack;
            var = textscan(st(1).name, 'MAPPatient.get.%s', 1);
            t = obj.getField(var{1}{1});
        end
    end
    
    methods (Static)
        function cst = getUnchangedFields()
            cst = { ...
                'stepTime', ...
                'startDateLimit', ...
                'endDateLimit', ...
                'timestampString', ...
                'time', ...
                'startTime', ...
                'duration', ...
                'startDate', ...
                'endDate', ...
                'isEmpty', ...
                'iobBolus', ...
                'iobBasal', ...
                'glucoseInterp', ...
                'hypoEvents', ...
                'studyName', ...
                'ID', ...
                'name', ...
                'tddPerWeight', ...
                'bmi', ...
                'female', ...
                'durationDiabetes', ...
                'age', ...
                'ageToday', ...
                'arms', ...
                };
        end
        
        function cst = getStaticFields()
            cst = { ...
                'nested', ...
                'studyName', ...
                'ID', ...
                'TDD', ...
                'TDB', ...
                'weight', ...
                'height', ...
                'hba1c', ...
                'gender', ...
                'diagnosticYear', ...
                'birthDate', ...
                'admissionDate', ...
                'studyStartDate', ...
                'sensorType', ...
                'sensorDuration', ...
                'therapyType', ...
                'therapyParam', ...
                'startDate', ...
                'endDate', ...
                'duration', ...
                'startTime', ...
                'isEmpty', ...
                'stepTime', ...
                'name', ...
                'tddPerWeight', ...
                'bmi', ...
                'female', ...
                'durationDiabetes', ...
                'age', ...
                'ageToday', ...
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
                'onlyCountBiggerMeal', ...
                'computeIOB', ...
                'outcomeStartTime', ...
                'outcomeEndTime', ...
                'softMerge', ...
                };
        end
        
        function signals = getTimeFields()
            signals = { ...
                'timeStamp', ...
                'sensorGlucose', ...
                'glucoseInterp', ...
                'sensorCalibration', ...
                'sensorScan', ...
                'bloodGlucose', ...
                'basalInsulin', ...
                'basalInjection', ...
                'basalOverride', ...
                'basalExternal', ...
                'bolusRecommend', ...
                'bolusInsulin', ...
                'usualBolus', ...
                'mealBolus', ...
                'corrBolus', ...
                'bolusGlucagon', ...
                'bolusExternal', ...
                'bolusOverride', ...
                'bolusUser', ...
                'carbs', ...
                'exercise',...
                'treats', ...
                'mealType', ...
                'carbFactors', ...
                'pumpBasals', ...
                'closedLoopActive', ...
                'time', ...
                'basalPramlintide', ...
                'bolusPramlintide', ...
                'steps', ...
                'heartrate', ...
                'tdd',...
                'tddBasal',...
                'tddBolus',...
                };
        end
    end
    
    methods (Access = public)
        
        %% Constructor
        function obj = MAPPatient(varargin)
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'patient'
                        obj.ID = varargin{nVar+1}.ID;
                        obj.fromStruct(varargin{nVar+1});
                    case 'name'
                        allNumbers = regexp(varargin{nVar+1}, '\d+.?', 'match');
                        if ~isempty(allNumbers) && ~isnan(str2double(allNumbers(end)))
                            obj.ID = str2double(allNumbers(end));
                        else
                            obj.ID = 1;
                        end
                    case 'id'
                        obj.ID = varargin{nVar+1};
                    otherwise
                        error('[MAPPatient] Unkown option %s, use doc MAPPatient for more information.', varargin{nVar});
                end
            end
        end
        
        function param = getProperties(obj)
            param = struct();
            for fn = {'ID', 'TDD', 'TDB', 'weight', 'height', 'hba1c', 'gender', 'diagnosticYear', 'birthDate', 'admissionDate'}
                param.(fn{1}) = obj.(fn{1});
            end
        end
        
        function add(obj, mapData)
            if ~isa(mapData, 'MAPData')
                error('[MAPPatient][addMAPData] Only supports adding a MAPData object.');
            end
            
            idxData = find(strcmpi({obj.data.name}, mapData.name), 1);
            
            if isempty(idxData)
                obj.data(end+1) = mapData.copy;
                
                obj.data(end).patientName = obj.name;
                obj.data(end).patientHandle = obj;
                
                obj.data(end).studyName = obj.studyName;
            else
                obj.data(idxData) = mapData.copy;
                
                obj.data(idxData).patientName = obj.name;
                obj.data(idxData).patientHandle = obj;
                
                obj.data(idxData).studyName = obj.studyName;
            end
        end
        
        function removeData(obj, name)
            if isnumeric(name) || islogical(name)
                obj.data(name) = [];
            else
                obj.data(strcmpi({obj.data.name}, name)) = [];
            end
        end
        
        function res = hasData(obj, name)
            if numel(obj) > 1
                res = arrayfun(@(c)(c.hasData(name)), obj, 'UniformOutput', false);
                res = [res{:}];
                return;
            end

            res = ~isempty(obj.getData(name));
        end
        
        function mapData = getData(obj, name)
            if isempty(obj)
                mapData = MAPData.empty;
                return;
            end
            
            if nargin < 2
                disp(strcat(num2str((1:length(obj.data))', '%d -> '), {obj.data.name}'));
                return;
            end
            
            if numel(obj) > 1
                mapData = arrayfun(@(c)(c.getData(name)), obj, 'UniformOutput', false);
                mapData(cellfun(@isempty, mapData)) = [];
                mapData = [mapData{:}];
                return;
            end
            
            dataNames = {obj.data.name};
            
            if isempty(obj) || isempty(dataNames)
                mapData = MAPData.empty;
                return;
            end
            
            if ischar(name)
                if contains(name, '*') % handles wildcard
                    mapData = obj.data(cellfun(@(c1)all(cellfun(@(c2)contains(c1, c2, 'IgnoreCase', true), strsplit(name, '*'))), dataNames));
                else
                    mapData = obj.data(strcmpi(dataNames, name));
                end
            else % name is cell array or array
                for idx = numel(name):-1:1
                    if isnumeric(name(idx))
                        mapData(idx) = obj.data(name(idx));
                    else
                        if any(strcmpi(dataNames, name{idx}))
                            if contains(name{idx}, '*') % handles wildcard
                                mapData(idx) = obj.data(cellfun(@(c1)all(cellfun(@(c2)contains(c1, c2, 'IgnoreCase', true), strsplit(name{idx}, '*'))), dataNames));
                            else
                                mapData(idx) = obj.data(strcmpi(dataNames, name{idx}));
                            end
                        else
                            mapData = MAPData.empty;
                            return;
                        end
                    end
                end
            end
        end
        
        function data_ = makeData(obj, name_, varargin)
            startDate_ = NaT;
            endDate_ = NaT;
            duration_ = [];
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'startdate' % this is inclusive
                        startDate_ = varargin{nVar+1};
                    case 'enddate' % this is exclusive
                        endDate_ = varargin{nVar+1};
                    case 'duration'
                        duration_ = varargin{nVar+1};
                    otherwise
                        error('[MAPPatient][makeData] Unkown option %s, use doc MAPData for more information.', varargin{nVar});
                end
            end
            
            if all(isnat(startDate_)) && all(~isempty(duration_)) && all(~isnat(endDate_))
                startDate_ = endDate_ - duration_ + minutes(obj.stepTime);
            elseif all(isnat(endDate_)) && all(~isempty(duration_)) && all(~isnat(startDate_))
                endDate_ = startDate_ + duration_ - minutes(obj.stepTime);
            elseif all(isnat(startDate_)) && all(isnat(endDate_))
                error('[MAPPatient][makeData] You need to specify at least a start date or an end date with duration.');
            end
            
            s = struct();
            minStartDate_ = min(startDate_);
            
            s.time = [];
            for n = 1:numel(startDate_)
                s.time = [s.time; (minutes(startDate_(n)-minStartDate_):obj.stepTime:(minutes(endDate_(n)-minStartDate_)))'];
            end
            if isempty(s.time)
                data_ = MAPData.empty;
                return;
            end
            
            s.time = sort(s.time);
            
            idxR = any(abs(obj.time-(s.time + minutes(minStartDate_-obj.startDate))') < obj.stepTime/2, 2);
            idxL = any(abs((s.time + minutes(minStartDate_-obj.startDate))-obj.time') < obj.stepTime/2, 2);
            if all(idxR == 0) || all(idxL == 0)
                data_ = MAPData.empty;
                return;
            end
            
            for fn = setdiff(obj.getTimeFields, obj.getUnchangedFields)
                if strcmp(fn{1}, 'time')
                    continue;
                end
                
                s.(fn{1}) = nan(size(s.time));
                if ~isempty(obj.(fn{1}))
                    s.(fn{1})(idxL) = obj.(fn{1})(idxR);
                end
            end
            
            data_ = MAPData('name', name_, ...
                'stepTime', obj.stepTime, ...
                'startDate', min(startDate_), ...
                'endDate', max(endDate_), ...
                'data', s);
            
            data_.patientName = obj.name;
            data_.patientHandle = obj;
            
            data_.studyName = obj.studyName;
            
            if nargout == 0
                obj.add(data_);
            end
        end
        
        %% Functions to load data
        function fromStruct(obj, s)
            for fn = setdiff(obj.getStaticFields, obj.getUnchangedFields)
                if isfield(s, fn{1})
                    obj.(fn{1}) = s.(fn{1});
                end
            end
            
            if isfield(s, 'data')
                nbrData = numel(s.data);
                if nbrData > 0
                    for d = nbrData:-1:1
                        obj.data(d) = MAPData('data', s.data(d));
                        
                        obj.data(d).patientName = obj.name;
                        obj.data(d).patientHandle = obj;
                        
                        obj.data(d).studyName = obj.studyName;
                    end
                end
            end
        end
        
        %% Functions to save data
        function s = toStruct(obj, varargin)
            average_ = false;
            arms_ = {obj.data.name};
            days_ = [];
            startTime_ = [];
            outcomes_ = false;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'average'
                        average_ = varargin{nVar+1};
                    case 'days'
                        days_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'starttime'
                        startTime_ = varargin{nVar+1};
                    case 'outcomes'
                        outcomes_ = varargin{nVar+1};
                end
            end
            
            arms_ = intersect(arms_, {obj.data.name});
            
            s = struct();
            
            for fn = obj.getStaticFields()
                s.(fn{1}) = obj.(fn{1});
            end
            
            for d = 1:numel(arms_)
                if average_
                    s.data(d) = obj.getData(arms_(d)).toStruct('format', 'average', ...
                        'outcomes', outcomes_, ...
                        'days', days_, ...
                        'starttime', startTime_);
                else
                    s.data(d) = obj.getData(arms_(d)).toStruct('format', 'array', ...
                        'outcomes', outcomes_, ...
                        'days', days_, ...
                        'starttime', startTime_);
                end
            end
        end
        
        function sheet = toSheet(obj, varargin)
            sheet = cell(0);
            
            for d = 1:numel(obj.data)
                if ~obj.data(d).outcome
                    continue;
                end
                
                sheet = MAPUtils.concatSheets(sheet, obj.data(d).toSheet(varargin{:}));
                if d < numel(obj.data)
                    sheet = MAPUtils.concatSheets(sheet, cell(size(sheet, 1), 1));
                end
            end
        end
        
        function toPNG(obj, path, varargin)
            inScreen_ = false;
            inFolder_ = true;
            daily_ = false;
            nested_ = obj.nested;
            arms_ = obj.arms;
            summary_ = true;
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
                    case 'arms'
                        arms_ = varargin{nVar+1};
                    case 'daily'
                        daily_ = varargin{nVar+1};
                    case 'average'
                        average_ = varargin{nVar+1};
                    case 'summary'
                        summary_ = varargin{nVar+1};
                end
            end
            
            if nested_ > 0
                dPath = [path, obj.name, filesep];
            else
                dPath = [path, obj.name, '_'];
            end
            for d = 1:numel(obj.data)
                if ~any(strcmpi(arms_, obj.data(d).name))
                    continue;
                end
                
                obj.data(d).toPNG(dPath, ...
                    varargin{:}, ...
                    'inscreen', inScreen_, ...
                    'infolder', inFolder_, ...
                    'daily', daily_, ...
                    'average', average_, ...
                    'nested', nested_-1);
            end
            
            if summary_ && numel(arms_) > 1
                fig = obj.toFigure( ...
                    'inscreen', inScreen_, ...
                    'daily', daily_, ...
                    'average', average_, ...
                    'arms', arms_);
                
                if ~isempty(fig)
                    dPath = [path, obj.name, 'Summary'];
                    obj.ensureFolder(dPath);
                    print(fig, [dPath, '.png'], '-dpng');
                    
                    if ~inScreen_
                        close(fig);
                    end
                end
            end
        end
        
        function fig = toFigure(obj, varargin)
            inScreen_ = true;
            average_ = false;
            arms_ = obj.arms;
            title_ = '';
            days_ = [];
            stats_ = true;
            summaryFct_ = @obj.plotSummary;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'inscreen'
                        inScreen_ = varargin{nVar+1};
                    case 'average'
                        average_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = intersect(varargin{nVar+1}, {obj.data.name});
                    case 'title'
                        title_ = varargin{nVar+1};
                    case 'days'
                        days_ = varargin{nVar+1};
                    case 'stats'
                        stats_ = varargin{nVar+1};
                    case {'summaryfun', 'summaryfct', 'summaryfunction'}
                        summaryFct_ = varargin{nVar+1};
                end
            end
            
            data_ = cell(numel(arms_), 1);
            for n = 1:numel(arms_)
                if ~any(strcmp(arms_{n}, {obj.data.name}))
                    continue;
                end
                
                if average_
                    dStruct = obj.getData(arms_{n}).toStruct( ...
                        'days', days_, ...
                        'outcomes', true, ...
                        'format', 'average');
                else
                    dStruct = obj.getData(arms_{n}).toStruct( ...
                        'days', days_, ...
                        'outcomes', true, ...
                        'format', 'array');
                end
                
                data_{n} = [data_{n}, dStruct];
                if any(mode([data_{n}.duration]) ~= [data_{n}.duration])
                    warning('[MAPPatient][toFigure] Data duration for %s is different between interventions ... can''t show summary.', obj.name);
                    fig = gobjects(0);
                    return;
                end
            end
            
            if sum(~cellfun(@isempty, data_)) == 1
                warning('[MAPPatient][toFigure] Only one arm exist for %s. Nothing to show.', obj.name);
                fig = gobjects(0);
                return;
            end
            
%             if numel(unique(round(minutes([cell2mat(data_).duration])))) > 1
%                 warning('[MAPPatient][toFigure] Data duration is different ... can''t show summary.');
%                 fig = gobjects(0);
%                 return;
%             end

            duration_ = minutes(data_{1}.duration);
            
            if inScreen_
                fig = figure(1110);
            else
                fig = figure('Visible', 'Off');
            end
            clf(fig);
            set(fig, 'name', sprintf('MAPStudy::toFigure::%s', obj.name), ...
                'numbertitle', 'off', ...
                'units', 'normalized', ...
                'outerposition', [0, 0, 1, 1]);
            
            if stats_
                axisPositions = [0.04, 0.15, 0.75, 0.75];
            else
                axisPositions = [0.04, 0.15, 0.81, 0.75];
            end
            
            ax = subplot('Position', axisPositions);
            
            if isempty(title_)
                title_ = sprintf(['%s Summary of ', repmat('%s, ', 1, numel(arms_)), '%5.2f hours'], obj.name, arms_{:}, hours(minutes(duration_)));
            end
            summaryFct_(ax, data_, ...
                'showbolus', false, ...
                'showhypos', true, ...
                'duration', duration_, ...
                'title', title_, ...
                varargin{:});
            ax.Position = axisPositions;
            if stats_
                MAPUtils.plotStatsSummary(fig, data_)
            end
        end
        
        function outcomeTab = outcomes(obj, varargin)
            interval_ = [];
            arms_ = obj.arms;
            
            % Get options
            for nVar = 1:2:length(varargin)
                switch lower(varargin{nVar})
                    case 'interval'
                        interval_ = varargin{nVar+1};
                    case 'arms'
                        arms_ = varargin{nVar+1};
                end
            end
            
            outcomeTab = table;
            for d = 1:numel(obj.data)
                if ~any(strcmpi(obj.data(d).name, arms_))
                    continue;
                end
                
                if isempty(outcomeTab)
                    outcomeTab = obj.data(d).outcomes(varargin{:}, 'interval', interval_);
                else
                    newOutcomeTab = obj.data(d).outcomes(varargin{:}, 'interval', interval_);
%                     a=newOutcomeTab.Properties.RowNames;
%                     b=outcomeTab.Properties.RowNames;
                    newOutcomeTab.Properties.RowNames = outcomeTab.Properties.RowNames;
                    temp = join(outcomeTab, newOutcomeTab, 'Keys', 'RowNames');
                    outcomeTab = temp;
                end
            end
        end
        
        function outcomeTab = glucoseOutcomes(obj, interval_)
            if nargin < 2
                interval_ = [];
            end
            
            outcomeTab = table;
            for d = 1:numel(obj.data)
                if ~obj.data(d).outcome
                    continue;
                end
                
                if isempty(outcomeTab)
                    outcomeTab = obj.data(d).glucoseOutcomes(interval_);
                else
                    outcomeTab = join(outcomeTab, obj.data(d).glucoseOutcomes(interval_), 'Keys', 'RowNames');
                end
            end
        end
        
        function outcomeTab = insulinOutcomes(obj, interval_)
            if nargin < 2
                interval_ = [];
            end
            
            outcomeTab = table;
            for d = 1:numel(obj.data)
                if ~obj.data(d).outcome
                    continue;
                end
                
                if isempty(outcomeTab)
                    outcomeTab = obj.data(d).insulinOutcomes(interval_);
                else
                    outcomeTab = join(outcomeTab, obj.data(d).insulinOutcomes(interval_), 'Keys', 'RowNames');
                end
            end
        end
        
        function outcomeTab = mealsOutcomes(obj, interval_)
            if nargin < 2
                interval_ = [];
            end
            
            outcomeTab = table;
            for d = 1:numel(obj.data)
                if ~obj.data(d).outcome
                    continue;
                end
                
                if isempty(outcomeTab)
                    outcomeTab = obj.data(d).mealsOutcomes(interval_);
                else
                    outcomeTab = join(outcomeTab, obj.data(d).mealsOutcomes(interval_), 'Keys', 'RowNames');
                end
            end
        end
        
        function outcomeTab = otherOutcomes(obj, interval_)
            outcomeTab = table;
            for d = 1:numel(obj.data)
                if ~obj.data(d).outcome
                    continue;
                end
                
                if isempty(outcomeTab)
                    outcomeTab = obj.data(d).otherOutcomes(interval_);
                else
                    outcomeTab = join(outcomeTab, obj.data(d).otherOutcomes(interval_), 'Keys', 'RowNames');
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
            
            detailsTab = table;
            for d = 1:numel(obj.data)
                if any(strcmp(arms_, obj.data(d).name))
                    if isempty(detailsTab)
                        detailsTab = obj.data(d).details();
                    else
                        detailsTab = join(detailsTab, obj.data(d).details(), 'Keys', 'RowNames');
                    end
                end
            end
        end
        
        function out_ = getField(obj, fn)
            out_ = vertcat(obj.data.(fn));
            if ~isempty(out_)
                timeStamp_ = vertcat(obj.data.timeStamp);
                if length(timeStamp_) ~= length(out_)
                    out_ = [];
                    for n = 1:numel(obj.data)
                        if ~isempty(obj.data(n).(fn))
                            out_ = [out_; obj.data(n).(fn)];
                        else
                            out_ = [out_; nan(size(obj.data(n).timeStamp))];
                        end
                    end
                end
                [sortedTimeStamp, iSort] = sort(timeStamp_);
                [~, iUnique] = unique(sortedTimeStamp);
                out_ = out_(iSort(any(iUnique == 1:length(timeStamp_))));
            end
        end
    end
end
