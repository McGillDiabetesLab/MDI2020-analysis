function lmeTab = lme(obj, varargin)

arms_ = obj.arms;
qqplot_ = ''; % path where qqplot is saved
interval_ = [];

% Get options
for nVar = 1:2:length(varargin)
    switch lower(varargin{nVar})
        case 'qqplot'
            qqplot_ = varargin{nVar+1};
        case 'arms'
            arms_ = varargin{nVar+1};
        case 'interval'
            interval_ = varargin{nVar+1};
    end
end

if ~isempty(qqplot_)
    qqplot_ = obj.ensureFolder(qqplot_);
end

data = obj.getData(arms_);

colName = {'ID', 'Intervention', 'Sequence', 'Period'};
colType = {'categorical', 'categorical', 'categorical', 'categorical'};

tbl = table('Size', [numel(data), numel(colName)], ...
    'VariableNames', colName, ...
    'VariableTypes', colType);
tbl.ID = categorical(arrayfun(@(d)(d.patientName), data, 'UniformOutput', false))';
tbl.Intervention = categorical(arrayfun(@(d)(d.name), data, 'UniformOutput', false))';

    function val = armsSequence(pMAP)
        [~, indices] = sort([pMAP.data.startDate]);
        orderedArms = pMAP.arms(indices);
        orderedArms([pMAP.data.outcome] == 0) = [];
        orderedArms(cellfun(@(c)any(strcmp(c, arms_)), orderedArms) == 0) = [];
        val = join(orderedArms, '_');
        val = val{1};
    end

tbl.Sequence = categorical(arrayfun(@armsSequence, [data.patientHandle], 'UniformOutput', false))';
tbl.Period = categorical(arrayfun(@(c)(find(strcmp(split(string(tbl.Sequence(c)), '_'), string(tbl.Intervention(c))))), 1:height(tbl)))';

warning('off', 'MATLAB:table:ModifiedVarnamesRows2Vars');
out = arrayfun(@(d)(d.outcomes(varargin{:})), data, 'UniformOutput', false);
outT = cellfun(@(o)rows2vars(o), out, 'UniformOutput', false);
tbl = [tbl, removevars(vertcat(outT{:}), 'OriginalVariableNames')];
tbl = sortrows(tbl, 'ID');

rowName = { ...
    'P-Value Sequence Effects', ... %1
    'P-Value Period Effects', ... %2
    'Shapiro Normality Test', ... %3
    '__Normal Distribution__', ...
    'P-Value Intervention Effects', ... %5
    'Paired t-test (non-adjusted)', ... %6
    'Fixed Effects', ... %7
    '2.5% Quantile', ... %8
    '97.5% Quantile', ... %9
    '__Non-Normal Distribution__', ...
    'P-Value Intervention Effects (SQRT)', ... %11
    'P-Value Wilcoxon Test', ... %12
    'Fixed Effects (Hodges-Lehmann)', ... %13
    '2.5% Quantile (Hodges-Lehmann)', ... %14
    '97.5% Quantile (Hodges-Lehmann)', ... %15
    };

warning('off', 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
tabRes = table('Size', [numel(rowName), numel(colName)], ...
    'RowNames', rowName, ...
    'VariableNames', colName, ...
    'VariableTypes', colType);

for n = (numel(colName) + 1):length(tbl.Properties.VariableNames)
    pValues = nan(size(tabRes, 1), 1);
    
    % Paired test p values
    out_ = arrayfun(@(c)(tbl.(tbl.Properties.VariableNames{n})(tbl.Intervention == c)), arms_, 'UniformOutput', false);
    out_ = [out_{:}];
    [~, pValues(6)] = ttest(out_(:, 1), out_(:, 2));
    
    % Shapiro normality test
    pValues(3) = MAPUtils.swtest(out_(:, 1)-out_(:, 2));
    
    % Wilcoxon test p values
    [pValues(12), stats] = MAPUtils.utest(out_(:, 1), out_(:, 2));
    pValues(13) = stats.median;
    pValues(14) = stats.ci(1);
    pValues(15) = stats.ci(2);
    
    % LME
    try
        lm = [];
        eval(sprintf('lm = fitlme(tbl, ''%s ~ Intervention + Sequence + Period + (1|ID)'', ''FitMethod'', ''REML'');', ...
            tbl.Properties.VariableNames{n}));
        an = anova(lm, 'dfmethod', 'satterthwaite');
        ci = coefCI(lm, 'dfmethod', 'satterthwaite');
        
        pValues(1) = an.pValue(strcmp(an.Term, 'Sequence'));
        pValues(2) = an.pValue(strcmp(an.Term, 'Period'));
        pValues(5) = an.pValue(strcmp(an.Term, 'Intervention'));
        
        if find(strcmp(cellfun(@(c)(['Intervention_', c]), arms_, 'UniformOutput', false), lm.CoefficientNames(contains(lm.CoefficientNames, 'Intervention')))) == 1
            pValues(7) = lm.Coefficients.Estimate(contains(lm.CoefficientNames, 'Intervention'));
            pValues(8) = ci(contains(lm.CoefficientNames, 'Intervention'), 1);
            pValues(9) = ci(contains(lm.CoefficientNames, 'Intervention'), 2);
        else
            pValues(7) = -lm.Coefficients.Estimate(contains(lm.CoefficientNames, 'Intervention'));
            pValues(8) = -ci(contains(lm.CoefficientNames, 'Intervention'), 2);
            pValues(9) = -ci(contains(lm.CoefficientNames, 'Intervention'), 1);
        end
        
        % qqplot
        if ~isempty(qqplot_)
            h = figure('Visible', 'off');
            clf;
            qqplot(lm.residuals);
            ax = gca;
            ax.YAxis.FontSize = 14;
            ax.XAxis.FontSize = 14;
            ax.FontWeight = 'bold';
            ax.LineWidth = 1.0;
            if isempty(interval_)
                print(h, [qqplot_, '/qqplot_', tbl.Properties.VariableNames{n}], '-dpng');
            else
                print(h, [qqplot_, '/qqplot_', tbl.Properties.VariableNames{n}, '_', num2str(interval_(1)/60), '_', num2str(interval_(2)/60)], '-dpng');
            end
        end
    catch
        % fitlme may fail if there is not enough data
    end
    
    % LME Skewed
    try
        tblSkewed = tbl;
        tblSkewed.(tbl.Properties.VariableNames{n}) = sqrt(tblSkewed.(tbl.Properties.VariableNames{n}));
        lmSkewed = [];
        eval(sprintf('lmSkewed = fitlme(tblSkewed, ''%s ~ Intervention + Sequence + Period + (1|ID)'', ''FitMethod'', ''REML'');', ...
            tbl.Properties.VariableNames{n}));
        anSkewed = anova(lmSkewed, 'dfmethod', 'satterthwaite');
        
        pValues(11) = anSkewed.pValue(strcmp(anSkewed.Term, 'Intervention'));
        
        % qqplot
        if ~isempty(qqplot_)
            h = figure('Visible', 'off');
            clf;
            qqplot(lmSkewed.residuals);
            ax = gca;
            ax.YAxis.FontSize = 14;
            ax.XAxis.FontSize = 14;
            ax.FontWeight = 'bold';
            ax.LineWidth = 1.0;
            if isempty(interval_)
                print(h, [qqplot_, '/qqplot_', tbl.Properties.VariableNames{n}, '_sqrt'], '-dpng');
            else
                print(h, [qqplot_, '/qqplot_', tbl.Properties.VariableNames{n}, '_', num2str(interval_(1)/60), '_', num2str(interval_(2)/60), '_sqrt'], '-dpng');
            end
        end
    catch
        % fitlme may fail if there is not enough data
    end
    
    % Fill the table
    tabRes.(tbl.Properties.VariableNames{n}) = pValues;
end
lmeTab = [tbl; tabRes];

end
