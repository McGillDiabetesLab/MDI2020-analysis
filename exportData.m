clc,
clear,
%% Load Study
sMAP = MAPStudy();
sMAP.fromMAT('./Data/MDI2020.mat');

%%
data = sMAP.getData('*');
for k = 1:length(data)
    data(k).name = data(k).name(5:end);
end
sMAP.toExcel('summary.xlsx', 'individual', true, 'summary', true, 'compare', true)