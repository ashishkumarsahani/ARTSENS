function [ inflection_X, L_Fitted, error, trial_Count, time_per_Iteration] = Fit_L(well_Half)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Tot_Len = length(well_Half);
L = 1*[1:Tot_Len];
L_Fitted = L;
error = 100;
inflection_X = Tot_Len;
trial_Count =-1;
time_per_Iteration = -1;
if(Tot_Len <=1)
    return
end

well_Half = well_Half-min(well_Half);
well_Half = well_Half/(max(well_Half));

tic;
[inflection_X, error, L_Fitted, trial_Count] = fitL(well_Half);
time_per_Iteration = toc/trial_Count;
end
