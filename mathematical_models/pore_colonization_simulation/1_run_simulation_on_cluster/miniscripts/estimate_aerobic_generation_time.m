function t_generation = estimate_aerobic_generation_time()

%% Summary

% Estimates generation time at top of pore

% Assumes aerobic conditions
% Assumes growth rate is same as in the lab (probably not true)

%% Estimate

gamma_aerobic = 0.06*24; % /day; SOURCE: Cove 1983 (max aerobic growth rate = 0.06/hr)
t_generation = log(2)/gamma_aerobic; % days

% Note: Anaerobic growth rate reported in Cove aligns well with our 
% observation that colonies (~10^7 cells?) form within ~5 days.


end