function sebum_flow_rate = estimate_sebum_flow_speed

%% Summary

% Estimates sebum flow speed upward in a sebaceous follicle

% Note1: There are many sebum excretion measurements in the literature.
% However, these pertain to the amount of sebum excreted onto the skin
% surface shortly after removing oils from the face, NOT the amount of
% sebum produced by sebaceous glands at the bottom of the pore. Converting
% these estimates to sebum flow per pore results in a much higher estimate,
% but is probably not accurate because these excretion rates are not
% typically maintained over time, and appear to be more of a result the
% skin surface disturbance.

% Note2: This assumes a *constant* and *uniform* sebum production rate by
% sebaceous glands.


%% Estimate

length_of_pore = 1000; % um
time_to_refill = 20; % days; SOURCE: Plewig 1974; Acne Book 2019

sebum_flow_rate = length_of_pore/time_to_refill; % um/day


end