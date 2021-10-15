function diffusion_coefficient = estimate_diffusion_coefficient_of_cell_in_pore()

%% Summary

% Estimates diffusion coefficient of a bacterial cell in sebum

% Note that actual diffusion is probably slower due to obstacles (dead
% human cells, other bacterial cells if the pore is already colonized,
% etc.)


%% Estimate

% Boltzmann constant
kB_SIunits = 1.38*10^(-23); % J/K = kg*m^2/sec^2/Kelvin
seconds_per_day = 60*60*24; % sec/d; 24 hrs/day * 60 min/hr * 60 sec/min
microns_per_meter = 10^6; 
kB = kB_SIunits * seconds_per_day^2 * microns_per_meter^2; % kg*um^2/day^2/Kelvin

% Temperature
body_temperature_Celsius = 37; % Celsius
celsius_to_kelvin = 273.15; % additive
body_temperature = body_temperature_Celsius + celsius_to_kelvin; % Kelvin

% Viscosity of sebum
% Source: The Physical Properites of Sebum, Butcher and Coonin, JID 1949; https://core.ac.uk/download/pdf/82361958.pdf
sebum_viscosity_millipoises = 600; % millipoises 
millipoises_per_pascalsec = 10^4;
sebum_viscosity_pascalsec = sebum_viscosity_millipoises / millipoises_per_pascalsec; % Pa*sec = kg/m/sec
sebum_viscosity = sebum_viscosity_pascalsec * seconds_per_day / microns_per_meter; % kg/um/day^2

% Radius of a cell
cell_radius = 0.5; % um; assumes diameter of a bacterium is 1um and that cell is spherical

% Stokes-Einstein equation
diffusion_coefficient = kB*body_temperature/(6*pi*sebum_viscosity*cell_radius); % um^2/day


end



