%% The following script can be used for unit conversion AND contains important physical constants
function output = physics_convert(convertType,value)
%% Constants
g = 9.81; %m/s^2
c = 2.9979e8;  %m/s
G = 6.67408e-11;  %Nm^2/kg^2
mass_earth = 5.97e24; %kg
radius_earth = 6371e3; %m
h_js = 6.626e-34; %Js --> Planck's Constant
h_evs = 4.1357e-15; %eVs
hc_jm = 1.9864e-25; %Jm
hc_evnm = 1239.84197; %eVnm
h_bar_js = h_js/(2*pi); %Js
h_bar_evs = h_evs/(2*pi); %eVs
e0 = 8.854e-12; %C^2/Nm^2 --> Vacuum permitivity
u0 = 4*pi*10^-7; %Tm/A --> Magnetic Constant/Vacuum Permeability
e = 1.602e-19; %C --> Elementary Charge
mass_electron_kg = 9.11e-31; %kg
mass_electron_mev = .511; %MeVs
mass_electron_amu = .00055; %amu
mass_proton_kg= 1.673e-27; %kg
mass_proton_mev = 938.272; %MeV
mass_proton_amu = 1.00728; %amu
mass_neutron_kg = 1.6749e-27; %kg
mass_neutron_mev = 939.565; %MeV
mass_neutron_amu = 1.008665; %amu
amu_kg = 1.6605e-27; %kg
amu_mev = 939.5654; %MeV
amu = 1; %amu
avogadro = 6.022e23; %1/mol
k = 1.3806e-23; %J/K --> Boltzmann's Constant
R = 8.3145; %J/molK --> Gas Constant
sigma = 5.670367e-8; %W/m^2K^4 --> Stefan-boltzmann Constant
beta_mmK = 2.8978; %mm/K --> Wien Displacement Constant
beta_ghzk = 58.789; %GHz/K
atm = 101.3e3; %Pa --> Standard Atmospheric Pressure
H0_s = 2.25e-18; %s^-1 --> Hubble Constant
H0 = 69.3; %km/s/Mpc
planck_length = 1.616229e-35; %m
planck_mass = 2.176470e-8; %kg
planck_time = planck_length/c; %s
u_B = 9.274e-24; %J/T --> Bohr Magneton
bohr_radius = 5.291772e-11; %m --> Bohr Radius
electron_radius = 2.8179e-15; %m
fine_structure = 7.297352e-3; %unitless
kb_ev = 8.617e-5; %boltzmann constant in eV/K
radius_sun = 695508*1000; % m
mass_sun = 1.989e30; %kg
save all;


%% Conversions  
if strcmp('m2mm', convertType)
    output = value*1000;
elseif strcmp(convertType, 'j2ev')
    output = value*6.242e18;
elseif strcmp(convertType, 'ev2j')
    output = value*1.6022e-19;
elseif strcmp(convertType, 'ev2kg')
    output = value*e/(c^2);
elseif strcmp(convertType, 'kg2ev')
    output = value*c^2/e;
    
% Astronomical value to kilometers
elseif (strcmp(convertType,'AU2km'))
    output = value*1.495978707e8; 
    
% Astronomical value to meters
elseif (strcmp(convertType,'AU2m'))
    output = value*1.49597870700e11;

% Kilometers to Astronomical value
elseif (strcmp(convertType,'km2AU'))
    output = value/1.495978707e8;

% Meters to Astronomical value
elseif (strcmp(convertType,'m2AU'))
    output = value/1.495978707e11;

% Arcseconds to degrees
elseif (strcmp(convertType,'arcsec2deg'))
    output = value*3600; % check here
    
% Arcseconds to radians
elseif (strcmp(convertType,'arcsec2rad'))
    output = value*3600*pi/180;
    
% Degrees to arcseconds
elseif (strcmp(convertType,'deg2arcsec'))
    output = value/60;
    
% Radians to arcseconds
elseif (strcmp(convertType,'rad2arcsec'))
    output = value*180/pi/60;
    
% Lightyears to Astronomical value
elseif (strcmp(convertType,'ly2AU'))
    output = value*63241.077;
    
% Lightyears to Parsecs
elseif (strcmp(convertType,'ly2pc'))
    output = value/3.26156;
    
% Astronomical value to Light Years
elseif (strcmp(convertType,'AU2ly'))
    output = value/63241.077;
    
% Astronomical value to Parsecs
elseif (strcmp(convertType,'AU2pc'))
    output = value/2.06265e5;
    
% Parsecs to Light Years
elseif (strcmp(convertType,'pc2ly'))
    output = value*3.26156;
    
% Parsecs to Astronomical value
elseif (strcmp(convertType,'pc2AU'))
    output = value*2.06265e5;
elseif (strcmp(convertType, 'pc2km'))
    output = value*3.086e+13;
elseif (strcmp(convertType,'day2sec'))
    output = value*86400;
elseif (strcmp(convertType,'year2sec'))
    output = value*3.154e+7;
elseif (strcmp(convertType,'sec2year'))
    output = value*3.17098e-8;
elseif (strcmp(convertType,'amu2kg'))
    output = value*1.66054e-27;
    
else
    disp('Sorry, I dont know that one yet')
end
    

end
%% End of Script
