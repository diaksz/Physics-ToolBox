%% Astrophysics Toolbox %%
% Authors: Andrew Winter and Zachary Diaks
% Description: Allows for creating an object that, upon instantiation, adds
% many useful constants to the user's workspace. Also gives access to many
% useful astrophysical conversions.
%
% HOW TO USE:
% In desired script, just assign the desired name equal to this class.
% Ex: toolbox = astroToolbox; OR toolbox = astroToolbox();
% If you want all of the constants added to your workspace, add any arbitrary argument
% Ex: toolbox = astroToolbox(hello); or toolbox = astroToolbox(1);

%% Class Definition %%
classdef astroToolbox
    %ASTROPHYSICSTOOLBOX Useful astrophysical constants and conversions
    
    % Debugging: Default value set to false so we know when class isn't working
    properties
        classWorking = false;
    end
    
    % Non-static methods (constructor, conversions)
    methods
        % Constructor
        % Adds many useful constants to user's workspace
        function this = astroToolbox(varargin)
            if (isempty(varargin))
                this.classWorking = true;
            elseif (~isempty(varargin))
                this.classWorking = true;
                this.constants();
            end
        end
        
        % Conversion Function
        %%% PLEASE NOTE: When adding stuff to this, don't use constants. Be
        %%% explicit and type in the whole number (i.e. don't type beta =
        %%% value/c, type beta = value/3e8. This is because those values
        %%% are not necessarily declared prior to when this function is used
        %
        %%% How to use: objectname.convertUnits('str',value)
        function output = convertUnits(~,convertType,value)
            if strcmp('m2mm', convertType)
                output = value*1000;
            elseif strcmp(convertType, 'j2ev')
                output = value*6.242e18;
            elseif strcmp(convertType, 'ev2j')
                output = value/6.242e18;
            elseif strcmp(convertType, 'ev2kg')
                output = value*(1.602e-19)/((3e8)^2);
            elseif strcmp(convertType, 'kg2ev')
                output = value*(3e8)^2/(1.602e-19);

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
                output = value/3600;

            % Arcseconds to radians
            elseif (strcmp(convertType,'arcsec2rad'))
                output = value*pi/(180*3600);

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
            
            % Unknown conversion
            else
                disp('Sorry, I dont know that one yet')
            end
        end
    end
    
    % Static methods(adds the constants to the base workspace)
    methods (Static)
        function constants()
            %Adds all constants to workspace
            assignin('base','g',9.81); % [m/s^2]
            assignin('base','c',2.9979e8); % [m/s]
            assignin('base','G',6.67408e-11); % [Nm^2/kg^2]
            assignin('base','mEarth',5.97e24); % [kg]
            assignin('base','rEarth',6371e3); % [m]
            assignin('base','h_js',6.626e-34); % [Js] --> Planck's Constant
            assignin('base','h_evs',4.1357e-15); % [eVs]
            assignin('base','hc_jm',1.9864e-25); % [Jm]
            assignin('base','hc_evnm',1239.84197); % [eVnm] --> Placnk's constant
            assignin('base','hBar_js',6.626e-34/(2*pi)); %[Js] --> Planck's Constant / 2pi
            assignin('base','hBar_evs',4.1357e-15/(2*pi)); % [eVs] --> Planck's Constant / 2pi
            assignin('base','eps0',8.854e-12); % [C^2/Nm^2] --> Vacuum permitivity
            assignin('base','mu0',4*pi*10^-7); % [Tm/A] --> Magnetic Constant/Vacuum Permeability
            assignin('base','e',1.602e-19); % [C] --> Elementary Charge
            assignin('base','mElectron_kg',9.11e-31); % [kg]
            assignin('base','mElectron_mev',0.511); % [MeVs]
            assignin('base','mElectron_amu',0.00055); % [amu]
            assignin('base','mProton_kg',1.673e-27); % [kg]
            assignin('base','mProton_mev',938.272); % [MeV]
            assignin('base','mProton_amu',1.00728); % [amu]
            assignin('base','mNeutron_kg',1.6749e-27); % [kg]
            assignin('base','mNeutron_mev',939.565); % [MeV]
            assignin('base','mNeutron_amu',1.008665); % [amu]
            assignin('base','amu_kg',1.6605e-27); % [kg]
            assignin('base','amu_mev',939.5654); % [MeV]
            assignin('base','avogadro',6.022e23); % [1/mol]
            assignin('base','k',1.3806e-23); % [J/K] --> Boltzmann's Constant
            assignin('base','R',8.3145); % [J/molK] --> Gas Constant
            assignin('base','sigma',5.670367); % [W/m^2K^4] --> Stefan-Boltzmann Constant
            assignin('base','beta_mmK',2.8978); % [mm/K] --> Wien Displacement Constant
            assignin('base','beta_ghzk',58.789); % [GHz/K] --> Wien Displacement Constant
            assignin('base','atm',101.3e3); % [Pa] --> Standard Atmospheric Pressure
            assignin('base','H0_s',2.25e-18); % [s^-1] --> Hubble Constant
            assignin('base','H0',69.3); % [km/s/Mpc] --> Hubble Constant
            assignin('base','planck_length',1.616229e-35); % [m]
            assignin('base','planck_mass',2.176470e-8); % [kg]
            assignin('base','planck_time',1.616229e-35/3e8); % [s]
            assignin('base','mu_Bohr',9.274e-24); % [J/T] --> Bohr Magneton
            assignin('base','rBohr',5.291772e-11); % [m] --> Bohr Radius
            assignin('base','rElectron',2.8179e-15); % [m]
            assignin('base','fine_structure',7.297352e-3); % Unitless
            assignin('base','kb_ev',8.617e-5); % [eV/K] --> Boltzmann Constant
            assignin('base','rSun_m',695.508e6); % [m]
        end
    end
end

