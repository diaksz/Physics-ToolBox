%% Astro HW #3
load('physics_constants.mat') 

f = fopen('hw3_out.txt','w');
%% Problem 1
alphaSun = findAlpha(7e5,physics_convert('AU2km',1));

RSun = 695508; %km
d_bg_km = physics_convert('pc2km',(10^((.8 + 5.5 + 5)/5)));
alphaBetelgeuse = findAlpha(650*RSun, d_bg_km);

alphaM31 = findAlpha(physics_convert('pc2km',30e3), physics_convert('pc2km',.7e6));

alphaComaCluster = findAlpha(physics_convert('pc2km',3e6), physics_convert('pc2km',100e6));

fprintf(f,'----------- PROBLEM 1 --------------\n') 
fprintf(f,'Sun: %8.4f \n',alphaSun)
fprintf(f,'Betelgeuse: %8.4f \n',alphaBetelgeuse)
fprintf(f,'Coma Cluster: %8.4f \n',alphaComaCluster)

%% Problem 2
% IDEA: magnitudes don't add, but luminosities can
% Page 311 from the textbook shows that each magnitude is 10^.4 times the
% previous. This implies that 10^(.4*mtot) = 10^(.4*m1) + 10^(.4*m2) 

m1 = 12.5;
m2 = 12.9;
mTotal = -2.5*log10(10^(-.4*m1) + 10^(-.4*m2));

fprintf(f,'----------- PROBLEM 2 --------------\n')
fprintf(f,'Total Magnitude: %8.4f \n', mTotal)

%% Problem 3
% I think abs magnitudes add the same way?

MTotal_p3 = -2.5*log10(100*10^(-.4*0) + 1000*10^(-.4*3) + 10000*10^(-.4*6));
fprintf(f,'----------- PROBLEM 3 --------------\n')
fprintf(f,'Total Magnitude (absolute): %8.4f \n', MTotal_p3)
%% Problem 4 
% Find total mass of binary system

G = 6.67408e-11;  %Nm^2/kg^2
parallax = .4; % arcsec
maxSeparation = 6; % arcsec
period = 80; % years

d_parsec = 1/parallax;
a_parsec = maxSeparation*d_parsec;
a_km = physics_convert('pc2km',a_parsec);
period_sec = physics_convert('year2sec',period);

MTotal_Binary = findMBinary(1000*a_km,period_sec); % In kg

fprintf(f,'----------- PROBLEM 4 --------------\n')
fprintf(f,'Total Mass: %e kg \n',MTotal_Binary)
%% Problem 5 
% Find total mass of binary system
 
period = 3.96; % days
vA = 108e3; % m/s
vB = 111e3; % m/s

MTotal_Binary_P5 = (physics_convert('day2sec',period)/(2*pi*G)) * (vA + vB)^3 ; % kg

fprintf(f,'----------- PROBLEM 5 --------------\n')
fprintf(f,'Total Mass: %e kg \n',MTotal_Binary_P5)
%% Problem 6
% a) Find ratio of luminosities
% b) Find ratio of surface temperatures

sigma = 5.670367e-8; %W/m^2K^4
RSun = 695508*1000; % m
LSun = 3.86e26; %W
TA = 6530; %K
RA = 2.06*RSun; % m
RB = .0096*RSun; % m
MBolB = 12.9;

LB = LSun * 10^(.4*(4.74-MBolB));
LA = 4*pi*RA^2*sigma*TA^4;

LAtoLB = LA/LB;

TB = (LB / (4*pi*RB^2*sigma))^(1/4);

TAtoTB = TA/TB;

fprintf(f,'----------- PROBLEM 6 --------------\n')
fprintf(f,'Luminosity A:B = %e \n',LAtoLB)
fprintf(f,'Temperature A:B = %8.4f \n',TAtoTB)
%% Problem 7
% Find the distance to 9 Sagittarii

mV = 6;
MV = -5.7;

d_sag = physics_convert('pc2km',10*10^((mV - MV)/5)); % km

fprintf(f,'----------- PROBLEM 7 --------------\n')
fprintf(f,'Distance to 9 Sagittarii: %e km \n',d_sag)
%% Problem 8
% Find the central pressure

mass_sun = 1.989e30; %kg

P_K0V = findPressure(.8*mass_sun, .85*RSun); % N/m^2
P_K0III = findPressure(4*mass_sun, 16*RSun); % N/m^2
P_K0I = findPressure(13*mass_sun, 200*RSun); % N/m^2

fprintf(f,'----------- PROBLEM 8 --------------\n')
fprintf(f,'Central Pressure K0V: %e N/m^2 \n',P_K0V)
fprintf(f,'Central Pressure K0III: %e N/m^2 \n',P_K0III)
fprintf(f,'Central Pressure K0I: %e N/m^2 \n',P_K0I)
%% Problem 9
% see attached

fprintf(f,'----------- PROBLEM 9 --------------\n')
fprintf(f,'See attached \n')
%% Problem 10
% see attached

fprintf(f,'----------- PROBLEM 10 --------------\n')
fprintf(f,'See attached \n')
%% Functions

% used in problem 1
function output = findAlpha(R,d)
    output = 2*atan(R/d);
end

% used in problem 4
function output = findMBinary(a,P)
    G = 6.67408e-11;
    output = (4*pi^2/G)*(a^3 / P^2);
end

% used in problem 8
function output = findPressure(M,R)
    G = 6.67408e-11;
    output = 3*M^2*G / (8*pi*R^4);
end







