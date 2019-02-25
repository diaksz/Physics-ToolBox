%% Astro HW #3

%% Problem 1
alphaSun = findAlpha(7e5,physics_convert('AU2km',1)); 

radius_sun = 695508; %km
d_bg_km = physics_convert('pc2km',(10^((.8 + 5.5 - 5)/5)));
alphaBetelgeuse = findAlpha(650*radius_sun, d_bg_km)

alphaM31 = findAlpha(physics_convert('pc2km',30e3), physics_convert('pc2km',.7e6))

alphaComaCluster = findAlpha(physics_convert('pc2km',3e6), physics_convert('pc2km',100e6))

%% Problem 2
% IDEA: magnitudes don't add, but luminosities can
% Page 311 from the textbook shows that each magnitude is 10^.4 times the
% previous. This implies that 10^(.4*mtot) = 10^(.4*m1) + 10^(.4*m2) 

m1 = 12.5;
m2 = 12.9;
mTotal = -2.5*log(10^(-.4*m1) + 10^(-.4*m2))

%% Problem 3
% I think abs magnitudes add the same way?

MTotal = -2.5*log(100*10^(-.4*0) + 1000*10^(-.4*3) + 10000*10^(-.4*6))

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

MTotal_Binary = findMBinary(a_km,period_sec) % In kg

%% Problem 5 - check units
% Find total mass of binary system
 
period = 3.96; % days
vA = 108e3; % m/s
vB = 111e3; % m/s

MTotal_Binary_P5 = (physics_convert('day2sec',period)/(2*pi*G)) * (vA + vB)^3 % kg

%% Problem 6
% a) Find ratio of luminosities
% b) Find ratio of surface temperatures

sigma = 5.670367e-8; %W/m^2K^4
RSun = 695508; %km
LSun = 3.86e26; %W
TA = 6530; %K
RA = 2.06*RSun; %km
RB = .0096*RSun; %km
MBolB = 12.9;

LB = LSun * 10^(.4*(4.74-MBolB));
LA = 4*pi*RA^2*sigma*TA^4;

LAtoLB = LA/LB

TB = (LB / (4*pi*RB^2*sigma))^(1/4);

TAtoTB = TA/TB

%% Problem 7

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

