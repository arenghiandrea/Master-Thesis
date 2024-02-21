function param = parameters(with_tF)
    %% BUILD PARAMETER STRUCT
G = 6.67430e-20;        % Universal Gravity Constant, [km3/kg/s2]

muSun = cspice_bodvrd('SUN','GM',1); %[km^3/s^2]
muEarth = cspice_bodvrd('EARTH','GM',1); %[km^3/s^2]
muMoon = cspice_bodvrd('MOON','GM',1); %[km^3/s^2]

massApo = 4.3e10;      % Apophis' mass, [kg] from 'Heliotropic orbits...'
muApo = massApo*G;      % Apophis gravitational constant

% SRP parameters
P0 = 1367; % [W/m2]
c = 2.998e8; % speed of light, [m/s] NOTE: leave in meters/s otherwise SRP wrong
AU = 1.495978707e8; % Astronomical unit in [km]
Cr = 0.4; % reflectivity coefficient
A = 1; % equivalent surface of satellite, [m2]
massSC = 1000; % s/c mass, [kg]
AMratio = 0.015;
Psr = 4.6e-6; % N/m2
Tmax = 50; % [N]
g0 = 9.81; 
Isp = 500; % [s] 


% NORMALIZING FACTORS
VEL = sqtrt(muSun/AU); % [km/s]
TIME = AU/VEL; % [s]

MASS = massSC; 

% Build struct for parameters
param.muSun = muSun;
param.muEarth = muEarth;
param.muMoon = muMoon;
param.muApo = muApo;
param.G = G;
param.P0 = P0;
param.Psr = Psr;
param.c = c;
param.AU = AU;
param.Cr = Cr;
param.A = A;
param.massSC = massSC;
param.AMratio = AMratio;
param.J2 = 'yes';
param.J2 = 'no';
param.Tmax = Tmax; % [N]
param.Isp = Isp; % s
param.g0 = g0;
param.with_tF = with_tF;

end

