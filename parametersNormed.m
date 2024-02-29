function param = parametersNormed(with_tF,mass)
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
AU = 1.49597870691e8; % Astronomical unit in [km]
Cr = 0.4; % reflectivity coefficient
A = 10; % equivalent surface of satellite, [m2]
massSC = mass; % s/c mass, [kg]
AMratio = 0.015;
Psr = 4.6e-6; % N/m2
Tmax = 50; % [N]
g0 = 9.81; 
Isp = 400; % [s] 


% NORMALIZING FACTORS
VEL = sqrt(muSun/AU); % [km/s]
TIME = AU/VEL; % [s]
MASS = massSC; 
param.VEL = VEL;
param.TIME = TIME; 
param.MASS = MASS;
param.AU = AU;
% Build struct for parameters
param.muSun = muSun*TIME*TIME/(AU^3);
param.muEarth = muEarth*TIME*TIME/(AU^3);
param.muMoon = muMoon*TIME*TIME/(AU^3);
param.muApo = muApo*TIME*TIME/(AU^3);
param.G = G;
param.P0 = P0*TIME/MASS;
param.Psr = Psr*1e3*TIME*TIME*AU/MASS;
% param.c = c;
param.AU = AU;
param.Cr = Cr;
param.A = A/((AU*1e3)^2);
param.massSC = massSC/MASS;
% param.AMratio = AMratio;
param.J2 = 'yes';
param.J2 = 'no';
param.Tmax = Tmax*(1e-3)*TIME*TIME/MASS/AU; % [N]
param.Isp = Isp/TIME; % s
param.g0 = g0*(1e-3)*TIME*TIME/AU;
param.with_tF = with_tF;

end

