% EARTH SPHERICAL HARMONICS SNIPPET
% The purpose of this script is to recover the gravitational field of Earth
% accounting for the spherical harmonics contribution up to a degree
% selected by the user.

% Ideally the workflow should include:
%     1) Build/import the zonal harmonics values in an Earth fixed reference frame 
%         * specify the used R.F;
%         * confirm that in such R.F. the contributions are constant;
%     2) Derive/import how this RF moves across the epochs of CPO
%     3) Provide the obtained gravity contribution as RHS of the EOM:
%         * centered in APOPHIS
%         * accounting for Earth spin on its axis

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
state = [0,0, 35000];
% tesseralHarmonicsAcceleration(0,state, 398600)

% Set up a test scenario (you may replace this with actual satellite data)
t = 0; % time (not used in this simple example)
r = [10000, 0, 0]; % [km] example position vector in ECEF frame 
muEarth = 398600.440; % [km^3/s^2] Earth's gravitational parameter

accMain = -muEarth/(norm(r)^3)
% Call the tesseralHarmonicsAcceleration function
[aJ2,aJ3, totalAcc] = J2J3(0, r, 1);

% Display the result
disp('Gravitational acceleration in ECEF frame. J2 contribute:');
disp(aJ2);

disp('Gravitational acceleration in ECEF frame. J3 contribute:');
disp(aJ3);

disp('Gravitational acceleration in ECEF frame. Total contribute:');
disp(totalAcc);
accMain+totalAcc
function [aJ2,aJ3,aTH] = J2J3(t, r, logicJ3)
    % Inputs:
    %   t: time (not used in this simple example)
    %   r: position vector of the satellite in the ECEF (Earth-Centered Earth-Fixed) frame
    %   mu: gravitational parameter of the Earth

    % Constants
    J2 = 1.75553e10; % [km^5/s^2]
    J3 = -2.61913e11; % [km^6/s^2]
    R = norm(r);
    
    x= r(1);
    y = r(2);
    z= r(3);
    
%     aJ2 = zeros(3,1);
    aJ3 = zeros(3,1);
    % Compute contribution for J2 term
    
    aJ2 = [J2*x/(R^7)*(6*z^2-1.5*(x^2+y^2));
       J2*y/(R^7)*(6*z^2-1.5*(x^2+y^2));
       J2*z/(R^7)*(3*z^2-4.5*(x^2+y^2))];

   if logicJ3 == 1
    aJ3 = [J3*x*z/(R^9)*(10*z^2-15/2*(x^2+y^2));
       J3*y*z/(R^9)*(10*z^2-15/2*(x^2+y^2));
       J3/(R^9)*(4*z^2*(z^2-3*(x^2+y^2))+1.5*(x^2+y^2)^2)];
   end
   aTH = aJ2+aJ3;
end
