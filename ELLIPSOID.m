% ELLIPSOID AND GRAVITATIONAL FIELD
clc
clear, close all
IbIcRatio = 0.96;
IaIcRatio = 0.64;
mass = 4.3*1e10;
a = 350; %m 


[a, b, c] = semiaxis(a,IbIcRatio, IaIcRatio,mass);

fprintf('Semi-axes: a = %.2f, b = %.2f, c = %.2f\n', a, b, c);





function [a, b, c] = semiaxis(a,Ib_Ic, Ia_Ic, mass)
    % Ib_Ic: Ratio of I_b to I_c
    % Ia_Ic: Ratio of I_a to I_c
    % mass: Mass of the ellipsoid
    
    % Check for valid input
    if Ib_Ic <= 0 || Ia_Ic <= 0 || mass <= 0
        error('Ratios and mass must be positive.');
    end

    % Solve for semi-axes using the given ratios and mass
    % The equations are derived from the formulas for moments of inertia
    % Ia = (1/5) * m * (b^2 + c^2)
    % Ib = (1/5) * m * (a^2 + c^2)
    % Ic = (1/5) * m * (a^2 + b^2)

    % Calculate semi-axes
    
    k1 = Ia_Ic/Ib_Ic;
    k2 = (k1+Ia_Ic*(k1-1))/(k1-Ia_Ic*(k1-1));
    
    % Assuming a as the known parameter
    b = sqrt(k2*a*a);
    c = sqrt((b*b-k1*a*a)/(k1-1));
end
