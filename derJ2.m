function [derivative] = derJ2(rel,kJ2,main,dir)
% generate component to add to the derivatives

% INPUT: r, s/c position vector (ECLIPJ2000 centered in Apopohis)
%        b, celestial body position vector (ECLIPJ2000 centered in Apopohis)  
%        mu, constant of the body
%        dir, defines the direction (x,y,z) of derivation as (1,2,3)

    if dir == 3
        N = 3;
    else
        N =1;

    end
    nj = norm(rel);

    if (main == 1 || main == 2) && dir == main

        derivative = kJ2*(...
            5*(rel(3)^2)/(nj^5)*((nj^2)-3*(rel(main)^2))+...
            -(nj-rel(main)^2)/(nj^2));
    elseif (main == 1 || main == 2) && dir ~= main
        derivative = kJ2*rel(dir)*(5*rel(3)^2*(-3*rel(main)/(nj^5))+N*rel(main)/nj);
    elseif main ==3 && dir == main
        derivative = 3*kJ2/(nj^2)*(5/(nj^3)*((nj^2)*rel(3)^2-rel(3)^4)-nj+rel(3)^2);
    else
        derivative = kJ2*rel(dir)*(5/(nj^5)*(2*rel(main)*nj^2-3*rel(main)^3)+rel(main)/nj);
    end
end

