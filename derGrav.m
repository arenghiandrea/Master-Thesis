function [derivative] = derGrav(r,b, mu, main, dir)
% generate component to add to the derivatives

% INPUT: r, s/c position vector (ECLIPJ2000 centered in Apopohis)
%        b, celestial body position vector (ECLIPJ2000 centered in Apopohis)  
%        mu, constant of the body
%        dir, defines the direction (x,y,z) of derivation as (1,2,3)


    rel = r-b;
    nrel = norm(rel);
    nr = norm(r);
    nb = norm(b);

    SQRT = sqrt(nb^2+nr^2-2*dot(r,b));
    
    if main == dir
        derivative = -mu/(nrel^5)*...
            (nrel^2 -3*(rel(dir))^2+3*b(dir)*rel(dir)/(nb^3)*(nrel^3-nb^2-nr^2+2*dot(r,b)));
%         derivative = -mu/(nrel^5)*(nrel^2
%         -3*(rel(dir))^2+3*b(dir)*rel(dir)/(nb^3)*(nrel^3-SQRT^2));CHECK
%         FOR CORRECTNESS BEFORE IMPLEMENTING
    else
        % example main = 1 (x) --> d/dx
        %         dir = 2 (y) or 3 (z)
        
        % build case 
        derivative = -3*mu*rel(main)/(nrel^5)*...
          (-r(dir)+b(dir)*(3*SQRT*(nrel^2)/(nb^3)-fEnckeDirect(r,b)));
    end
    
    end

