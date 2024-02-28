function [derivative] = derApophis(r,muApo, main, dir)
% PURPOSE: compute the lambda_dot RHS contribution of the SRP accoerding to
%           derivation variable (1 = x, 2 = y, 3 = z)

% INPUT: r, position s/C wrt apophis
%        muAPo = constant planetary constant for Apophis
%        main, specify main derivation
%        dir, specify component of derivation 

        nr = norm(r);

        if main == dir
            derivative = -muApo*(nr^2-3*r(main)^2)/(nr^5);
        else
            derivative = 3*muApo*r(dir)*r(main)/(nr^5);
        end
end

