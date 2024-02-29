function [derivative] = derSRP(d,KSRP,m, main, dir)
% PURPOSE: compute the lambda_dot RHS contribution of the SRP accoerding to
%           derivation variable (1 = x, 2 = y, 3 = z)

% INPUT: d = r+apophis
%        KSRP = constant terms factor
%        s/c mass
%        main, specify main derivation
%        dir, specify component of derivation 

        nd = norm(d);

        if main == dir
            derivative = KSRP/(nd^5)*(nd^2-3*(d(main)^2));
        else
            derivative = KSRP*d(dir)*(-3*d(main)/(nd^5));
        end
end

