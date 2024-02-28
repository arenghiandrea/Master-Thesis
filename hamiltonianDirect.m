function [H, hAve] = hamiltonianDirect(state,param)
    % PURPOSE: compute the Hamiltonian funcitonal for each instant of the TOF
    % INPUT: - state, [14x1] [position (3), velocity (3), mass, lambdas (7)]
    %        - param,  includes several parameters (mass, AU, VEL, Tmax, N, Isp, g0..)
    % OUTPUT: - H [nx1], values of Hamiltonian functional for each instant
    %           of time;
    %         - hAve [1x1], avarage value of Hamiltonian functional
    L = max(size(state));

    H = zeros(L,1);
    for i = 1:L
        pos = state(i,1:3);
        vel = state(i,4:6);
        mass = state(i,7);
        lamR = state(i,8:10);
        lamV = state(i,11:13);
        lamM = state(i,14);
        H(i) = 1 + dot(lamR,vel)-param.muSun/(norm(pos)^3)*dot(lamV,pos)+ param.Tmax/param.Isp/param.g0*param.N*(-norm(lamV)*param.Isp*param.g0/mass-lamM);
    end
    hAve = mean(H);
end