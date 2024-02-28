function [dState,acc] = closeProxODEseparate(t,state,param)
    % PURPOSE: define multi-body problem dynamics controlled with low thrust
    %          propulsion
    
    % INPUTS:
        %
        %
        %
        %
    % OUTPUTS: derivatives of state vector dState
    
    % Extract constants for easier use
    muSun = param.muSun;
    muEarth = param.muEarth;
    muMoon = param.muMoon;
    muApo = param.muApo;
    G = param.G;
    P0 = param.P0;
    c = param.c;
    AU = param.AU;
    Cr = param.Cr;
    Psr = param.Psr;
    A = param.A;
    massSC = param.massSC;
    AMratio = param.AMratio;
    J2 = param.J2;
    
%   State
    % Extract state elements
    r = state(1:3);
    nr = norm(r);
    
    x  = state(1);
    y  = state(2);
    z  = state(3);
    u  = state(4);
    v  = state(5);
    w  = state(6);
    m = state(7); % s/c mass
    
    % Extract costate
    lambda = state(8:14);
    
    
    % Retrieve ephemeris @(t) wrt SSB (ECLIPJ2000)
    apophis = cspice_spkpos('20099942',t,'ECLIPJ2000','NONE','10'); %[km; km/s] rAS
    earth = cspice_spkpos('Earth',t,'ECLIPJ2000','NONE','10'); %[km; km/s] rES
    moon = cspice_spkpos('Moon',t,'ECLIPJ2000','NONE','10');   %[km; km/s] rMS
    
    
    % Relative vectors: r_j
    rSun = -apophis;
    rEarth = earth-apophis;
    rMoon = moon - apophis;
    d = r+apophis;
    nd = norm(d);
    
    % Compute the perturbing acceleration contribution
    accApo = -muApo*r/(nr^3); % gravity field of Apophis [km/s2]
    
    % Perform computation for Encke method
    qSun = (nr^2)/(norm(rSun)^2)-2*dot(r,rSun)/(norm(rSun)^2);
    qEarth = (nr^2)/(norm(rEarth)^2)-2*dot(r,rEarth)/(norm(rEarth)^2);
    qMoon = (nr^2)/(norm(rMoon)^2)-2*dot(r,rMoon)/(norm(rMoon)^2);
    
    accSun = -muSun/((norm(r-rSun))^3)*(r+fEncke(qSun)*rSun); % [km/s2]
    accEarth = -muEarth/((norm(r-rEarth))^3)*(r+fEncke(qEarth)*rEarth); % [km/s2]
    accMoon = -muMoon/((norm(r-rMoon))^3)*(r+fEncke(qMoon)*rMoon); % [km/s2]
    accSRP = P0*AU^2*Cr*AMratio*d/(c*nd^3)*1e-3; % last coefficient to obtain [km/s3]
    
    % Earth J2 effect
    % Build relative vector
    RJ2000 = r+apophis-earth; % [km]
    nRJ2000 = norm(RJ2000);
    R_earth = 6378.137; % [km]
    J2 = 0.00108263; % [-]
    kJ2 = 1.5*J2*muEarth*R_earth^2/nRJ2000^4;
    
    accJ2 = [kJ2*RJ2000(1)/nRJ2000*(5*RJ2000(3)^2/(nRJ2000^2)-1);
            kJ2*RJ2000(2)/nRJ2000*(5*RJ2000(3)^2/(nRJ2000^2)-1);
            kJ2*RJ2000(3)/nRJ2000*(5*RJ2000(3)^2/(nRJ2000^2)-3)];
    
%         accJ2 =kJ2*(1-5*RJ2000(3)^2/nRJ2000^2)*RJ2000/nRJ2000;
        
    totAcc = accApo+accSun+accEarth+accMoon+accSRP+accJ2;

    acc = [norm(accApo), norm(accSun), norm(accEarth), norm(accMoon), norm(accSRP), norm(accJ2)];
    
    dState = [
        u;
        v;
        w;
        totAcc(1);
        totAcc(2);
        totAcc(3)];