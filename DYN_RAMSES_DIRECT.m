function [dY,acc] = DYN_RAMSES_DIRECT(t,state,param)
    % PURPOSE: define multi-body problem dynamics controlled with low thrust
    %          propulsion. The vector includes: spacecraft state (7x1),
    %          costate (7x1 lagrange multipliers) --> 14x1 
    
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
%     G = param.G;
    P0 = param.P0;
    c = param.c;
    AU = param.AU;
    Cr = param.Cr;
%     Psr = param.Psr;
%     A = param.A;
%     massSC = param.massSC;
    AMratio = param.AMratio;
%     J2 = param.J2;
    
    % Thrusting parameters
    Tmax = param.Tmax;
    Isp = param.Isp;
    g0 = param.g0;
    
%   State
    % Extract state elements
    r = state(1:3);
    nr = norm(r);
    
%     x  = state(1);
%     y  = state(2);
%     z  = state(3);
    u  = state(4);
    v  = state(5);
    w  = state(6);
    m = state(7); % mass
    
    % extract lambda values
    lambdaR = state(8:10);
    lambdaV = state(11:13);
    lambdaM = state(14);
    
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
    
    
    %% Compute phase angle
%     phaseAngle = acos(dot(r,rSun)/(nr*norm(rSun))); % radians
%     phaseAngle = rad2deg(phaseAngle); % degrees
    
    
    %% PERURBATIONS (accelerations)
    
    % Apophis gravity field
    accApo = -muApo*r/(nr^3); % gravity field of Apophis [km/s2]
    
   % Sun gravity field
    accSun = -muSun/((norm(r-rSun))^3)*(r+fEnckeDirect(r,rSun)*rSun); % [km/s2]
    % Earth gravity field
    accEarth = -muEarth/((norm(r-rEarth))^3)*(r+fEnckeDirect(r,rEarth)*rEarth); % [km/s2]
    % Moon gravity field
    accMoon = -muMoon/((norm(r-rMoon))^3)*(r+fEnckeDirect(r,rMoon)*rMoon); % [km/s2]
    
    % Solar Radiation Pressure
    KSRP = P0*AU^2*Cr*AMratio/c*1e-3;
    accSRP = P0*AU^2*Cr*AMratio*d/(c*nd^3)*1e-3; % [3x1] last coefficient to obtain [km/s3]
    
    %% check consistency
%     accSRP = KSRP*d/(nd^3);
    
    % Earth J2 effect
    % Build relative vector
    rJ2 = r+apophis-earth; % [km]
    nRJ2000 = norm(rJ2);
    R_earth = 6378.137; % [km]
    J2 = 0.00108263; % [-]
    kJ2 = 1.5*J2*muEarth*R_earth^2/nRJ2000^4;
    
    accJ2 = [kJ2*rJ2(1)/nRJ2000*(5*rJ2(3)^2/(nRJ2000^2)-1);
            kJ2*rJ2(2)/nRJ2000*(5*rJ2(3)^2/(nRJ2000^2)-1);
            kJ2*rJ2(3)/nRJ2000*(5*rJ2(3)^2/(nRJ2000^2)-3)];
    
%         accJ2 =kJ2*(1-5*RJ2000(3)^2/nRJ2000^2)*RJ2000/nRJ2000;
        
    totAcc = accApo+accSun+accEarth+accMoon+accSRP+accJ2;
  
%     
%     %% testing if the accelerations are included / IN SCOPE
%     totAcc = zeros(3,1);
%     %% testing

    acc = [norm(accApo), norm(accSun), norm(accEarth), norm(accMoon), norm(accSRP), norm(accJ2)];
    
    
    % INCLUDE SWITCHING FUNCTION CONTROL HERE!!!
    SF = - norm(lambdaV)*Isp*g0/m-lambdaM;
    if SF < 0
        
        % LATER INSERT FURTHER CONTROL OVER STATE/PHASE ANGLE
%         if nr <1 || nr>9.5 || phaseAngle<20 || phaseAngle>90 %range in [km] // later insert other conditions
            throttle = 1;
    else
            throttle = 0;
%         end
    end
    
    
    RHS_lambda = zeros(7,1);
    % lambdaDOT component by component
    for kkk = 1:3
        RHS_lambda(kkk) = lambdaV(1)*(derApophis(r,muApo, kkk, 1)+... % apophis
                            +derGrav(r,rSun, muSun, kkk, 1)+... % sun
                            +derGrav(r,rEarth, muEarth, kkk, 1)+... % earth
                            +derGrav(r,rMoon, muMoon, kkk, 1)+... % moon
                            +derSRP(d,KSRP, kkk, 1)+... % SRP
                            +derJ2(rJ2,kJ2,kkk,1))+...
            lambdaV(2)*(derApophis(r,muApo, kkk, 2)+... % apophis
                            +derGrav(r,rSun, muSun, kkk, 2)+... % sun
                            +derGrav(r,rEarth, muEarth, kkk, 2)+... % earth
                            +derGrav(r,rMoon, muMoon, kkk, 2)+... % moon
                            +derSRP(d,KSRP, kkk, 2)+... % SRP
                            +derJ2(rJ2,kJ2,kkk,2))+...
            lambdaV(3)*(derApophis(r,muApo, kkk, 3)+... % apophis
                            +derGrav(r,rSun, muSun, kkk, 3)+... % sun
                            +derGrav(r,rEarth, muEarth, kkk, 3)+... % earth
                            +derGrav(r,rMoon, muMoon, kkk, 3)+... % moon
                            +derSRP(d,KSRP, kkk, 3)+... % SRP
                            +derJ2(rJ2,kJ2,kkk,3));
    end
    
    for kkkk = 4:6
       RHS_lambda(kkkk) = lambdaR(kkkk-3); 
    end
    
    RHS_lambda(7) = norm(lambdaV)*Tmax*throttle/(m^2);
      %% RICORDA DI CAMBIARE IL SEGNO A TUTTA LA RHS DI LAMBDA_DOT COME DA EQUAZIONE (LDOT = - dH/dx)
        % prima scrivi le semplici derivate, assembla e cambia il segno poi
    RHS_lambda = -RHS_lambda;
    
    
    % assemble RHS
    dY = zeros(14,1);
    dY(1:7,1) = [
         u;
        v;
        w;
        totAcc(1)- throttle/m*Tmax*lambdaV(1)/norm(lambdaV);
        totAcc(2)- throttle/m*Tmax*lambdaV(2)/norm(lambdaV);
        totAcc(3)- throttle/m*Tmax*lambdaV(3)/norm(lambdaV);
        -throttle*Tmax/Isp/g0;
    ];

    dY(8:14) = RHS_lambda(:,1);
    
end