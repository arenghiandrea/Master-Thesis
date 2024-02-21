%% testing snippet
clc
clear


Y = [1,2,3,4,5,6,7,8,9,10,11]';

m0 = 1000; % initial s/c mass, [kg]
n = length(Y);
PIX = 'deltav';

S = objectiveFun(Y,PIX);


% itial guess shall include
%  x,y,z,u,v,w,m,t0,tF  (9x1)
xGuess = [X0; tIn; tF]; 

param = 5;

% inline (auxiliary) functions
% DEFINE OBJECTIVE FUNCTION (FUEL, TOF,...)
fMinimize = @(X) OBJFUN(X, param);

% DEFINE NON LINEAR CONSTRAINT
fNLCon = @(X) nonlinearConst(X, param);

% linear constraint: say final time should be larger than intial time
A = [0,0,0,0,0,0,0, -1,1, 0,0];
b = zeros(n,1);

% define lower bound for variables
lb = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf,... % s/c state
    600,... % mass, [kg]
    initEpoch, initiEpoch+1, ... % initial time (as initial epoch in the unit of measure proper as input for ephemeris)
    0,0]; % thrust control module and direction
% define upper bound for variables
ub = [+Inf, +Inf, +Inf, +Inf, +Inf, +Inf,...% s/c state
    m0,... % mass, [kg]
    finalEpoch-1, finalEpoch, ... % initial time (as initial epoch in the unit of measure proper as input for ephemeris)
    1,2*pi]; % thrust control module and direction


optionsFMC = optimoptions('fmincon','SpecifyObjectiveGradient',false, ...
    'SpecifyConstraintGradient',false,...
    'CheckGradient',false,...
    'algorithm','active-set');

% state is Y = {x,y,z, u,v,w, m, t0,tF, mod,alpha}
[~,fval,~, ~] = fmincon(fMinimize,xGuess, ... 
    A,b,... % linear inequalities OK
    [],[],... % linear equalities OK
    lb,ub,... & lower/upper bound OK
    fNLCon,...
    optionsFMC); 

%% FUNCTIONS


function [dY,acc] = DYN_RAMSES(t,Y,param)
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
    
    % Thrusting parameters
    Tmax = param.Tmax;
    
    
%   State
    % Extract state elements
    r = Y(1:3);
    nr = norm(r);
    
    x  = Y(1);
    y  = Y(2);
    z  = Y(3);
    u  = Y(4);
    v  = Y(5);
    w  = Y(6);
    
    m = Y(7);
    t0 = Y(8); tF = Y(9);
    h = Y(10); alpha= Y(11);
    
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
    
    dY = [
        u;
        v;
        w;
        totAcc(1)+h*Tmax/m*alphaX;
        totAcc(2)+h*Tmax/m*alphaY;
        totAcc(3)++h*Tmax/m*alphaZ;
        -h*Tmax/Isp/g0];
end

function [obj] = OBJFUN(Y, param)

    % extract (augmented) state variables
    pos = Y(1:3); % position
    vel = Y(4:6); % velocity
    mass = Y(7); % mass
    t0 = Y(8); % initial time
    tF = Y(9); % final time
    tModule = Y(10); % [0,1] control command
    tAlpha = Y(11); % control direction

    aux = param;
    options = odeset('Abstol', 1e-13,'Reltol', 1e-13);
    
    % insert DYN_ODE integration 
    
    
    
    
    %
    
    obj = tF-t0;
    
    
    % BEGIN WITH A  TIME OPTIMAL
end


function [cIneq,cEq] = nonlinearConst(t, Y, param)
    
  % extract (augmented) state variables
    pos = Y(1:3); % position
    vel = Y(4:6); % velocity
    mass = Y(7); % mass
    t0 = Y(8); % initial time
    tF = Y(9); % final time
    tModule = Y(10); % [0,1] control command
    tAlpha = Y(11); % control direction

% define range as the norm of position vector
range = norm(yOut(1:3));
cIneq = [range - 10;
         -range+0.5
         %add phase angle constraints eventually
         ];
     
     options = odeset('Abstol', 1e-14,'Reltol', 1e-14);
     % perform integration of DYN_ODE
     [out1, dY] = ode113(@DYN_RAMSES, [t0 tF], Yguess, options, param)
     
cEq = []
end