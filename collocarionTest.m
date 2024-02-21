%% testing snippet
clc
clear
close all

cspice_furnsh('MasterThesis.tm');

m0 = 1000; % initial s/c mass, [kg]
mF = 900; % estimate of final mass

% grid intervals
nGrid = 70; % nGrid +1 nodes
nNodes = nGrid+1;
var = 9;

  t0s = '2029-02-08-00:00:00.000 UTC';
  tFs = '2029-02-20-23:59:59.000 UTC';
  
  t0 = cspice_str2et(t0s); %[s]
  tF = cspice_str2et(tFs); %[s]
    
    
  h = (tF-t0)/nGrid; % Time step of discretization, OK
  timeVect = (t0:h:tF)'; % OK, vertical time vector (0) = initial time, t(end) = final time
% itial guess shall include
% The full vector include all the points of the state (x, )

% BUILD INITIAL GUESS VECTOR (9XnNodes,1)
%   position bounds: 0.5 - 10 km
    posLB = 0.5; posUB = 10;
%   velocities bounds: 0.1 m/s - 1 m/s
velLB = 1e-4; velUB = 1e-3; %[km/s]
%   mass bounds: linspace()
mass0 = linspace(m0,mF,nNodes)';
% angles will be randomly selected from 0 to 2*pi radians




% Initialize guess vector
with_tF = 1;
switch with_tF
    case 1
        x0 = zeros(var*nNodes+1,1);
        x0(end) = tF;
    case 0
        x0 = zeros(var*nNodes,1);
end
param = parameters(with_tF);
for k = 1:nNodes
    x0((var*(k-1)+1):(var*(k-1)+3),1) = posLB + (posUB-posLB).*rand(3,1); % positions
    x0((var*(k-1)+4):(var*(k-1)+6),1) = velLB + (velUB-velLB).*rand(3,1); % velocities
    x0((var*(k-1)+7),1) = mass0(k);                                        % mass
    x0((var*(k-1)+8):var*k,1) = 1.9 *pi.*rand(2,1);                          % angles (rad)
end

% x0 = [linspace(0.8,11, 3*nNodes)'; ...% positions [km]
%       linspace(1,50,3*nNodes)'*1e-3; ... % velocitiees
%       linspace(1000,900, nNodes)';... % mass [kg]
% %       ones(nGrid,1)*0.5; ... % control module
%       rand(nNodes,1); ... % angle between r projection over xy-plane and x vector;
%       rand(nNodes,1);... % r-(xy)plane angle;
%       t0+24*3600 %1-day TOF
%       ];
x0;

%%

% inline (auxiliary) functions
% DEFINE OBJECTIVE FUNCTION (FUEL, TOF,...)
fMinimize = @(X) OBJFUN(X, param);

% DEFINE NON LINEAR CONSTRAINT

% [ahah, hihi] = nonlinearConst(x0, t0, param, nGrid,s<tateStart, stateTarget);
% linear constraint: say final time should be larger than intial time
% A = [0,0,0,0,0,0,0, -1,1, 0,0];
A = [];
b = [];
ones(6,1)
lowElem = [-Inf*ones(6,1); mF; 0; 0];
upElem = [+Inf*ones(6,1); m0; 2*pi;2*pi];
lb = []; ub = [];
for el = 1:nNodes
    lb = [lb; lowElem];
    ub = [ub; upElem];
end

if with_tF ==1
    % define lower bound for variables
    lb = [lb; t0+1]; % thrust control module and direction
    % define upper bound for variable
    ub = [ub; tF]; % thrust control module and direction
end
% lb = [];
% ub = [];
optionsFMC = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
    'CheckGradient',true,...
    'MaxIterations',1000,...
    'Display','iter',...
    'Algorithm','active-set',...
    'MaxFunctionEvaluations',1e6);

% % 
% % % test subfunctions: OK
% % 
% % for kk = 1:nNodes
% %     timeVect(kk)
% %     from = var*(kk-1)+1;
% %     to = var*kk;
% %     DYN_RAMSES(timeVect(kk),x0(from:to,1),param)
% % %     [dy(:,1),~] = DYN_RAMSES(timeVect(kk),x0(var*(kk-1):var*kk,1),param)
% % end

% % [ahah, hihi] = nonlinearConst(x0, t0, param, nGrid,stateStart, stateTarget);
% define starting state of the s/c in ECLIPJ200 centered in Apophis
stateStart = [5;5;5]; % [km]
stateTarget = [5;-5; 5; 1.5e-3; 1.5e-3;1.5e-3]; % [km]

fNLCon = @(X) nonlinearConst(X, t0, param, nGrid,stateStart, stateTarget);



%% NLP SOLVER

% state is Y = {x,y,z, u,v,w, m, t0,tF, mod,alpha}
[~,fval,~, ~] = fmincon(fMinimize,x0, ... 
    A,b,... % linear inequalities OK
    [],[],... % linear equalities OK
    lb,ub,... & lower/upper bound OK
    fNLCon,...
    optionsFMC); 

%% FUNCTIONS


function [dY,acc] = DYN_RAMSES(t,state,param)
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
    alpha= state(8); % alpha angle
    beta = state(9); %beta angle
    
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
    
    if nr <1 || nr>9.5 %range in [km] // later insert other conditions
        throttle = 1;
    else
        throttle = 0;
    end
    
    dY = [
        u;
        v;
        w;
        totAcc(1)+throttle*Tmax/m*cos(alpha)*cos(beta);
        totAcc(2)+throttle*Tmax/m*sin(alpha)*cos(beta);
        totAcc(3)++throttle*Tmax/m*sin(beta);
        -throttle*Tmax/Isp/g0];
end

function [obj,grad] = OBJFUN(Y, param)
% PURPOSE: define the objective function (performance index)to be minimized
%          in the NLP problem

% INPUT:   Y (10*nGrid+1,1)--> - positions at mesh point (3*nGrid),
%                              - velocities at mesh point (3*nGrid),
%                              - mass at mesh point (3*nGrid),

%                              - contol module (nGrid)
%                              - control angle alpha (nGrid)
%                              - control angle beta (nGrid)
%                              - final time (1)
   

% Define output: performance index

% Time optimal problem
    if param.with_tF ==1
        obj = Y(end);
        if nargout >1
            grad = zeros(length(Y)-1,1);
            grad = [grad;1];
        end
    else 
        obj = 1/Y(end-2,1);
    end
end

function [cIneq,c] = nonlinearConst(Y, t0, param, nInterv,stateStart, stateTarget)
    % PURPOSE: In the context of direct methods/collocation the dynamics
    %           equations are trasnformed into constraints between successive nodes
    %           (Euler or Hermite-Simpson schemes). The purpose of this function is
    %           to define such constraints and all those associated to path and
    %           boundary constraints.
    
    % INPUT: Y, meshed grid states and controls (9*nGrid)
    %        param, generic parameters useful for the problem considered
    
    % OUTPUT: cIneq: inequality constraints in the shape cIneq <0
    %          cEq: equality constraint
    
    %___________________%
    
    % Extract variables
    tF = Y(end);
    % define time step
    h = (tF-t0)/nInterv;
    nodes = nInterv+1;
    % tune this parameter according to the number of variables included in
    % the vector: 3 positions, 3 velocities, 1 mass, 2 angles
    var = 9;
    t = linspace(t0,tF,nodes);
%     nodes = length(t);
    cIneq = [];
    c = [];
    dy = zeros((var-2)*nodes,1);
    
     
    for kk = 1: nodes
%         fprintf('Inside here at time instant %1.f \n',kk)
       from = var*(kk-1)+1;
       to = var*kk;
       [dy((var-2)*(kk-1)+1:(var-2)*kk,1),~]= DYN_RAMSES(t(kk),Y(from:to,1),param);
%        [dy(from:to,1),~] = DYN_RAMSES(t(kk),Y(from:to,1),param)
    end
    % Euler method to itroduce dynamics equation
    for i = 1:(nodes-1)
        c(((i-1)*(var-2)+1):(var-2)*i,1) = ...
            Y(i*var+1:i*var+(var-2),1)-Y((i-1)*var+1:(i-1)*var+(var-2),1)-h*dy((var-2)*(i-1)+1:(var-2)*(i-1)+(var-2),1);
    end
    size(c);
    
    % add constraints on initial condition
    c = [c;
        Y(1:3,1)-stateStart(1:3,1);
        ];
    % add constraints on target condition
    if param.with_tF == 1
        
        c= [c;...
            Y(end-9:end-4,1)-stateTarget(1:6,1)];
    else
        c= [c;...
            Y(end-8:end-3,1)-stateTarget(1:6,1)];
    end
end