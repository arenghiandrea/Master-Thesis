%% SPACECRAFT GUIDANCE AND NAVIGATION
% PROFs: FRANCESCO TOPPUTO, PIERLUIGI DI LIZIA
% A.Y. 2022-23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%   AUTHOR: ANDREA ARENGHI   %
%   PERSONAL CODE: 10523901  %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ASSIGNMENT 1: EX2
% prepare environment
clc
clear
close all
warning('off','all')
%% TASK 1: produce first guess solution

% define four scalar paramters for initial guesses to the problem
alpha = 1.5*pi; % [rad] Angle defined on the Earth circular parking orbit
beta = 1.41; % [-] Initial-to-circular velocity ratio
delta = 7; % [day] Initial time
tIn = 0; % [day] Transfer duration
% compute final time
tF = tIn+delta;

% USEFUL PARAMETERS
mu =0.0121506683; % Earth-Moon mass parameter
omegaSun = -9.25195985e-1; % Scaled angular velocity of the Sun  
DU = 3.84405e5; % [km] - intial circular orbit scaled radius at Earth 
Rearth = 6378; % [km] - Earth's mean radius
Rmoon = 1738; % [km] - Moon's mean radius
hi = 167; % [km]- orbit altitude at Earth
hf = 100; % [km] - orbit altitude at Moon

% COMPUTE OTHER USEFUL VARIABLES
r0 = (Rearth+hi)/DU; % scaled radius at Earth
rF = (Rmoon+hf)/DU; % scaled radius at Moon
v0 = beta*sqrt((1-mu)/r0); % initial velocity aligned

% INITIAL STATE FORMALIZATION
x0 = r0*cos(alpha)-mu;
y0 = r0*sin(alpha);
v0x = -(v0-r0)*sin(alpha);
v0y = (v0-r0)*cos(alpha);
% assemble initial condition vector
X0 = [x0, y0, v0x, v0y]';

% orbit @ Earth
optionsInteg = odeset('Abstol', 1e-14,'Reltol', 1e-14);
a = reshape(eye(6),1,[]);

% integrate dynamics
[tt, intdX] = ode113(@(t,X) PBCR4BP(t,X,mu),[tIn tF],X0, optionsInteg);

% PLOTS
% 1) Synodic reference frame
figure()
subplot(1,2,1)
plot(intdX(:,1),intdX(:,2),'Color','#ff9933','Linewidth',1.4); % Orbit
hold on
plot((1-mu),0,'ok','MarkerFaceColor', '#ffff66', 'MarkerSize',6); % Moon
grid on
plot(-mu,0,'ok','MarkerFaceColor','#00cc66', 'MarkerSize',6); % Earth
% plot(0,0, 'ok','MarkerFaceColor', 'black', 'MarkerSize',1);
plot(intdX(1,1), intdX(1,2),'xk','MarkerSize',8); % Initial point
plot(intdX(end,1), intdX(end,2),'dk','MarkerFaceColor','#cc5200','MarkerSize',7); %Final point
xlabel('x [DU]','FontSize',14);
ylabel('y [DU]','FontSize',14);

title('Initial guess orbit propagation','FontSize',16);
subtitle('Synodic reference frame','FontSize',12)

% 2) Earth-centered Inertial reference frame
% get ECI r.f.
stateECI = inertialTrans(intdX, tt, mu, 'M');

subplot(1,2,2)
plot(stateECI(:,1), stateECI(:,2),'Color','#ff9933','Linewidth',1.4);
hold on
plot(0,0,'ok','MarkerFaceColor','#00cc66', 'MarkerSize',6);% Earth
plot(1,0,'ok','MarkerFaceColor', '#ffff66', 'MarkerSize',6);% Moon
plot(stateECI(1,1), stateECI(1,2),'xk','MarkerSize',8); % Initial point
plot(stateECI(end,1), stateECI(end,2),'dk','MarkerFaceColor','#cc5200','MarkerSize',7); % Final point
xlabel('x [DU]','FontSize',14);
ylabel('y [DU]','FontSize',14);
axis padded
grid on
subtitle('Earth-centered Inertial reference frame','FontSize',12)
legend('Orbit propagation','Earth','Moon','Initial state','Final state','FontSize',13, 'Location', 'best')

%% TASK 2 - SIMPLE SHOOTING PROBLEM
% rename initial guess value
xGuess = [X0; tIn; tF]; 
% xGuess is either 6x1 or 22x1 (in this case we have state(1:4,1), STM(5:20),tIn(22), tF(22)

% define function to be inputted in in fmincon
fMinimize = @(X) costFun(X,r0,rF, mu);

% define function that summarizes the non linear constraints of the problem
fNLCon = @(X) distAndTan(X, r0,rF, mu);

% set fmincon options
optionsFMC = optimoptions('fmincon','SpecifyObjectiveGradient',false, ...
    'SpecifyConstraintGradient',false,...
    'CheckGradient',false,...
    'algorithm','active-set');

% define constraint on time
A = [0 0 0 0 1 -1];

% Case a) GRADIENTS/STM NOT PROVIDED
tic;
% perform fmincon
[~,fval,~, ~] = fmincon(fMinimize,xGuess, ... 
    A,0,... % linear inequalities
    [],[],... % linear equalities
    [-Inf, -Inf, -Inf, -Inf, 0,0],[+Inf,+Inf,+Inf,+Inf, abs(2*pi/omegaSun),23],...
    fNLCon,...
    optionsFMC); 
timeElapsed = toc;
fprintf('The search for optimal Delta V (not providing derivatives) took %3.2f seconds', timeElapsed);

% b) GRADIENTS/STM PROVIDED
% Now the gradients of the cost function and constraint shall be fed to the
% algortihm, after being defined analytically (see FUNCTIONS section)

% redefine options including gradients
optionsFMCgrad = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true,...
    'CheckGradient',false,...
    'algorithm','active-set');
% perform fmincon
tic;
[xMin,fvalGrad,~, ~] = fmincon(fMinimize,xGuess, ...
    A,0,[],[],...
    [-Inf, -Inf, -Inf, -Inf, 0,0],[+Inf,+Inf,+Inf,+Inf, abs(2*pi/omegaSun),23],...
    fNLCon,optionsFMCgrad); 
timeElapsedGrad = toc;
fprintf('The search for optimal Delta(v) (providing derivatives) took %3.2f seconds.\n', timeElapsedGrad);
relErr = abs(fvalGrad-fval)/fvalGrad*100;
fprintf('The relative error between the two cases is %.3f%%.\n', relErr)

% PLOT SOLUTION
% retrieve final state for plotting purposes
[ttFFFF, FINALORBgrad] = ode113(@PBCR4BP,[xMin(5), xMin(6)],xMin(1:4), optionsInteg,mu);
figure()
subplot(2,5,[1,2,3,6,7,8])
plot(FINALORBgrad(:,1),FINALORBgrad(:,2),'Color','#ff9933','Linewidth',1.4);
hold on
grid on
plot((1-mu),0,'ok','MarkerFaceColor', '#ffff66', 'MarkerSize',6); % Moon
plot(-mu,0,'ok','MarkerFaceColor','#00cc66', 'MarkerSize',6); % Earth
circle((1-mu),0, rF,'#ffcc00'); % starting parking orbit
circle(-mu,0, r0,'#669900'); % arrival parking orbit
title('Simple shooting trajectory','FontSize',14)
xlabel('x [DU]','FontSize',12);
ylabel('y [DU]','FontSize',12);
legend('Orbit propagation','Moon','Earth', 'Moon parking orbit','Earth parking orbit','Location','bestoutside','FontSize',13)
subplot(2,5,[4,5])
plot(FINALORBgrad(:,1),FINALORBgrad(:,2),'Color','#ff9933','Linewidth',1.4);
hold on
plot((1-mu),0,'ok','MarkerFaceColor', '#ffff66', 'MarkerSize',6); % Moon
grid on
plot(-mu,0,'ok','MarkerFaceColor','#00cc66', 'MarkerSize',12); % Earth
circle((1-mu),0, rF,'#ffcc00'); % starting parking orbit
circle(-mu,0, r0,'#669900'); % arrival parking orbit
ylim([-0.022 0.042])
xlim([-0.04 0.04])
title('Earth departure','FontSize',14)
xlabel('x [DU]','FontSize',10);
ylabel('y [DU]','FontSize',10);

subplot(2,5,[9,10])
plot(FINALORBgrad(:,1),FINALORBgrad(:,2),'Color','#ff9933','Linewidth',1.4);
hold on
plot((1-mu),0,'ok','MarkerFaceColor', '#ffff66', 'MarkerSize',8); % Moon
grid on
plot(-mu,0,'ok','MarkerFaceColor','#00cc66', 'MarkerSize',6); % Earth
circle((1-mu),0, rF,'#ffcc00'); % starting parking orbit
circle(-mu,0, r0,'#669900'); % arrival parking orbit
ylim([-0.01 0.01])
xlim([0.98 1.001])
title('Moon arrival','FontSize',14)
xlabel('x [DU]','FontSize',10);
ylabel('y [DU]','FontSize',10);

%% TASK 2.3 - MULTIPLE SHOOTING
% number of firings
N = 4;

% initialize state vector
xGuessMult = [];
% initial state
I = reshape(eye(4), 16,1);
xGuessMult(:,1) = [xGuess(1:4,1); I];

% time grid
tGrid = [];
tGrid(1) = tIn;
for iGrid = 2:N
    tGrid(iGrid) = tIn + (iGrid-1)/(N-1)*(tF-tIn);
    tspanMult = [tGrid(iGrid-1), tGrid(iGrid)];
    
    [~, xMult] = ode113(@(t,X) PBCR4BP(t,X,mu), tspanMult, xGuessMult(:,iGrid-1),optionsInteg);
    aux = xMult(end,1:20);
    aux = aux';
    xGuessMult(:,iGrid) = aux;
end

% NLP
% generate variables
stateMulti = reshape(xGuessMult(1:4,:),16,1);
yFmin = [stateMulti; tIn; tF];

% define options
optionsFMCgradMulti = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true,...
    'CheckGradient',false,...
    'algorithm','active-set');

% SET FUNCTIONS
fMulti = @(x) costFunMulti(x,r0,rF, mu);
nonLinMulti = @(X) multiDistAndTan(X, r0,rF, mu,N);

% SOLVE
tic
[xMulti, fMulti, exflag, output] = fmincon(fMulti,yFmin, ...
    [],[],...
    [],[],[],[],...
    nonLinMulti,optionsFMCgradMulti); 
multiToc = toc;
fprintf('The optimal Delta(v) in the multiple shooting case is: %1.5f km/s.\n', fMulti);
fprintf('The elapsed time for the Multiple Shooting was %3.1f seconds.\n', multiToc)

% final solution
XFINALMULTI = reshape(xMulti(1:16),4,4);
tInMulti=xMulti(17);
tFinMulti=xMulti(18);

% PLOT Multiple Shooting solution

optionsOdeM = odeset('Abstol', 1e-14,'Reltol', 1e-14);
tjM(1)= tInMulti;

figure()
for j = 2:N
    tjM(j)=tInMulti+(j-1)/(N-1)*(tFinMulti-tInMulti);
    tspan = [tjM(j-1),tjM(j)];
    [~, xxPlot] = ode113(@PBCR4BP, tspan, XFINALMULTI(:,j-1),optionsInteg,mu);
    plot(xxPlot(:,1), xxPlot(:,2),'LineWidth',2);
    hold on
end
plot((1-mu),0,'ok','MarkerFaceColor', '#ffff66', 'MarkerSize',6); % Moon
grid on
plot(-mu,0,'ok','MarkerFaceColor','#00cc66', 'MarkerSize',8); % Earth
circle((1-mu),0, rF,'#ffcc00'); % starting parking orbit
circle(-mu,0, r0,'#669900'); % arrival parking orbit
plot(xMulti(1), xMulti(2), '*', 'MarkerSize', 10,'Color','r');
plot(xMulti(5), xMulti(6), '*', 'MarkerSize', 10,'Color','r');
plot(xMulti(9), xMulti(10), '*', 'MarkerSize', 10,'Color','r');
plot(xMulti(13), xMulti(14), '*', 'MarkerSize', 10,'Color','r');
xlabel('x [DU]','FontSize',10);
ylabel('y [DU]','FontSize',10);
title('Multiple shooting solution in the rotating frame','FontSize',16)
subtitle('Synodic reference frame','FontSize',14)
legend('Initial leg','Second leg','Final leg','Moon','Earth','','','Shooting',...
    'Location','southoutside',...,
    'Orientation','horizontal','NumColumns',6,...
    'FontSize',14)
%% FUNCTIONS definitions
% DYNAMIC PROBLEM
function dxdt = PBCR4BP(t,X,mu)
% PURPOSE: Define PLANAR BICIRCULAR RESTRICTED 4BP in state form
% formulation
% INPUT: 
%     t, time vector [nx1]
%     X, state vector (4x1), including x, y, u, v
%     mu, mass parameter [scalar]
%     isSTM, logic operator to understand if STM is present
% OUTPUT
%     dxdt, RHS of PBCR4BP[nx4]
% mu =0.0121506683;
% constants
rho = 388.811143; % scaled Sun-E-Moon distance
omegaSun = -0.925195985; %  scaled angular velocity of the Sun
mSun = 3.28900541*1e5; % scaled mass of the Sun

% PLANAR BICIRCULAR RESTRICTED 4BP
x = X(1);
y = X(2);
u = X(3);
v = X(4);

% define auxiliary quantities
r1 = sqrt((x+mu)^2+y^2);
r2 = sqrt((x+mu-1)^2+y^2);
r3 = sqrt((x-rho*cos(omegaSun*t))^2+(y-rho*sin(omegaSun*t))^2);

% % definition of 2D potential
% U3 = 0.5*(x^2+y^2)+(1-mu)/r1+mu/r2+0.5*mu*(1-mu);
% U4 = U3 + mSun/r3-mSun/(rho^2)*(x*cos(omegaSun*t)+y*sin(omegaSun*t));

% define potential U derivatives
U4x = x - (1-mu)*(mu+x)/(r1^3)+mu*(1-mu-x)/(r2^3)...
    -mSun*(x-rho*cos(omegaSun*t))/(r3^3)...
    -mSun/(rho^2)*cos(omegaSun*t);
U4y = y - (1-mu)*y/(r1^3)-mu*y/(r2^3)...
    -mSun/(r3^3)*(y-rho*sin(omegaSun*t))...
    -mSun/(rho^2)*sin(omegaSun*t);

dxdt = zeros(4,1);
% generate differential vector for s/c state
dxdt(1,1) = u;
dxdt(2,1) = v;
dxdt(3,1) = 2*v+U4x; 
dxdt(4,1) = -2*u+U4y; 

%% this section takes care of task 2.2b which requirese the integration of State Transition Matrix
if(length(X) > 4)
    % -- Unpack STM --
    STM = reshape(X(5:20), [4,4]);
    
    % -- Variational equations --
    a31 = 1 - (1-mu)/(r1^3) + 3*(1-mu)*(x+mu)^2/(r1^5) - mu/(r2^3) + 3*mu*(x-1+mu)^2/r2^5 ...
          - mSun/(r3^3) + 3*mSun*(x-rho*cos(omegaSun*t))^2/(r3^5);
    a32 = 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5 ...
          + 3*mSun*(x-rho*cos(omegaSun*t))*(y-rho*sin(omegaSun*t))/r3^5;
    a41 = a32;
    a42 = 1 - (1-mu)/r1^3 + 3*(1-mu)*y^2/r1^5 - mu/r2^3 + 3*mu*y^2/r2^5 ...
          - mSun/r3^3 + 3*mSun*(y-rho*sin(omegaSun*t))^2/r3^5;

    A_stm = [0     0    1   0;
            0     0    0   1;
            a31  a32   0   2;
            a41  a42  -2   0];
    
    dMdt = A_stm*STM;
    diff = [dxdt; reshape(dMdt,16, 1)];
    else
    diff = dxdt;
end
    dxdt = diff;
end

% SIMPLE SHOOTING
function [DELTAV, grad] = costFun(X,r0,rF, mu)
    % PURPOSE: Define DELTAV function to be minimized
    % INPUT: 
    %     t, time vector [nx1]
    %     X, state vector (4x1), including x, y, u, v
    %     mu, mass parameter [scalar]

    % extract variables
    x = X(1);
    y = X(2);
    u = X(3);
    v = X(4);
    t0 = X(5);
    tFin = X(6);
    options = odeset('Abstol', 1e-13,'Reltol', 1e-13);

    I_v = reshape(eye(4), 16,1);
    % solve differential equation system
    [~,dxdtInt] = ode113(@PBCR4BP,[t0 tFin],[X(1:4); I_v],options, mu);

    % retrieve final state 
    XF = dxdtInt(end,:);

    % retrieve final STM and reshape as matrix
    STM = reshape(XF(5:20),[4,4]);

    % this is the simple shooting case without gradients
    DELTAV = abs(sqrt((u-y)^2+(v+x+mu)^2)-sqrt((1-mu)/r0))+...
        +abs(sqrt((XF(3)-XF(2))^2+(XF(4)+XF(1)+mu-1)^2)-sqrt(mu/rF));
    
    % case with gradient
        if nargout > 1
            % define all the derivatives of the cost function

            % dDV/d(initial time)
            v1N = sqrt((u-y)^2+(v+x+mu)^2);
            dynTi = PBCR4BP(t0,[x,y,u,v], mu);
            DV1dxi = (1/v1N)*[v+x+mu,  y-u, u-y, v+x+mu];

            %dDV/d(final time)
            v2N = sqrt((XF(3)-XF(2))^2+(XF(4)+XF(1)+mu-1)^2);
              dynTf = PBCR4BP(tFin,[XF(1),XF(2),XF(3),XF(4)], mu);
              DV2dxf = (1/v2N)*[(XF(4)+XF(1)+mu-1), ...
                         (XF(2)-XF(3)), ...
                         (XF(3)-XF(2)), ...
                         (XF(4)+XF(1)+mu-1)];

              DV1 = sqrt((u-y)^2+(v+x+mu)^2) - sqrt((1-mu)/r0);
              DV2 = sqrt((XF(3)-XF(2))^2+(XF(4)+XF(1)+mu-1)^2) - sqrt(mu/rF);

              grad =[(DV1/abs(DV1)*DV1dxi + DV2/abs(DV2)*DV2dxf*STM)';
                                      -DV2/abs(DV2)*DV2dxf*STM*dynTi;
                                           DV2/abs(DV2)*DV2dxf*dynTf];
        end
end

function [c, cEq, DC, DCeq] = distAndTan(X, r0,rF, mu)
% PURPOSE: define constraints matrices for fmincon
% INPUT: 
%     X, state vector (4x1), including x, y, u, v
%     r0, reduced initial parking orbit radius
%     rF, reduced final parking orbit radius
%     mu, mass parameter [scalar]
 c = [];
 
 % extract state
x = X(1);
y = X(2);
u = X(3);
v = X(4);
t0 = X(5);
tF = X(6);
options = odeset('Abstol', 1e-14,'Reltol', 1e-14);
I_v = reshape(eye(4), 16,1);

    % solve dynamics to retrieve final state
    [~,dxdtInt] = ode113(@PBCR4BP,[t0 tF],[x;y;u;v;I_v],options,mu);
    % store final state in variables and reshape
    xf = dxdtInt(end,:);
    STMf = reshape(xf(5:20),[4,4]);

   
 cEq = [(x+mu)^2+y^2-r0^2;
     (x+mu)*(u-y)+y*(v+x+mu);
     (xf(1)+mu-1)^2+xf(2)^2-rF^2;
     (xf(1)+mu-1)*(xf(3)-xf(2))+xf(2)*(xf(4)+ xf(1)+mu-1)];
 
 if nargout > 2
        % initialize gradient of constraint function
        DC = [];
        % get dynamics at initial and final time
        dynTi = PBCR4BP(t0,[x,y,u,v],mu);
        dynTf = PBCR4BP(tF,[xf(1),xf(2),xf(3),xf(4)], mu);
           
        dceq11 = [2*(x+mu),2*y,0,0]';
        dceq12 = [u,v,x+mu,y]';
        dceq13 = ([2*(xf(1)+mu-1),2*xf(2),0,0]*STMf)';
        dceq14 = ([xf(3),xf(4),(xf(1)+mu-1),xf(2)]*STMf)';
        dceq23 = 2*(xf(1)+mu-1)*(-STMf(1,:)*dynTi) + 2*xf(2)*(-STMf(2,:)*dynTi);
        dceq33 = 2*(xf(1)+mu-1)*xf(3) + 2*xf(2)*xf(4);
        dceq24 = xf(3)*(-STMf(1,:)*dynTi) + xf(4)*(-STMf(2,:)*dynTi) + (xf(1)+mu-1)*(-STMf(3,:)*dynTi) + xf(2)*(-STMf(4,:)*dynTi);
        dceq34 = xf(3)*dynTf(1) + xf(4)*dynTf(2) + (xf(1)+mu-1)*dynTf(3) + xf(2)*dynTf(4);

        DCeq = [dceq11,dceq12, dceq13,    dceq14;
                0,    0,     dceq23,  dceq24;
                0,    0,     dceq33,  dceq34];
  end
end

% MULTIPLE SHOOTING FUNCTIONS
function [DELTAV, g] = costFunMulti(X,r0,rF, mu)
    % PURPOSE: Define DELTAV function to be minimized
% INPUT: 
%     t, time vector [nx1]
%     X, state vector (4x1), including x, y, u, v
%     mu, mass parameter [scalar]

x = X(1);
y = X(2);
u = X(3);
v = X(4);

% this is the simple shooting case without gradients
DELTAV = abs(sqrt((u-y)^2+(v+x+mu)^2)-sqrt((1-mu)/r0))+...
    +abs(sqrt((X(15)-X(14))^2+(X(16)+X(13)+mu-1)^2)-sqrt(mu/rF));
    if nargout > 1
        % define all the derivatives of the cost function
        % wrt initial time
        v1N = sqrt((u-y)^2+(v+x+mu)^2);
        DV1dxi = (1/v1N)*[v+x+mu,  y-u, u-y, v+x+mu];
        %wrt final time
        v2N = sqrt((X(15)-X(14))^2+(X(16)+X(13)+mu-1)^2);
          DV2dxf = (1/v2N)*[(X(16)+X(13)+mu-1), ...
                     (X(14)-X(15)), ...
                     (X(15)-X(14)), ...
                     (X(16)+X(13)+mu-1)];

          DV1 = sqrt((u-y)^2+(v+x+mu)^2) - sqrt((1-mu)/r0);
          DV2 = sqrt((X(15)-X(14))^2+(X(16)+X(13)+mu-1)^2) - sqrt(mu/rF);

       g = [DV1/abs(DV1)*DV1dxi';...
           0;0;0;0;0;0;0;0;...
           DV2/abs(DV2)*DV2dxf';
           0;0];
    end
end
function [c, cEq, DC, DCeq] = multiDistAndTan(X, r0,rF, mu,N)
% PURPOSE: define constraints matrices for fmincon (multishoot case)
% INPUT: 
%     X, state vector (4x1), including x, y, u, v
%     r0, reduced initial parking orbit radius
%     rF, reduced final parking orbit radius
%     mu, mass parameter [scalar]
%  For multiple shooting case, given the states of initial guess, returns the
%  nonlinear constraints functions and its gradients as output.
% Inputs:
%    X   : [-]           The initial guess
% Constants:
%   mu   : [-]           Gravitational parameter of Sun
%   RE   : [DU]          Mean Earth's radius
%   RM   : [DU]          Mean Moon's radius
%   ri   : [DU]          Departure radius
%   rf   : [DU]          Arrival radius
%   rho  : [-]           Scaled Sun-(Earth+Moon) distance
% omega_s: [-]           Scaled angular velocity of the Sun
%   ms   : [-]           Scaled mass of the Sun
% options: [-]           Integrator options
% Outputs:
%   c    : [-]           Nonlinear inequality constraints functions values
%  ceq   : [-]           Nonlinear equality constraints functions values
%  DC    : [-]           Nonlinear inequality constraints functions gradients
%  DCeq  : [-]           Nonlinear equality constraints functions gradients

% constants
DU = 3.84405e5;
rEarth = 6378/DU; % [km] - Earth radius
rMoon = 1738/DU; % [km]  - Moon radius

c = [];
cEq = zeros(16,1);
 
options = odeset('Abstol', 1e-14,'Reltol', 1e-14);

 tj = []; tj(1) = X(17);
 I_v = reshape(eye(4), 16,1);
 
 guessXode = [];
 guessXode(:,1) = [X(1:4);I_v];
 guessXode(:,2) = [X(5:8);I_v];
 guessXode(:,3) = [X(9:12);I_v];
 
 solX = [];
 solX(:,1) = guessXode(:,1);
 
 % perform integration at each step of time grid

    for j=2:N
        tj(j)=tj(1)+(j-1)/(N-1)*(X(18)-tj(1));
        tspan = [tj(j-1),tj(j)];
        [~, xx] = ode113(@PBCR4BP, tspan, guessXode(:,j-1),options,mu);
        solX(:,j) = xx(end,:);
        cEq(4*(j-2)+(1:4)) = xx(end,1:4)'-X(4*(j-1)+(1:4));

        if nargout > 2
            fIn(1:4,j-1) = PBCR4BP(tj(j-1),[xx(1,1),xx(1,2),xx(1,3),xx(1,4)], mu);
            fF(1:4,j-1) = PBCR4BP(tj(j),[xx(end,1),xx(end,2),xx(end,3),xx(end,4)], mu);
        end
    end
% NL inequality constraints
c = [rEarth^2 - (X(1)+mu)^2 - X(2)^2;
     rMoon^2 - (X(1)+mu-1)^2 - X(2)^2;
     rEarth^2 - (X(5)+mu)^2 - X(6)^2;
     rMoon^2 - (X(5)+mu-1)^2 - X(6)^2;
     rEarth^2 - (X(9)+mu)^2 - X(10)^2;
     rMoon^2 - (X(9)+mu-1)^2 - X(10)^2;
     rEarth^2 - (X(13)+mu)^2 - X(14)^2;
     rMoon^2 - (X(13)+mu-1)^2 - X(14)^2;
     X(17)-X(18)];

% NL equality constraints
cEq(13:16) = [(X(1)+mu)^2 + X(2)^2 - r0^2;
              (X(1)+mu)*(X(3)-X(2)) + X(2)*(X(4)+X(1)+mu);
              (X(13)+mu-1)^2 + X(14)^2 - rF^2;
              (X(13)+mu-1)*(X(15)-X(14)) + X(14)*(X(16)+X(13)+mu-1)];
    if nargout > 2
        % Gradients of constraints
        % initial time
        c11 = zeros(18,1);
        c12 = zeros(18,1);
        c11(1) = -2*(X(1)+mu);
        c11(2) = -2*X(2);
        c12(1) = -2*(X(1)+mu-1);
        c12(2) = -2*X(2);
        % step 2
        c21 = zeros(18,1);
        c22 = zeros(18,1);
        c21(5) = -2*(X(5)+mu);
        c21(6) = -2*X(6);
        c22(5) = -2*(X(5)+mu-1);
        c22(6) = -2*X(6);
        % step 3
        c31 = zeros(18,1);
        c32 = zeros(18,1);
        c31(9) = -2*(X(9)+mu);
        c31(10) = -2*X(10);
        c32(9) = -2*(X(9)+mu-1);
        c32(10) = -2*X(10);
        % step 4
        c41 = zeros(18,1);
        c42 = zeros(18,1);
        c41(13) = -2*(X(13)+mu);
        c41(14) = -2*X(14);
        c42(13) = -2*(X(13)+mu-1);
        c42(14) = -2*X(14);
        % final time
        cfinal  = zeros(18,1);
        cfinal(17) = 1;
        cfinal(18) = -1;

        % gradient for equality constraints
        % initial time
        dceq11 = zeros(18,1);
        dceq11(1:4) = [2*(X(1)+mu);2*X(2);0;0];
        dceq12 = zeros(18,1);
        dceq12(1:4) = [X(3);X(4);X(1)+mu;X(2)];
        % step2
        STM2 = reshape(solX(5:20,2),[4,4]);
        STM2 = transpose(STM2);
        dceq2 = zeros(18,4);
        dceq2(1:4,1:4) = STM2;
        dceq2(5:8,1:4) = -eye(4);
        dceq2(17,1:4) = (-fIn(1:4,1)'*STM2(:,:))*(N-2+1)/(N-1) + fF(1:4,1)'*(N-2)/(N-1);
        dceq2(18,1:4) = (-fIn(1:4,1)'*STM2(:,:))*(1-1)/(N-1) + fF(1:4,1)'*1/(N-1);
        % step 3
        STM3 = reshape(solX(5:20,3),[4,4]);
        STM3 = transpose(STM3); 
        dceq3 = zeros(18,4);
        dceq3(5:8,1:4) = STM3;
        dceq3(9:12,1:4) = -eye(4);
        dceq3(17,1:4) = (-fIn(1:4,2)'*STM3(:,:))*(N-3+1)/(N-1) + fF(1:4,2)'*(N-3)/(N-1);
        dceq3(18,1:4) = (-fIn(1:4,2)'*STM3(:,:))*(2-1)/(N-1) + fF(1:4,2)'*2/(N-1);
        % step 4
        STM4 = reshape(solX(5:20,4),[4,4]);
        STM4 = transpose(STM4);
        dceq4 = zeros(18,4);
        dceq4(9:12,1:4) = STM4;
        dceq4(13:16,1:4) = -eye(4);
        dceq4(17,1:4) = (-fIn(1:4,3)'*STM4(:,:))*(N-4+1)/(N-1) + fF(1:4,3)'*(N-4)/(N-1);
        dceq4(18,1:4) = (-fIn(1:4,3)'*STM4(:,:))*(3-1)/(N-1) + fF(1:4,3)'*3/(N-1);
        % final time
        dceqf1 = zeros(18,1);
        dceqf1(13:16) = [2*(X(13)+mu-1);2*X(14);0;0];
        dceqf2 = zeros(18,1);
        dceqf2(13:16) = [X(15);X(16);(X(13)+mu-1);X(14)];
        
        % assemble gradients
        DC = [c11,c12,c21,c22,c31,c32,c41,c42,cfinal];
        DCeq = [dceq2,dceq3,dceq4,dceq11,dceq12,dceqf1,dceqf2];
    end
end

% REFERENCE FRAME TRANSFORMATION
function xTrans = inertialTrans(X,time, mu, body)
% PURPOSE: Transform reference frame from synodic to ECI
% INPUT: 
%     X, state vector (4x1), including x, y, u, v
%     time, time vector [nx1]
%     mu, mass parameter [scalar]
%     body, 'string' operator to identify the body of the r.f. orign.
%     
% OUTPUT
%     xTrans, transformed state vector

    t = length(time);
    if body == 'M'
       MU = mu;
    elseif body == 'm'
       MU = mu-1; 
    end
     xTrans = zeros(size(X));
    for i = 1:t  
       x = X(i,1);
       y = X(i,2);
       u = X(i,3);
       v = X(i,4);
       xTrans(i,1) = (x+MU)*cos(time(i))-y*sin(time(i));
       xTrans(i,2) = (x+MU)*sin(time(i))+y*cos(time(i));
       xTrans(i,3) = (u-y)*cos(time(i))-(v+x+MU)*sin(time(i));
       xTrans(i,4) = (u-y)*sin(time(i))+(v+x+MU)*cos(time(i));
    end
end
function circle(x,y,r,color)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'--','Color',color);
end