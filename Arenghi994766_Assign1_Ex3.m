%% SPACECRAFT GUIDANCE AND NAVIGATION
% PROFESSORS: FRANCESCO TOPPUTO, PIERLUIGI DI LIZIA
% A.Y. 2022-23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%   AUTHOR: ANDREA ARENGHI   %
%   PERSONAL CODE: 10523901  %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ASSIGNMENT 1: EX3
% prepare environment
clc
clear
close all
warning('off','all')

%% TASK 1

% Load ephemeris from SPICE
cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\de425s.bsp');
cspice_furnsh('kernels\gm_de431.tpc');
cspice_furnsh('kernels\pck00010.tpc');

% Extract parameters
muSun = cspice_bodvrd('SUN','GM',1); %[km^3/s^2]

% Define scaled factors accordin to AU (distance) and muSun
AU = 149597870.691; % Astronomical Unit in [km]
VEL = sqrt(muSun/ AU); %[km/s]
TIME = AU/VEL;
MASS = 1500;

% Normalize/adimensionalize constants
m0 = 1500/MASS; % [-] - s/c mass
Tmax = (150e-6)*TIME*TIME/MASS/AU; %[-] - max thrust
g0 = (9.81e-3)*TIME*TIME/AU; % [-] - strandard grav. acceleration 
Isp = 3000/TIME; %[-] - specific impulse
muSun = muSun*TIME*TIME/(AU^3);
N = 4; % [-] - thrusters number
% Unify parameters in struct variable to simplify functions' inputs
param.m0 = m0;
param.Tmax = Tmax;
param.g0 = g0;
param.Isp = Isp;
param.muSun = muSun;
param.AU = AU;
param.VEL = VEL;
param.TIME = TIME;
param.MASS = MASS;

% Initialize TOF and angles components' vectors in 4 and 3 thrusters case
dateVect = zeros(4,2);
INplane4th = [];
OUTplane4th = [];
INplane3th = [];
OUTplane3th = [];

% Open cases as nominal (4) or degraded (3) configuration applies
for thrusters = 4:-1:3
   % Update number of thrusters according to case
   param.N = thrusters;

    % Retrieve time of launch and Earth ephemeris
    tLaunch = cspice_str2et('2022-08-03-12:45:20.000 UTC');
    Earth0 = cspice_spkezr('Earth',tLaunch,'ECLIPJ2000','NONE','SUN'); % [km km/s]
    % Adimensionalize with respect to scale factors
    Earth0(1:3) = Earth0(1:3)/AU;
    Earth0(4:6) = Earth0(4:6)/VEL;

    tLaunch = tLaunch/TIME; % adimensional launch time
    
    % Provide initial guess for final time (adimensionalized): 1 year TOF
    tArrivalGuess = tLaunch + 2*pi;
    EXITFLAG = 0;
    % Find the lambda0 parameter solving the equation of the TBP.
    while EXITFLAG ~= 1
        lambdaGuess = randi([-10 10],7,1);
        lambdaGuess(7) = abs(lambdaGuess(7));
        % Make the resolution more quicker
        if thrusters == 4
            iter = 100;
        else
            iter = 200;
        end
        
        guess = [lambdaGuess; tArrivalGuess];
        options = optimoptions('fsolve','Display','none');
        % Solve optimization problem
        [STATE,FVAL,EXITFLAG,OUTPUT] = fsolve(@(x) optimalProb(x, tLaunch,  Earth0, param), guess, options);
    end
    
    %% POSTPRODUCTION
    % Print solution state
    fprintf('The costate solution for the %d thrusters case was: ', thrusters)
    STATE %#ok<NOPTS>

    % Compute arrival time and TOF
    tArrival = STATE(end);
    TOF = (tArrival-tLaunch)*TIME/3600/24; % [days]
    % store variables for TOF and arrivalDate
    dateVect(thrusters,:) = [thrusters, TOF];
    % Print arrival date
    arrivalDate = cspice_et2utc(tArrival*param.TIME,'C',4)%#ok<NOPTS>
    
    % Recover Mars state at arrival
    MarsArrival = cspice_spkezr('Mars',tArrival*TIME,'ECLIPJ2000','NONE','SUN')'; % [km km/s]
    % Adimensionalize Mars state
    MarsArrival(1:3) = MarsArrival(1:3)/param.AU;
    MarsArrival(4:6) = MarsArrival(4:6)/param.VEL;

    % Compute spacecraft orbit
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    launchState = [Earth0; param.m0; STATE(1:7)];
    [timeOut, stateFull]= ode113(@TBP_ODE, [tLaunch tArrival],launchState, options, param);
    
    % Compute Hamiltonian
    [Hvec,hMean] = hamiltonian(stateFull,param);
    
    % Compute in-plane and out-of-plane angles and store variables for
    % later plotting
    [INplane, OUTplane] = angles(stateFull);
    if thrusters == 4
        INplane4th(:,1) = (timeOut-ones(length(timeOut),1)*tLaunch)*param.TIME/3600/24;
        INplane4th(:,2) = INplane/pi*180;
        OUTplane4th(:,1) = (timeOut-ones(length(timeOut),1)*tLaunch)*param.TIME/3600/24;
        OUTplane4th(:,2) = OUTplane/pi*180;
    else
        INplane3th(:,1) = (timeOut-ones(length(timeOut),1)*tLaunch)*param.TIME/3600/24;
        INplane3th(:,2) = INplane/pi*180;
        OUTplane3th(:,1) = (timeOut-ones(length(timeOut),1)*tLaunch)*param.TIME/3600/24;
        OUTplane3th(:,2) = OUTplane/pi*180;
    end
    % Print consumed mass
    fprintf('\n The consumed propellant mass is m_p = %6.2f kg.',param.MASS - stateFull(end,7)*param.MASS)
    
    % Plot hamiltonian
    figure()
    plot((timeOut-ones(length(timeOut),1)*tLaunch)*param.TIME/3600/24,Hvec, 'Linewidth', 1.4);
    title('Hamiltonian functional', 'FontSize',16)
    subtitle([ num2str(thrusters) ' thrusters case'],'FontSize',14)
    xlabel('Time of flight [days]','FontSize',14)
    ylabel('H(t)','FontSize',14)
    ylim([hMean-1e-10 hMean+1e-10])
    yticks(linspace(hMean-1e-10,hMean+1e-10,3))
    ytickformat('%.10f');
    grid on
    hold off
    
    % Compute error in position and velocity at arrival wrt the Mars
    % final state
    errPos = norm(MarsArrival(end,1:3)-stateFull(end,1:3))*AU; % [km]
    errVel = norm(MarsArrival(end,4:6)-stateFull(end,4:6))*VEL*1000; % [m/s]
    fprintf('\nThe error in position at arrival is %.3f km. \n', errPos);
    fprintf('\nThe error in velocity at arrival is %.3f x 10^-3 m/s. \n', errVel*1000);
    
    % Plot Earth and Mars yearly revolution
    EarthRev = 366; % [days]
    MarsRev = 688; % [days]
    nPlot = 1000;
    tEOrb = linspace(0,EarthRev,nPlot);
    tMOrb = linspace(0,MarsRev,2*nPlot);
    EarthOrbit = zeros(nPlot,3);
    MarsOrbit = zeros(2*nPlot,3);
    rowE =1;
    rowM =1;
    for i = 1:nPlot
        instantE = tEOrb(i)*24*3600;
        ephemE = cspice_spkezr('Earth',instantE,'ECLIPJ2000','NONE','SUN'); % [km km/s]
        EarthOrbit(rowE,:) = ephemE(1:3)/AU;
        rowE = rowE+1;
    end

    for j = 1:2*nPlot
        instantM = tMOrb(j)*24*3600;
        ephemM = cspice_spkezr('Mars',instantM,'ECLIPJ2000','NONE','SUN'); % [km km/s]
        MarsOrbit(rowM,:) = ephemM(1:3)/AU;
        rowM = rowM+1;
    end

    % PLOT TRAJECTORIES
    figure()
    plot3(stateFull(:,1),stateFull(:,2),stateFull(:,3),'LineWidth',1.3); % spacecraft
    hold on
    grid on
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    plot3(Earth0(1),Earth0(2),Earth0(3),'ok','MarkerFaceColor','#009900','MarkerSize',8);
    plot3(EarthOrbit(:,1),EarthOrbit(:,2),EarthOrbit(:,3),'Color','#009900');
    plot3(0,0,0,'ok','MarkerSize',10,'MarkerFaceColor','#ffcc00')
    plot3(MarsArrival(1),MarsArrival(2),MarsArrival(3),'ok','MarkerFaceColor','#cc2900','MarkerSize',8);
    plot3(MarsOrbit(:,1),MarsOrbit(:,2),MarsOrbit(:,3),'Color','#cc2900');

    legend('Spacecraft','Earth at Launch','Earth Orbit','Sun','Mars at Arrival','Mars Orbit',...
        'Location','bestoutside',...
        'Orientation','horizontal',...
        'NumColumns',3,...
        'FontSize',12)
    title('Time optimal continuous guidance problem','FontSize',17)
    subt = ['TOF = ', int2str(TOF),' days. Thrusters: ',int2str(thrusters)];
    subtitle(subt,'FontSize',14);
end

% Plot angles comparisons
figure()
subplot(1,2,1)
plot(INplane4th(:,1),INplane4th(:,2),'--','Color','#29a329','Linewidth',1.5)
hold on
plot(INplane3th(:,1),INplane3th(:,2),'Linewidth',1.5)
legend('4 thrusters','3 thrusters','Fontsize',13)
subtitle('IN PLANE','Fontsize',13)
grid on
xlabel('TOF [days]','Fontsize',13)
ylabel('Angle [deg]','Fontsize',13)

subplot(1,2,2)
plot(OUTplane4th(:,1),OUTplane4th(:,2),'--','Color','#29a329','Linewidth',1.5);
hold on
plot(OUTplane3th(:,1),OUTplane3th(:,2),'Linewidth',1.5)
legend('4 thrusters','3 thrusters','Fontsize',13)
subtitle('OUT OF PLANE','Fontsize',13)
grid on
xlabel('TOF [days]','Fontsize',13)
ylabel('Angle [deg]','Fontsize',13)


%% FUNCTIONS
% Optimal problem resolution
function [BC] = optimalProb(costateOpt, tLaunch,  Earth0, param)
    %   PURPOSE: solve the optimal problem
    %   INPUT: -stateOpt, [8x1] include 7 components of the costate and guessed
    %                   time of arrival;
    %           - tLaunch, launch time [s] in ET
    %           - Earth0, state of Earth (hence of the spacecraft) at launch
    %           - param, includes several parameters (mass, AU, VEL, Tmax, N, Isp, g0..)
    %   OUPUT: - BC, [3x1] vector including boundary condition and
    %           transversality condition to be verified for optimal solution

    % Extract variables
    lam0 = costateOpt(1:7);
    tArrival = costateOpt(8);
    
    % Unify in guess vector to feed into TBP
    stateGuess = [Earth0; param.m0; lam0];
    
    % Set tolerances
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    
    % Integrate two body problem
    [~, stateOut] = ode113(@TBP_ODE, [tLaunch tArrival], stateGuess, options, param);


    % Mars state at final time
    MarsArrival = cspice_spkezr('Mars',tArrival*param.TIME,'ECLIPJ2000','NONE','SUN')'; % [km km/s]
    % adimensionalized Mars state
    MarsArrival(1:3) = MarsArrival(1:3)/param.AU;
    MarsArrival(4:6) = MarsArrival(4:6)/param.VEL;

    % retrieve final state to compute final conditions
    posArrival = stateOut(end,1:3);
    velArrival = stateOut(end,4:6);
    massArrival = stateOut(end,7);
    lambdaArrival = stateOut(end, 8:14);

    % where
    lRArr = lambdaArrival(1:3);

    lVArr = lambdaArrival(4:6);

    lMArr = lambdaArrival(7);
    
    % Compute switching function to set value of U (bang-bang control)
    SF = -norm(lVArr)*param.Isp*param.g0/massArrival-lMArr;
    if SF >= 0
        U = 0;
    else
        U = 1;
    end

    % compute Hamiltonian at arrival time
    HArr = 1+dot(lRArr,velArrival)-...
        param.muSun/(norm(posArrival)^3)*dot(posArrival,lVArr)+...
        param.N*param.Tmax/(param.g0*param.Isp)*U*SF;

    % compute function psi (velocity and acceleration) of Mars
    psiArr = [MarsArrival(4:6)'; -param.muSun/(norm(MarsArrival(1:3))^3)*MarsArrival(1:3)'];

    % BC
    BC = [ (stateOut(end,1:6)-MarsArrival(1:6))';
           lMArr;
           HArr- dot(lambdaArrival(1:6),psiArr)];
end

% Define ODE of the sun centered two body problem
function [dState] = TBP_ODE(t,state,param)
    % PURPOSE: define two body problem dynamics controlled with low thrust
    %          propulsion
    
    % INPUTS:
        %     -t, time vector
        %     -state, [14x1] [position (3), velocity (3), mass, lambdas (7)]
        %     -param, includes several parameters (mass, AU, VEL, Tmax, N, Isp, g0..)
    % OUTPUTS: derivatives of state vector dState
    
    % Extract constants for easier use
    mu = param.muSun;
    g0 = param.g0;
    N = param.N;
    Tmax = param.Tmax;
    Isp = param.Isp;
    
    % Define position and velocity vector
    X = state(1:3);
    
    % Extract state elements
    x  = state(1);
    y  = state(2);
    z  = state(3);
    u  = state(4);
    v  = state(5);
    w  = state(6);
    m = state(7);
    
    % Extract costate
    lambda = state(8:14);
    lR = lambda(1:3);
    lV = lambda(4:6);
    lM = lambda(7);
    
    % Define switching function
    SF = - norm(lV)*Isp*g0/m-lM;
    
    % Define switch in bang-bang control strategy
    if SF >= 0
        U = 0;
    else
        U = 1;
    end
    
    % compute postion vector norm
    r = sqrt(x^2+y^2+z^2);
    % Define derivatives vector
    dState = [
        u;
        v;
        w;
        -mu/(r^3)*x - N*U/m*Tmax*lV(1)/norm(lV);
        -mu/(r^3)*y - N*U/m*Tmax*lV(2)/norm(lV);
        -mu/(r^3)*z - N*U/m*Tmax*lV(3)/norm(lV);
        -N*U*Tmax/Isp/g0;
        -(3*mu/(r^5))*dot(X,lV)*x + mu/(r^3)*lV(1);
        -(3*mu/(r^5))*dot(X,lV)*y + mu/(r^3)*lV(2);
        -(3*mu/(r^5))*dot(X,lV)*z + mu/(r^3)*lV(3);
        -lR(1);
        -lR(2);
        -lR(3);
        - U*norm(lV)*Tmax*N/(m^2)];
end

% Compute Hamiltonian
function [H, hAve] = hamiltonian(state,param)
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
 % Compute angles
function [INplane, OUTplane] = angles(state)
    % PURPOSE: compute the Hamiltonian funcitonal for each instant of the TOF
    % INPUT: - state, [14x1] [position (3), velocity (3), mass, lambdas (7)]
   
    % OUTPUT: - INplane [nx1], inplane angle across TOF
    %         - OUTplane [nx1], out of plane angle across TOF
    
    L = max(size(state));
    INplane=  zeros(L,1);
    OUTplane = zeros(L,1);
    for i = 1:L
        lamV = state(i,11:13);
        % Define vector alpha
        alpha = - lamV/norm(lamV);
        INplane(i) =  atan2(alpha(2),alpha(1));
        OUTplane(i) =  asin(alpha(3));
    end
end

    %%% EXERCISE 3 END %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%   AUTHOR: ANDREA ARENGHI   %
%   PERSONAL CODE: 10523901  %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%