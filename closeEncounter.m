%% CLOSE ENCOUNTER
% Purpose: Plot the geometry of the close encounter including Earth, Moon
%          and Apophis.


% Date of creation: November 10th, 2023
% Author: Andrea Arenghi
% Context: Master Thesis: "CPO strategies at Apophis"

%% 
% Settings
clc
clear
close all

%% WORKFLOW
    % 1) Import SPICE ✓
    % 2) Define targets of interest (their code in catalogue)
    % 3) Identify period of interest ✓
    % 4) Identify suitable reference frame ✓
    % 5) Apply time/space frame conversion (✓)
    % 6) Recover ephemeris ✓
    % 7) Plot ephemeris ✓

%% ENVIRONMENT SETTING
% The loaded metakernel will include all the required kernels for the
% purpose of this thesis work. Please, be aware to first include the kernels in the
% proper folder 'kernels', then to add it to the list in the metakenrnel.

cspice_furnsh('MasterThesis.tm');

% The ephemeris for Apophis (99942)asteroid were generated from JPL Online Catalogue 
% with the following settings:
%     - Type: Small-Body SPK File
%     - Coordinate Center: Solar System Barycenter SSB [500@0]
%     - From January 1st, 2025 to December 31st, 2035
    
    
%% BUILD PARAMETER STRUCT
G = 6.67430e-20;        % Universal Gravity Constant, [km3/kg/s2]

muSun = cspice_bodvrd('SUN','GM',1); %[km^3/s^2]
muEarth = cspice_bodvrd('EARTH','GM',1); %[km^3/s^2]
muMoon = cspice_bodvrd('MOON','GM',1); %[km^3/s^2]

massApo = 4.3e10;      % Apophis' mass, [kg] from 'Heliotropic orbits...'
muApo = massApo*G;      % Apophis gravitational constant

% SRP parameters
P0 = 1367; % [W/m2]
c = 2.998e8; % speed of light, [m/s] NOTE: leave in meters/s otherwise SRP wrong
AU = 1.495978707e8; % Astronomical unit in [km]
Cr = 0.4; % reflectivity coefficient
A = 1; % equivalent surface of satellite, [m2]
massSC = 1000; % s/c mass, [kg]
AMratio = 0.015;
Psr = 4.6e-6; % N/m2

% Build struct for parameters
param.muSun = muSun;
param.muEarth = muEarth;
param.muMoon = muMoon;
param.muApo = muApo;
param.G = G;
param.P0 = P0;
param.Psr = Psr;
param.c = c;
param.AU = AU;
param.Cr = Cr;
param.A = A;
param.massSC = massSC;
param.AMratio = AMratio;
param.J2 = 'yes';
param.J2 = 'no';

% BUILD TIME VECTOR
% Define time window as UTC string
% An interval of +/- 3months from closest encounter are chosen for the
% analysis

steps = 3000;
t0s = '2029-01-01-00:00:00.000 UTC';
tFs = '2029-07-31-23:59:59.000 UTC';
% t0s = '2029-04-01-00:00:00.000 UTC';
% tFs = '2029-04-30-23:59:59.000 UTC';
t0 = cspice_str2et(t0s); %[s]
tF = cspice_str2et(tFs); %[s]
tVect = linspace(t0, tF, steps);
tVectNet = (tVect-ones(steps,1)'*tVect(1))/3600/24;

% Closest approach will be on April 13th
% Set moon time frame
t0Moon = cspice_str2et('2029-04-13-00:00:00.000 UTC'); %[s]
tFMoon = cspice_str2et('2029-04-13-23:59:00.000 UTC'); %[s]
tMoon = linspace(t0Moon,tFMoon,500);

% Observer: 10, the Sun
ID_Sun = '10';
apophis = cspice_spkezr('20099942',tVect,'ECLIPJ2000','NONE',ID_Sun);
earth = cspice_spkezr('Earth',tVect,'ECLIPJ2000','NONE',ID_Sun);
moon = cspice_spkezr('Moon',tVect,'ECLIPJ2000','NONE',ID_Sun);

distApophisSun = zeros(1,steps);
distApophisEarth = zeros(1,steps);
distApophisMoon = zeros(1,steps);
dE = earth-apophis;
dM = moon-apophis;
distance = zeros(1,steps);
for j = 1:steps
    distApophisSun(j) = norm(apophis(1:3,j));
    distApophisEarth(j)= norm(dE(1:3,j)) ;
    distApophisMoon(j) =  norm(dM(1:3,j)) ;
end
[minDist, indMin] = min(distApophisEarth);
% Identify epoch of closest encounter between Apophis and Earth
closestEncounterET = tVect(indMin); % ET


tVectRel = (tVect-ones(steps,1)'*tVect(indMin))/3600/24;



%% Perform RHS computation
options = odeset('reltol', 1e-13, 'abstol', 1e-13);

state0 = [50; 50; 10; 4; 4; 1];
[timeOut, fullState]= ode15s(@closeProxODE, tVect,state0, options,param);
    

%% COMPUTE ACCELERATIONS
% Retrieve data for plots:
%     A: acc(time) @ discrete ranges from Apophis
%     B: acc(range) @ discrete time instants/epochs

%% Type A
% Discretize ranges
rangeApophis = [5 15 500]; % [km]
nR = length(rangeApophis);
% NOTE: 5 shall be changed if more terms (J2, ellipsoid..will be included)
% accTypeA = zeros(nR,steps,5);
accTypeA = zeros(nR,steps,6);
% span vector of ranges
for rr = 1:nR
    
    % generate random position vector with given norm
    RADIUS = rangeApophis(rr);
    randoms = 2*rand(3,1)-1;
    randomState = (RADIUS/norm(randoms))*randoms;
    % initialize state (velocity is not relevant in this case)
    state = [randomState(1);randomState(2);randomState(3); 0; 0; 0];
    
    %     % Plot position of s/c if needed
    %     [X,Y,Z] = sphere();
    %     figure()
    %     quiver3(0,0,0,randomState(1),randomState(2),randomState(3))
    %     hold on
    %     axis equal
    %     surf(X*0.2,Y*0.2,Z*0.2)
    %     nameAxis(rangeApophis(rr))
    for jjj = 1:steps
       [~,accTypeA(rr,jjj,:)]= closeProxODE(tVect(jjj),state,param);
    end    
end
% Call function for plotting
plotAccelA(tVectRel,accTypeA,rangeApophis);


%% Type B
% Discretize times
epochs = [
    '2029-01-13-00:00:00.000 UTC';
    '2029-03-13-00:00:00.000 UTC';
    '2029-04-13-21:42:11.000 UTC';
    '2029-05-13-00:00:00.000 UTC';
    '2029-07-13-00:00:00.000 UTC'];

KP = min(size(epochs));

sizeRanges = 400;
ranges = linspace(2,30000, sizeRanges);
accTypeB = zeros(KP,sizeRanges,7);
for tt = 1:KP
    
    discreteEpoch = cspice_str2et(epochs(tt,:));
    
    for hh = 1:sizeRanges
        RADIUS = ranges(hh);
        randoms = 2*rand(3,1)-1;
        randomState = (RADIUS/norm(randoms))*randoms;
        % initialize state (velocity is not relevant in this case)
        state = [randomState(1);randomState(2);randomState(3); 0; 0; 0];
        [~,accTypeB(tt,hh,2:end)]= closeProxODE(discreteEpoch,state,param);
        accTypeB(tt,hh,1) = norm(state);
    end
end
plotAccelB(epochs,accTypeB);

%% PLOTS
% % % %% Plot relative distance between apophis and Sun and Earth
% % % figure()
% % % hdl = gcf;
% % % set(hdl, 'DefaultLineLineWidth', 2);
% % % % plot(tVectNet,distApophisSun/AU)
% % % grid on 
% % % hold on
% % % plot(tVectNet,distApophisEarth/AU)
% % % hold on
% % % xline(tVectNet(indMin),'--')
% % % xlabel('Days')
% % % ylabel('Distance [AU]')
% % % legend('Sun-Apophis','Earth-Apophis','Earth closest approach','Location','best')
% % % title('Apophis distance wrt Sun and Earth')
% % % startTime = cspice_et2utc(tVect(1),'C',0);
% % % sTime = startTime(1:end-9);
% % % finalTime = cspice_et2utc(tVect(end),'C',0);
% % % fTime = finalTime(1:end-9);
% % % subtitle(['Time window from ', sTime,' to ',fTime])
% % % axis tight

% % % %% Plot trajectories
% % % figure()
% % % hdl = gcf;
% % % set(hdl, 'DefaultLineLineWidth', 1.5);
% % % plot3(apophis(1,:),apophis(2,:),apophis(3,:),'Linewidth',1.3)
% % % hold on
% % % grid on
% % % plot3(earth(1,:),earth(2,:),earth(3,:),'Linewidth',1.1)
% % % plot3(moon(1,:),moon(2,:),moon(3,:),'Linewidth',1.1)
% % % plot3(0,0,0,'oy','MarkerFaceColor','#ffcc00','MarkerSize',9)
% % % axis equal
% % % plot3(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin),'dk','MarkerFaceColor','#3333cc','MarkerSize',6)
% % % scaleFact = 40000000;
% % % nameAxis(scaleFact);
% % % xlabel('x_{ECLIPJ2000} [km]')
% % % ylabel('y_{ECLIPJ2000} [km]')
% % % zlabel('z_{ECLIPJ2000} [km]')
% % % plot3(apophis(1,1),apophis(2,1),apophis(3,1),'ob','MarkerSize',2)
% % % plot3(earth(1,1),earth(2,1),earth(3,1),'or','MarkerSize',2)
% % % legend('Apophis','Earth','Moon','Sun','Closest approach','Location','best')
% % % title('Apophis closest approach')
% % % subtitle('Heliocentric orbits')
% % % 
% % % 
% % % %% Plot relative positions at closest encounter
% % % figure()
% % % plot3(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin),'ok','MarkerFaceColor','#3333cc','MarkerSize',6)
% % % hold on
% % % grid on
% % % axis equal
% % % plot3(earth(1,indMin),earth(2,indMin),earth(3,indMin),'ok','MarkerFaceColor','#33cc33','MarkerSize',15)
% % % plot3(moon(1,indMin),moon(2,indMin),moon(3,indMin),'ok','MarkerFaceColor','#999966 ','MarkerSize',10)
% % % % Plot velocities of the bodies
% % % quiver3(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin),...
% % %     1e8*apophis(4,indMin)/minDist,1e8*apophis(5,indMin)/minDist,1e8*apophis(6,indMin)/minDist,...
% % %     1,'Color','#3333cc','LineWidth',1.2)
% % % quiver3(earth(1,indMin),earth(2,indMin),earth(3,indMin),...
% % %     1e8*earth(4,indMin)/minDist,1e8*earth(5,indMin)/minDist,1e8*earth(6,indMin)/minDist,...
% % %     1,'Color','#33cc33','LineWidth',1.2)
% % % quiver3(moon(1,indMin),moon(2,indMin),moon(3,indMin),...
% % %     1e8*moon(4,indMin)/minDist,1e8*moon(5,indMin)/minDist,1e8*moon(6,indMin)/minDist,...
% % %     1,'Color','#999966','LineWidth',1.2)
% % % % To Sun vector
% % % quiver3(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin),...
% % %     -30*apophis(1,indMin)/minDist,-30*apophis(2,indMin)/minDist,-30*apophis(3,indMin)/minDist,...
% % %     1,'Color','#ff6600','LineWidth',1.2)
% % % scaleFact = 50000;
% % % quiver3(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin),scaleFact,0,0,'k')
% % % text(apophis(1,indMin)+scaleFact,apophis(2,indMin),apophis(3,indMin), 'x')
% % % quiver3(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin),0,scaleFact,0,'k')
% % % text(apophis(1,indMin),apophis(2,indMin)+scaleFact,apophis(3,indMin), 'y')
% % % quiver3(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin),0,0,scaleFact,'k')
% % % text(apophis(1,indMin),apophis(2,indMin),apophis(3,indMin)+scaleFact, 'z')
% % % legend('Apophis','Earth','Moon','v_{Apo}','v_E','v_M','To Sun','Location','best')
% % % title('Relative geometry at closest encounter')
% % % closeEnc = cspice_et2utc(closestEncounterET,'C',0);
% % % epochCE = closeEnc(1:end);
% % % subtitle(['Epoch: ',epochCE])
% % % xlabel('x [km]')
% % % ylabel('y [km]')
% % % zlabel('z [km]')

% % % 

%% WIP Plot close proximity trajectory
% figure(4)
% plot3(fullState(:,1),fullState(:,2),fullState(:,3))
% hold on
% grid on
% plot3(fullState(1,1),fullState(1,2),fullState(1,3),'*g','MarkerSize',5)
% plot3(fullState(end,1),fullState(end,2),fullState(end,3),'xr')
% %Plot initial velocity
% quiver3(fullState(1,1),fullState(1,2),fullState(1,3),...
%     fullState(1,4)*1e6,fullState(1,5)*1e6,fullState(1,6)*1e6,'-x')
% quiver3(0,0,0,-apophis(1,1)/2,-apophis(2,1)/2,-apophis(3,1)/2)
% quiver3(0,0,0,-apophis(1,indMin)/2,-apophis(2,indMin)/2,-apophis(3,indMin)/2)
% plot3(0,0,0,'ok','MarkerFaceColor','#3333cc','MarkerSize',3)
% scaleFactor = 6e6;
% quiver3(0,0,0,scaleFactor,0,0,'k')
% quiver3(0,0,0,0,scaleFactor,0,'k')
% quiver3(0,0,0,0,0,scaleFactor,'k')
% legend('Traj','Start','End','v_0','To Sun')
% axis equal
% xlabel('x_{ECLIPJ2000} [km]')
% ylabel('y_{ECLIPJ2000} [km]')
% zlabel('z_{ECLIPJ2000} [km]')

%% WIP: Plot perturbing accelerations 

% distances = linspace(1, 20, steps);
% for k= 1:step
%     accelerations(k,1) = distances(k);
%     
% end

% figure(5)
% hdl = gcf;
% set(hdl, 'DefaultLineLineWidth', 2);
% 
% semilogy(accelerations(:,1),accelerations(:,2));
% hold on
% grid on
% semilogy(accelerations(:,1),accelerations(:,3));
% semilogy(accelerations(:,1),accelerations(:,4));
% semilogy(accelerations(:,1),accelerations(:,5));
% semilogy(accelerations(:,1),accelerations(:,6));
% semilogy(accelerations(:,1),accelerations(:,7));
% % semilogy(accelerations(:,1),accelerations(:,7));
% ylabel('Accelerations [km/s^2]')
% xlabel('S/c-Apophis distance [AU]')
% legend('acc_{Apo}','acc_{Sun}','acc_{Earth}','acc_{Moon}','acc_{SRP}','acc_{J2}')
%% FUNCTIONS

function [dState,acc] = closeProxODE(t,state,param)
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
        
        
    totAcc = accApo+accSun+accEarth+accMoon+accSRP+accJ2;

    acc = [norm(accApo), norm(accSun), norm(accEarth), norm(accMoon), norm(accSRP), norm(accJ2)];
    
    dState = [
        u;
        v;
        w;
        totAcc(1);
        totAcc(2);
        totAcc(3)];
end

function fEncke = fEncke(q)
    fEncke = ((1+q)^1.5)-1;
end


function [] = nameAxis(scaleFact)

quiver3(0,0,0,scaleFact,0,0,'k')
text(scaleFact,0,0, 'x')
quiver3(0,0,0,0,scaleFact,0,'k')
text(0,scaleFact,0, 'y')
quiver3(0,0,0,0,0,scaleFact,'k')
text(0,0,scaleFact, 'z')
end

function [] = plotAccelA(t,m3D,ranges)
    [a,~,c] = size(m3D);
        
    for i = 1:a
      ACC(:,:) = m3D(i,:,:);
      figure()
      for j = 1:c
         h5 = gcf;
         set(h5, 'DefaultLineLineWidth', 2);
         semilogy(t,ACC(:,j)*1e3)
         hold on
         grid on
         axis tight
      end
      
    
    title('Accelerations')
    subtitle(['r = ',num2str(ranges(i)),' km'])
    legend('a_{Apo}','a_{Sun}','a_{Earth}','a_{Moon}','a_{SRP}','a_{J2}')
    ylabel('a(t) [m/s^2]')
    xlabel('Days from Closest Encounter')
    hold off
    end
    
    grid on
end

function [] = plotAccelB(epoch,accTypeB)
    a = min(size(epoch));
    for p = 1:a
        t = cspice_str2et(epoch(p,:));
        tPlot = cspice_et2utc(t,'C',0);
       M(:,:) = accTypeB(p,:,:);
        
        figure()
        h1 = gcf;
         set(h1, 'DefaultLineLineWidth', 1.4);
        semilogy(M(:,1),M(:,2)*1e3)
        hold on
        grid on
        axis tight
        semilogy(M(:,1),M(:,3)*1e3)
        semilogy(M(:,1),M(:,4)*1e3)
        semilogy(M(:,1),M(:,5)*1e3)
        semilogy(M(:,1),M(:,6)*1e3)
        semilogy(M(:,1),M(:,7)*1e3)
        xlabel('S/c - Apophis distance [km]')
        ylabel('Acceleration [m/s^2]')
        legend('a_{Apo}','a_{Sun}','a_{Earth}','a_{Moon}','a_{SRP}','a_{J2}')
        title('a(Range)')
        subtitle(['Epoch = ',tPlot(1:end-9)])
    end
end