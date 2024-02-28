%% testing snippet for indirect method
clc
clear
close all
warning('off','all')
cspice_furnsh('MasterThesis.tm');

m0 = 1000; % initial s/c mass, [kg]
mF = 900; % estimate of final mass

param = parameters(1);

  t0s = '2029-04-13-00:00:00.000 UTC';
  tFs = '2029-04-15-23:59:59.000 UTC';
  
  t0 = cspice_str2et(t0s); %[s]
  tF = cspice_str2et(tFs); %[s]

%% Problem solver

state0_SC = [5;5;5; ... %s/c position [km]
            3e-3;3e-3;1e-4;...
            m0]; % sé/c velocity [km/s]
        
%         state0_SC = [50;15;15; ... %s/c position [km]
%             1;1;1;...
%             m0]; % sé/c velocity [km/s]
tLaunch = t0;
tArrivalGuess = t0+10*3600;

targetState = [5;-5;-5;...
                0;1e-3;-1e-3]; % sé/c velocity [km/s]

EXITFLAG = 0;
    % Find the lambda0 parameter solving the equation of the TBP.
    while EXITFLAG ~= 1
        lambdaGuess = randi([-10 10],7,1);
        lambdaGuess(7) = abs(lambdaGuess(7));
        
        guess = [lambdaGuess; tArrivalGuess];
        options = optimoptions('fsolve','Display','none');
        % Solve optimization problem
        [STATE,FVAL,EXITFLAG,OUTPUT] = fsolve(@(x) optimalProbDirect(x, tLaunch,  state0_SC, param, targetState), guess, options);
    end