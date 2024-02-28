function [BC] = optimalProbDirect(costateOpt, tLaunch,  state0, param, targetState)
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
    stateOrigin = [state0; lam0];
    
    % Set tolerances
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    
    % Integrate two body problem
    [~, stateOut] = ode15s(@DYN_RAMSES_DIRECT, [tLaunch tArrival], stateOrigin, options, param);

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
    SF_fin = -norm(lVArr)*param.Isp*param.g0/massArrival-lMArr;
    if SF_fin >= 0
        U_fin = 0;
    else
        U_fin = 1;
    end
    
    %% TRASVERSALITY CONDITION??
% RHS_end = DYN_RAMSES_DIRECT(tArrival,stateOut,param);

statePsi = [targetState; state0(7); zeros(7,1)];
psiArr = DYN_RAMSES_DIRECT(tArrival,statePsi,param);

psiArr_small = psiArr(1:6);


[~,accFin] = DYN_RAMSES_DIRECT(tArrival,stateOut(end,:)',param);

H_tF = 1 + dot(lRArr,velArrival)+...
    dot(lVArr,accFin)+...
    param.Tmax/param.Isp/param.g0*U_fin*SF_fin;
    


TR_COND = H_tF-dot(lambdaArrival(1:6),psiArr_small);
    
% BC 
    BC = [ stateOut(end,1:6)'-targetState;
        lMArr;
%            lMArr+norm(lVArr)*param.Isp*param.g0/massArrival;
           TR_COND];

end