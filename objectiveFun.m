function P = objectiveFun(state,PIX)
    % This function shall allow the selection of the performance index (PIX) to
    % consider according to the input given
    %   Detailed explanation goes here

    % Extract state
    x = state(1);
    y = state(2);
    z = state(3);
    u = state(4);
    v = state(5);
    w = state(6);
    m = state(7);
    
    
    switch PIX
        case 'deltav'
            fprintf('The selected performance index is the Delta V.');

        case 'tof'
            fprintf('The selected performance index is the time of flight.');
    
    end
   


    P = state.*2;



end

